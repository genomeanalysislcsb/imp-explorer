"""imp-explorer fraction: add a second fraction of the same samples to a release, e.g. the extracellular DNA / RNA
next to the intracellular IMP3 runs that the release was exported from.

  imp-explorer fraction --data release_v0.3 --imp-root /path/imp_ext_mg --suffix ext
  imp-explorer fraction --data release_v0.3 --imp-root /path/imp_ext_mt --suffix ext --remap-only

Per sample folder (named like the release's sample ids, translated with --id-map / --id-regex / --id-case):

  IMP3 run of the fraction (Stats/, Analysis/, Binning/)   -> omics <omic>_<suffix>: mg_ext (read QC, Kraken2 /
      Bracken and mOTUs taxa, MAGs and their abundance, Mantis functions), as for the main export
  Analysis/taxonomy/bracken_<db>/<omic>/...  (--kraken-db)  -> taxa of <omic>_<suffix> from a Kraken2 rerun, also
      with --remap-only / without an IMP3 run (e.g. the single-end extracellular mt)
  remapping/annotation/<omic>.mantis.<db>.counts.tsv        -> omics <omic>_<suffix>_remap: the fraction's reads
      counted on the annotation of the main (e.g. intracellular) assembly (featureCounts tables, as IMP3 writes them)

Rows of these omics already in the release are replaced, so the command can be run again. The samples table gets
has_<omic>_<suffix> columns; manifest['fractions'] records what was added. Afterwards run 'imp-explorer qc' (the new
omics become QC layers) and 'imp-explorer cohort' (layers such as kraken.genus.mg_ext or kegg_ko.mt_ext_remap).
"""
import argparse
import concurrent.futures as cf
import json
import multiprocessing as mp
import os
import sys

import re

import pandas as pd

from .export import derive_ids, read_id_map
from .readers import imp3
from .release import read_manifest, read_table, write_table

TABLES_ = ("qc", "taxonomy", "bins", "bin_abundance", "functions")


def say(msg):
    sys.stderr.write("[fraction] %s\n" % msg)


def discover(root):
    out = []
    for d in sorted(os.listdir(root)):
        p = os.path.join(root, d)
        if os.path.isdir(p) and any(os.path.isdir(os.path.join(p, x)) for x in ("Stats", "status", "remapping")):
            out.append((d, p))
    return out


def remap_functions(base, sample, dbs, suffix):
    rows = []
    for o in imp3.OMICS:
        for db in dbs:
            p = imp3.find(base, "remapping/annotation/%s.mantis.%s.counts.tsv" % (o, db))
            fc = imp3._featurecounts(p) if p else {}
            rx = imp3.ID_PATTERN.get(db)
            if rx is not None:
                fc = {re.sub(r"\.\d+$", "", f): v for f, v in fc.items()
                      if rx.match(re.sub(r"\.\d+$", "", f))} or fc
            for f, v in fc.items():
                rows.append((sample, "%s_%s_remap" % (o, suffix), db, f, v))
    return pd.DataFrame(rows, columns=["sample_id", "omic", "db", "feature", "reads"])


def rename(res, sample, suffix):
    """omics of a fraction's own IMP3 run -> <omic>_<suffix>; bin ids <sample>:<suffix>:<bin>"""
    ren = lambda o: "%s_%s" % (o, suffix)
    out = {}
    for k in TABLES_:
        d = res[k].copy()
        if "omic" in d.columns:
            d["omic"] = d["omic"].astype(str).map(ren)
        if "bin_id" in d.columns:
            d["bin_id"] = d["bin_id"].astype(str).str.replace("%s:" % sample, "%s:%s:" % (sample, suffix), n=1)
        out[k] = d
    return out


def _one(args):
    sid, base, dbs, suffix, remap_only, kdb, ctax = args
    try:
        out = {k: pd.DataFrame() for k in TABLES_}
        omics = []
        if not remap_only and os.path.isdir(os.path.join(base, "Stats")):
            res = imp3.read_sample(base, sid, dbs, kdb, ctax)
            out = rename(res, sid, suffix)
            omics = ["%s_%s" % (o, suffix) for o in res["_omics"]]
            for n in res.get("_notes", []):
                say("  %s: %s" % (sid, n))
        elif kdb:                                  # no IMP3 run (e.g. single-end ext mt): only the Kraken2 rerun
            tx = imp3.taxonomy(base, sid, kraken_db=kdb)
            if len(tx):
                out["taxonomy"] = rename({**{k: pd.DataFrame() for k in TABLES_}, "taxonomy": tx}, sid, suffix)["taxonomy"]
                omics += ["%s_%s" % (o, suffix) for o in sorted(set(tx["omic"]))]
        rf = remap_functions(base, sid, dbs, suffix)
        if len(rf):
            out["functions"] = pd.concat([out["functions"], rf], ignore_index=True) if len(out["functions"]) else rf
            omics += sorted(set(rf["omic"]))
        return sid, out, omics, None
    except Exception as e:
        return sid, None, [], "%s: %s" % (type(e).__name__, e)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="imp-explorer fraction", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True, help="release directory (imp-explorer export)")
    ap.add_argument("--imp-root", required=True, help="folder with one sub-folder per sample of the fraction")
    ap.add_argument("--suffix", default="ext", help="added to the omic names: mg -> mg_<suffix>  [ext]")
    ap.add_argument("--label", default="", help="name of the fraction in the app  [extracellular for 'ext']")
    ap.add_argument("--remap-only", action="store_true", help="only the remapping/ counts (no IMP3 run of the "
                    "fraction, or not wanted)")
    ap.add_argument("--id-map", help="TSV/CSV: folder name, sample id")
    ap.add_argument("--id-regex", help="regular expression; its first group (or match) is the sample id")
    ap.add_argument("--id-case", choices=("keep", "upper", "lower"), default="keep")
    ap.add_argument("--function-dbs", default="", help="Mantis databases  [those of the release]")
    ap.add_argument("--kraken-db", default=None, help="Kraken2 database of a rerun (as for export)  [the release's "
                    "kraken_db]")
    ap.add_argument("--contig-taxonomy", default="", help="mmseqs2 LCA table of the fraction's contigs (as for "
                    "export), e.g. '/path/imp_ext_mg/contigs/{folder}ext/{folder}ext.{omic}.GTDB_lca.tsv'")
    ap.add_argument("--format", choices=("auto", "parquet", "tsv"), default="auto")
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args(argv)
    man = read_manifest(a.data)
    if not man:
        sys.exit("fraction: %s is not a release (no manifest.json)" % a.data)
    dbs = tuple(x for x in (a.function_dbs or ",".join(man.get("sources", {}).get("function_dbs") or
                                                       imp3.DEFAULT_FUNCTION_DBS)).split(",") if x)
    found = discover(a.imp_root)
    if not found:
        sys.exit("fraction: no sample folders with Stats/ or remapping/ in %s" % a.imp_root)
    try:
        ids = derive_ids([n for n, _ in found], read_id_map(a.id_map) if a.id_map else None, a.id_regex, a.id_case)
    except ValueError as e:
        sys.exit("fraction: %s" % e)
    smp = read_table(a.data, "samples")
    smp["sample_id"] = smp["sample_id"].astype(str)
    rel_ids = set(smp["sample_id"])
    kdb = a.kraken_db if a.kraken_db is not None else man.get("sources", {}).get("kraken_db", "")
    kdb = None if kdb in ("", "IMP3 output") else kdb
    jobs = [(sid, base, dbs, a.suffix, a.remap_only, kdb, a.contig_taxonomy or None) for sid, (_, base) in zip(ids, found) if sid in rel_ids]
    unknown = sorted(set(ids) - rel_ids)
    if unknown:
        say("%d folders are not samples of the release (left out; check --id-*): %s"
            % (len(unknown), ", ".join(unknown[:10])))
    if not jobs:
        sys.exit("fraction: no folder matches a sample of the release")
    parts = {k: [] for k in TABLES_}
    omics, errors = {}, {}
    with cf.ProcessPoolExecutor(max_workers=max(1, a.threads), mp_context=mp.get_context("spawn")) as ex:
        for i, (sid, res, om, err) in enumerate(ex.map(_one, jobs), 1):
            if err:
                errors[sid] = err
                say("%s: FAILED %s" % (sid, err))
                continue
            for k in TABLES_:
                if len(res[k]):
                    parts[k].append(res[k])
            omics[sid] = om
            say("%d/%d %s: %s" % (i, len(jobs), sid, ", ".join("%s %d" % (k, len(res[k])) for k in TABLES_)))
    new = sorted({o for v in omics.values() for o in v})
    if not new:
        sys.exit("fraction: nothing read (no IMP3 outputs or remapping counts found)")
    for k in TABLES_:
        old = read_table(a.data, k)
        if "omic" in old.columns:
            old = old[~old["omic"].astype(str).isin(new)]
        if k == "bins" and len(old):                    # bins have no omic: those of this fraction by their id
            old = old[~old["bin_id"].astype(str).str.contains(":%s:" % a.suffix, regex=False)]
        add = pd.concat(parts[k], ignore_index=True) if parts[k] else None
        df = pd.concat([old, add], ignore_index=True) if add is not None and len(old) else (add if add is not None
                                                                                         else old)
        write_table(df, a.data, k, a.format)
    for o in new:
        smp["has_" + o] = smp["sample_id"].map(lambda s: o in omics.get(s, []))
    write_table(smp, a.data, "samples", a.format)
    label = a.label or {"ext": "extracellular", "int": "intracellular"}.get(a.suffix, a.suffix)
    m = read_manifest(a.data)
    fr = m.setdefault("fractions", {})
    for o in new:
        base = o.split("_")[0]
        fr[o] = {"label": "%s %s%s" % (label, {"mg": "metagenome", "mt": "metatranscriptome"}.get(base, base),
                                       " (reads on the main assembly)" if o.endswith("_remap") else ""),
                 "source": os.path.abspath(a.imp_root), "samples": sum(1 for v in omics.values() if o in v),
                 "tables": [k for k in TABLES_ if any(len(p) and "omic" in p.columns and (p["omic"] == o).any()
                                                      for p in parts[k])],
                 "kraken_db": kdb or "IMP3 output",
                 "contig_taxonomy": a.contig_taxonomy or ""}
    m.setdefault("failed_fraction_samples", {}).update(errors)
    with open(os.path.join(a.data, "manifest.json"), "w") as fh:
        json.dump(m, fh, indent=2, default=str)
    for o in new:
        say("%s: %d samples (%s)" % (o, fr[o]["samples"], fr[o]["label"]))
    return 0 if not errors else 2


if __name__ == "__main__":
    sys.exit(main())
