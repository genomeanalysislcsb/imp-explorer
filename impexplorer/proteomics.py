"""imp-explorer proteomics: add metaproteomics results to a release (optional layer).

Input: a folder with generic tables, each as .parquet, .tsv or .tsv.gz
  proteomics_taxonomy    sample_id, fraction, rank, taxon, auc, rel, n_proteins          (required)
  proteomics_functions   sample_id, fraction, db, feature, auc, rel, n_proteins          (required)
  proteomics_samples     sample_id, fraction, n_proteins, n_peptides, total_auc          (optional)
  proteomics_proteins    protein_id, taxon, rank, lineage, product, kegg_ko, ...         (optional, catalogue)
  proteomics_abundance   sample_id, fraction, protein_id, auc, n_peptides                (optional)
rel = share of the total protein abundance (AUC) of the sample and fraction. A fraction is a protein extract, e.g.
extra (extracellular) and intra (pellet); each becomes an omic 'mp_<fraction>' next to mg and mt. The Expobiome
export (proteomics_export.py) writes this format.

Writes into the release: prot_taxonomy, prot_functions, prot_samples, prot_proteins, prot_abundance (sample ids mapped
with --id-map / --id-regex / --id-case) and manifest['proteomics']. The app shows the protein profiles next to mg and
mt (Taxonomy, Functions, Search) and on the Proteomics page; 'imp-explorer cohort' analyses them as layers
prot.<rank or database>.<fraction>."""
import argparse
import json
import os
import sys

import pandas as pd

from .release import read_manifest, read_table, write_table

TABLES = {"taxonomy": ("sample_id", "fraction", "rank", "taxon", "auc"),
          "functions": ("sample_id", "fraction", "db", "feature", "auc"),
          "samples": ("sample_id", "fraction"), "proteins": ("protein_id",),
          "abundance": ("sample_id", "fraction", "protein_id", "auc")}
REQUIRED = ("taxonomy", "functions")
FRACTION_LABEL = {"extra": "extracellular", "intra": "intracellular (pellet)"}
RANK_ALIAS = {"superkingdom": "domain", "kingdom": "domain"}


def say(msg):
    sys.stderr.write("[proteomics] %s\n" % msg)


def read_input(folder, name):
    for ext in (".parquet", ".tsv.gz", ".tsv"):
        p = os.path.join(folder, "proteomics_" + name + ext)
        if os.path.exists(p):
            return pd.read_parquet(p) if ext == ".parquet" else pd.read_csv(p, sep="\t", low_memory=False,
                                                                              keep_default_na=False, na_values=[""])
    return None


def omic(fraction):
    return "mp_" + str(fraction).lower()


def omic_label(o, manifest=None):
    lab = (manifest or {}).get("proteomics", {}).get("fractions", {}).get(o)
    if lab:
        return lab
    f = str(o)[3:] if str(o).startswith("mp_") else str(o)
    return "metaproteome, " + FRACTION_LABEL.get(f, f)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="imp-explorer proteomics", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True, help="release directory (imp-explorer export)")
    ap.add_argument("--prot", required=True, help="folder with proteomics_taxonomy / _functions / ...")
    ap.add_argument("--no-proteins", action="store_true", help="leave out the protein catalogue and abundances "
                    "(smaller release; no taxon x function view on the Proteomics page)")
    ap.add_argument("--id-map", help="TSV/CSV: sample id of the proteomics tables, release sample id")
    ap.add_argument("--id-regex", help="regular expression; its first group (or match) is the sample id")
    ap.add_argument("--id-case", choices=("keep", "upper", "lower"), default="keep")
    ap.add_argument("--format", choices=("auto", "parquet", "tsv"), default="auto")
    a = ap.parse_args(argv)
    from .export import derive_ids, read_id_map
    folder = os.path.abspath(a.prot)
    t = {}
    for name, cols in TABLES.items():
        if a.no_proteins and name in ("proteins", "abundance"):
            continue
        d = read_input(folder, name)
        if d is None:
            if name in REQUIRED:
                sys.exit("proteomics: %s/proteomics_%s.(parquet|tsv|tsv.gz) not found" % (folder, name))
            continue
        miss = [c for c in cols if c not in d.columns]
        if miss:
            sys.exit("proteomics: proteomics_%s lacks the columns %s" % (name, ", ".join(miss)))
        t[name] = d
        say("read proteomics_%s: %d rows" % (name, len(d)))
    id_map = read_id_map(a.id_map) if a.id_map else None
    names = sorted({str(x) for n in t if "sample_id" in t[n].columns for x in t[n]["sample_id"]})
    ids = dict(zip(names, derive_ids(names, id_map, a.id_regex, a.id_case)))
    for n, d in t.items():
        if "sample_id" in d.columns:
            d["sample_id"] = d["sample_id"].astype(str).map(ids)
        if "fraction" in d.columns:
            d.insert(1, "omic", d.pop("fraction").map(omic))
        if "rank" in d.columns:
            d["rank"] = d["rank"].replace(RANK_ALIAS)
    # shares of the total AUC per sample and omic, when not given
    tot = None
    if "samples" in t and "total_auc" in t["samples"].columns:
        tot = t["samples"].set_index(["sample_id", "omic"])["total_auc"]
    for n in ("taxonomy", "functions"):
        d = t[n]
        if "rel" not in d.columns or d["rel"].isna().all():
            base = tot if tot is not None else d.groupby(["sample_id", "omic"])["auc"].sum()
            d["rel"] = d["auc"] / d.set_index(["sample_id", "omic"]).index.map(base).values
            say("%s: rel computed from the AUC" % n)
    if "samples" not in t:
        s = t["taxonomy"][t["taxonomy"]["rank"] == t["taxonomy"]["rank"].iloc[0]]
        t["samples"] = s.groupby(["sample_id", "omic"], as_index=False).agg(total_auc=("auc", "sum"))
    rel = read_table(a.data, "samples")
    rel_ids = set(rel["sample_id"].astype(str)) if len(rel) else set()
    pids = set(t["taxonomy"]["sample_id"])
    both = pids & rel_ids
    say("%d samples, %d in the release%s" % (len(pids), len(both),
                                             "" if both else " - check --id-case / --id-regex / --id-map"))
    if pids - rel_ids:
        say("not in the release: %s" % ", ".join(sorted(pids - rel_ids)[:20]))
    omics = sorted(set(t["taxonomy"]["omic"]))
    for o in omics:
        say("%s: %d samples" % (o, t["samples"].loc[t["samples"]["omic"] == o, "sample_id"].nunique()))
    for n, d in t.items():
        d = d.copy()
        for c in d.columns:
            if d[c].dtype == object:
                d[c] = d[c].where(d[c].notna(), None).astype("string")
        write_table(d, a.data, "prot_" + n, a.format)
    m = read_manifest(a.data)
    m["proteomics"] = {"source": folder, "samples": len(pids), "samples_in_release": len(both),
                       "fractions": {o: omic_label(o) for o in omics},
                       "ranks": sorted(set(t["taxonomy"]["rank"])), "dbs": sorted(set(t["functions"]["db"])),
                       "proteins": int(len(t["proteins"])) if "proteins" in t else 0,
                       "tables": {"prot_" + k: int(len(v)) for k, v in t.items()}}
    with open(os.path.join(a.data, "manifest.json"), "w") as fh:
        json.dump(m, fh, indent=2, default=str)
    # previous proteomics tables that are not part of this import
    if a.no_proteins:
        for n in ("prot_proteins", "prot_abundance"):
            for ext in (".parquet", ".tsv.gz", ".tsv"):
                p = os.path.join(a.data, n + ext)
                if os.path.exists(p):
                    os.remove(p)
    say("done: prot_* tables and manifest['proteomics'] in %s" % a.data)
    return 0


if __name__ == "__main__":
    sys.exit(main())
