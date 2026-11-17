"""imp-explorer qc: QC status of every sample in every layer of a release (table sample_qc).

One row per sample and layer: status PASS, WARN, FAIL or missing, the reasons and where the flag comes from. Layers are
the omics of the release: mg, mt (IMP3), fractions added with 'imp-explorer fraction' (mg_ext, mt_ext_remap, ...),
the sRNA-IMP fractions (e.g. srna_int, srna_ext), the protein fractions
(mp_<fraction>) and the metabolomics assays (metab.<assay>). Sources, in this order:

  auto     rules on the release itself: no data in a layer -> missing; mg / mt read pairs after trimming below
           --min-mg-pairs / --min-mt-pairs -> FAIL; outliers flagged by 'imp-explorer cohort' (cohort_qc) -> WARN
  srna     the QC verdict of sRNA-IMP per library (PASS / WARN / FAIL, with its rules)
  manual   --flags FILE: sample_id, layer, status, reason (TSV/CSV). layer '*' = every layer of that sample; a manual
           row replaces the automatic status of that sample and layer (e.g. the exclusions of a paper, or PASS to
           clear an automatic flag)

All samples stay in the release. The app marks the flagged samples and has a switch to leave them out (all samples /
without FAIL / PASS only), per layer.

  imp-explorer qc --data release_v0.4 [--flags paper_exclusions.tsv] [--list-layers]
"""
import argparse
import json
import os
import sys

import pandas as pd

from .release import read_manifest, read_table, table_path, write_table

LEVELS = ("PASS", "WARN", "FAIL", "missing")
RANK = {"PASS": 0, "WARN": 1, "FAIL": 2, "missing": 3}


def say(msg):
    sys.stderr.write("[qc] %s\n" % msg)


def _have(release, name):
    return table_path(release, name) is not None


def layer_samples(release):
    """{layer: set(sample_id with data)} from the release tables"""
    out = {}
    smp = read_table(release, "samples")
    smp["sample_id"] = smp["sample_id"].astype(str)
    qc = read_table(release, "qc")
    for o in ("mg", "mt"):
        col = "has_" + o
        have = set(smp.loc[smp[col].astype(str).str.lower().isin(("true", "1", "yes")), "sample_id"]) \
            if col in smp.columns else set()
        if not have and not qc.empty:
            have = set(qc.loc[qc["omic"] == o, "sample_id"].astype(str))
        if have or col in smp.columns:
            out[o] = have
    for col in [c for c in smp.columns if str(c).startswith("has_") and c not in ("has_mg", "has_mt")]:
        o = str(col)[4:]                          # fractions added with 'imp-explorer fraction' (mg_ext, ...)
        out[o] = set(smp.loc[smp[col].astype(str).str.lower().isin(("true", "1", "yes")), "sample_id"])
    if _have(release, "srna_samples"):
        s = read_table(release, "srna_samples")
        for o, x in s.groupby("omic"):
            out[str(o)] = set(x["sample_id"].astype(str))
    if _have(release, "prot_samples"):
        p = read_table(release, "prot_samples")
        for o, x in p.groupby("omic"):
            out[str(o)] = set(x["sample_id"].astype(str))
    if _have(release, "metab_values"):
        v = read_table(release, "metab_values")
        v = v.dropna(subset=["value"]) if "value" in v.columns else v
        for a, x in v.groupby("assay"):
            out["metab.%s" % a] = set(x["sample_id"].astype(str))
    return out


def read_flags(path):
    sep = "," if path.endswith(".csv") else "\t"
    f = pd.read_csv(path, sep=sep, dtype=str, comment="#").fillna("")
    f.columns = [c.strip().lower() for c in f.columns]
    miss = [c for c in ("sample_id", "layer", "status") if c not in f.columns]
    if miss:
        sys.exit("qc: %s lacks the columns %s (need sample_id, layer, status[, reason])" % (path, ", ".join(miss)))
    if "reason" not in f.columns:
        f["reason"] = ""
    f["status"] = [x.upper() if x.upper() in ("PASS", "WARN", "FAIL") else x.lower() for x in f["status"].str.strip()]
    for c in ("sample_id", "layer", "reason"):
        f[c] = f[c].str.strip()
    bad = sorted(set(f["status"]) - set(LEVELS))
    if bad:
        sys.exit("qc: unknown status %s in %s (use %s)" % (bad, path, " / ".join(LEVELS)))
    return f


def build(release, flags=None, min_mg_pairs=1e6, min_mt_pairs=1e6):
    """DataFrame sample_qc (sample_id, layer, status, reasons, source) and a list of notes"""
    notes = []
    smp = read_table(release, "samples")
    ids = list(smp["sample_id"].astype(str))
    have = layer_samples(release)
    cells = {}                                  # (sid, layer) -> [status, [reasons], [sources]]

    def put(sid, layer, status, reason, source):
        c = cells.setdefault((sid, layer), ["PASS", [], []])
        if RANK[status] > RANK[c[0]]:
            c[0] = status
        if reason:
            c[1].append(reason)
        if source not in c[2] and (reason or status != "PASS"):
            c[2].append(source)

    extra = set()
    for layer, got in have.items():
        extra |= got - set(ids)
        for sid in ids:
            put(sid, layer, "PASS" if sid in got else "missing", "" if sid in got else "no data", "auto")
    if extra:
        notes.append("%d samples have layer data but are not in the release samples (left out): %s"
                     % (len(extra), ", ".join(sorted(extra)[:8])))
    # mg / mt read pairs after trimming
    qc = read_table(release, "qc")
    if not qc.empty:
        for o, lim in (("mg", min_mg_pairs), ("mt", min_mt_pairs)):
            v = qc[(qc["omic"] == o) & (qc["metric"] == "trimmed:pairs")].groupby("sample_id")["value"].first()
            for sid, x in v.items():
                sid = str(sid)
                if o in have and sid in have[o] and pd.notna(x) and float(x) < lim:
                    put(sid, o, "FAIL", "%s read pairs %.2g < %.2g" % (o, float(x), lim), "auto")
    # outliers from the cohort analysis (cohort_qc: one row per sample, reasons as text)
    if _have(release, "cohort_qc"):
        c = read_table(release, "cohort_qc")
        for _, r in c.iterrows():
            for part in str(r.get("reasons", "") or "").split("; "):
                for o in ("mg", "mt"):
                    if part.startswith("outlier %s_" % o) and o in have and str(r["sample_id"]) in have[o]:
                        put(str(r["sample_id"]), o, "WARN", part, "auto")
    # sRNA-IMP verdicts
    if _have(release, "srna_samples"):
        s = read_table(release, "srna_samples")
        for _, r in s.iterrows():
            lv = str(r.get("qc", "") or "").strip().upper()
            if lv in ("WARN", "FAIL"):
                why = str(r.get("rules", "") or r.get("details", "") or "sRNA-IMP QC %s" % lv)
                put(str(r["sample_id"]), str(r["omic"]), lv, why if why != "nan" else "sRNA-IMP QC %s" % lv, "srna")
    # manual flags replace the automatic status
    if flags:
        f = read_flags(flags)
        unknown_s = sorted(set(f["sample_id"]) - set(ids))
        unknown_l = sorted(set(f["layer"]) - set(have) - {"*"})
        if unknown_s:
            notes.append("%d manual flags for samples not in the release (ignored): %s"
                         % (len(unknown_s), ", ".join(unknown_s[:8])))
        if unknown_l:
            notes.append("manual flags for unknown layers (ignored): %s; layers: %s"
                         % (", ".join(unknown_l), ", ".join(sorted(have))))
        for _, r in f.iterrows():
            if r["sample_id"] not in ids:
                continue
            for layer in (sorted(have) if r["layer"] == "*" else [r["layer"]]):
                if layer not in have:
                    continue
                c = cells.setdefault((r["sample_id"], layer), ["PASS", [], []])
                if c[0] == "missing" and r["status"] != "missing":
                    continue                                   # no data stays missing
                c[0] = r["status"]
                c[1] = [x for x in [r["reason"] or "manual %s" % r["status"]] + c[1] if x]
                c[2] = ["manual"] + [x for x in c[2] if x != "manual"]
    rows = [{"sample_id": sid, "layer": layer, "status": c[0], "reasons": "; ".join(dict.fromkeys(c[1])),
             "source": ", ".join(c[2])} for (sid, layer), c in sorted(cells.items(), key=lambda kv: (kv[0][1], kv[0][0]))]
    return pd.DataFrame(rows, columns=["sample_id", "layer", "status", "reasons", "source"]), notes


def summary(t):
    if t.empty:
        return pd.DataFrame(columns=["layer"] + list(LEVELS))
    s = t.pivot_table(index="layer", columns="status", values="sample_id", aggfunc="count", fill_value=0)
    return s.reindex(columns=list(LEVELS), fill_value=0).reset_index()


def main(argv=None):
    ap = argparse.ArgumentParser(prog="imp-explorer qc", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True, help="release directory")
    ap.add_argument("--flags", help="manual flags: sample_id, layer, status (PASS/WARN/FAIL/missing), reason")
    ap.add_argument("--min-mg-pairs", type=float, default=1e6, help="mg read pairs after trimming  [1e6]")
    ap.add_argument("--min-mt-pairs", type=float, default=1e6, help="mt read pairs after trimming  [1e6]")
    ap.add_argument("--list-layers", action="store_true", help="only print the layers of the release and stop")
    ap.add_argument("--format", default="auto", choices=("auto", "parquet", "tsv"))
    a = ap.parse_args(argv)
    if not os.path.isdir(a.data):
        sys.exit("qc: release directory %s not found" % a.data)
    if a.list_layers:
        for layer, got in sorted(layer_samples(a.data).items()):
            print("%s\t%d samples" % (layer, len(got)))
        return 0
    t, notes = build(a.data, a.flags, a.min_mg_pairs, a.min_mt_pairs)
    for n in notes:
        say("warning: " + n)
    write_table(t, a.data, "sample_qc", a.format)
    s = summary(t)
    say("sample_qc: %d samples x %d layers" % (t["sample_id"].nunique(), t["layer"].nunique()))
    for _, r in s.iterrows():
        say("  %-24s %s" % (r["layer"], "  ".join("%s %d" % (lv, r[lv]) for lv in LEVELS)))
    m = read_manifest(a.data)
    m["sample_qc"] = {"flags": os.path.abspath(a.flags) if a.flags else "", "min_mg_pairs": a.min_mg_pairs,
                      "min_mt_pairs": a.min_mt_pairs,
                      "counts": {r["layer"]: {lv: int(r[lv]) for lv in LEVELS} for _, r in s.iterrows()}}
    with open(os.path.join(a.data, "manifest.json"), "w") as fh:
        json.dump(m, fh, indent=2, default=str)
    return 0
