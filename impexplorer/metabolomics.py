"""imp-explorer metabolomics: add metabolomics (or any other per-sample measurement) results to a release.

Input: a folder with three generic tables, each as .parquet, .tsv or .tsv.gz
  metabolomics_assays    assay, label, platform, matrix, acquisition_date, value_type, units, notes, ...
  metabolomics_features  assay, feature_id, name, annotation, class, ... (any further columns are kept)
  metabolomics_values    sample_id, assay, feature_id, value
(metabolomics_samples is optional: sample_id, assay, source_column.) One assay = one platform / panel; features are
the measured compounds; values are concentrations, normalised intensities or any other numbers per sample. The
Expobiome export (metabolomics_export.py) writes this format; other studies can write it from their own result files.

Writes into the release: metab_assays, metab_features, metab_values (sample ids mapped to the release with
--id-map / --id-regex / --id-case), and manifest['metabolomics']. The app shows them on the Metabolomics page and in
Search; 'imp-explorer cohort' analyses them as layers metab.<assay> (same tests, models and PERMANOVA as for the
other layers)."""
import argparse
import json
import os
import sys

import pandas as pd

from .release import read_manifest, read_table, write_table

NAMES = ("assays", "features", "values", "samples")
REQUIRED = {"assays": ("assay",), "features": ("assay", "feature_id"), "values": ("sample_id", "assay", "feature_id",
                                                                                   "value")}


def say(msg):
    sys.stderr.write("[metabolomics] %s\n" % msg)


def read_input(folder, name, prefix="metabolomics_"):
    for ext in (".parquet", ".tsv.gz", ".tsv"):
        p = os.path.join(folder, prefix + name + ext)
        if os.path.exists(p):
            return pd.read_parquet(p) if ext == ".parquet" else pd.read_csv(p, sep="\t", low_memory=False)
    return None


def load(folder):
    t = {n: read_input(folder, n) for n in NAMES}
    for n, cols in REQUIRED.items():
        if t[n] is None:
            sys.exit("metabolomics: %s/metabolomics_%s.(parquet|tsv|tsv.gz) not found" % (folder, n))
        miss = [c for c in cols if c not in t[n].columns]
        if miss:
            sys.exit("metabolomics: metabolomics_%s lacks the columns %s" % (n, ", ".join(miss)))
    for n in ("assays", "features", "values"):
        for c in ("assay", "feature_id", "sample_id"):
            if c in t[n].columns:
                t[n][c] = t[n][c].astype(str)
    t["values"]["value"] = pd.to_numeric(t["values"]["value"], errors="coerce")
    return t


def main(argv=None):
    ap = argparse.ArgumentParser(prog="imp-explorer metabolomics", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True, help="release directory (imp-explorer export)")
    ap.add_argument("--metab", required=True, help="folder with metabolomics_assays / _features / _values")
    ap.add_argument("--assays", default="", help="only these assays (comma-separated)  [all]")
    ap.add_argument("--id-map", help="TSV/CSV: sample id of the metabolomics tables, release sample id")
    ap.add_argument("--id-regex", help="regular expression; its first group (or match) is the sample id")
    ap.add_argument("--id-case", choices=("keep", "upper", "lower"), default="keep")
    ap.add_argument("--format", choices=("auto", "parquet", "tsv"), default="auto")
    a = ap.parse_args(argv)
    from .export import derive_ids, read_id_map
    t = load(os.path.abspath(a.metab))
    want = [x.strip() for x in a.assays.split(",") if x.strip()]
    if want:
        unknown = sorted(set(want) - set(t["assays"]["assay"]))
        if unknown:
            sys.exit("metabolomics: unknown assays %s (have %s)" % (unknown, ", ".join(t["assays"]["assay"])))
        for n in NAMES:
            if t[n] is not None and "assay" in t[n].columns:
                t[n] = t[n][t[n]["assay"].isin(want)]
    id_map = read_id_map(a.id_map) if a.id_map else None
    names = sorted(set(t["values"]["sample_id"]))
    ids = dict(zip(names, derive_ids(names, id_map, a.id_regex, a.id_case)))
    for n in ("values", "samples"):
        if t[n] is not None:
            t[n]["source_id"] = t[n]["sample_id"]
            t[n]["sample_id"] = t[n]["sample_id"].map(ids)
    v = t["values"]
    dup = v.duplicated(["sample_id", "assay", "feature_id"])
    if dup.any():
        sys.exit("metabolomics: %d values twice for the same sample, assay and feature after the id mapping "
                 "(e.g. %s)" % (int(dup.sum()), v[dup].iloc[0][["source_id", "assay", "feature_id"]].tolist()))
    orphan = set(zip(v["assay"], v["feature_id"])) - set(zip(t["features"]["assay"], t["features"]["feature_id"]))
    if orphan:
        say("%d features of metabolomics_values not in metabolomics_features: added with their id as name" % len(orphan))
        t["features"] = pd.concat([t["features"], pd.DataFrame(
            [{"assay": x, "feature_id": f, "name": f} for x, f in sorted(orphan)])], ignore_index=True)
    if "name" not in t["features"].columns:
        t["features"]["name"] = t["features"]["feature_id"]
    t["features"]["name"] = t["features"]["name"].fillna(t["features"]["feature_id"]).astype(str)
    rel = read_table(a.data, "samples")
    rel_ids = set(rel["sample_id"].astype(str)) if len(rel) else set()
    mids = set(v["sample_id"])
    both = mids & rel_ids
    say("%d samples, %d in the release%s" % (len(mids), len(both),
                                             "" if both else " - check --id-case / --id-regex / --id-map"))
    if mids - rel_ids:
        say("not in the release: %s%s" % (", ".join(sorted(mids - rel_ids)[:15]),
                                          " ..." if len(mids - rel_ids) > 15 else ""))
    # per assay: counts for the manifest and the app
    assays = t["assays"].copy()
    cnt = v.dropna(subset=["value"]).groupby("assay").agg(n_samples=("sample_id", "nunique"),
                                                          n_features=("feature_id", "nunique"))
    assays = assays.drop(columns=[c for c in ("n_samples", "n_features") if c in assays.columns]).merge(
        cnt, left_on="assay", right_index=True, how="left")
    assays["n_in_release"] = assays["assay"].map(
        v[v["sample_id"].isin(rel_ids)].groupby("assay")["sample_id"].nunique()).fillna(0).astype(int)
    if "label" not in assays.columns:
        assays["label"] = assays["assay"]
    assays["label"] = assays["label"].fillna(assays["assay"])
    out = {"metab_assays": assays, "metab_features": t["features"],
           "metab_values": v[["sample_id", "assay", "feature_id", "value", "source_id"]]}
    if t["samples"] is not None:
        out["metab_samples"] = t["samples"]
    for name, df in out.items():
        d = df.copy()
        for c in d.columns:
            if d[c].dtype == object:
                d[c] = d[c].where(d[c].notna(), None).astype("string")
        write_table(d, a.data, name, a.format)
    for _, r in assays.iterrows():
        say("%-22s %4s samples (%s in the release) x %5s features  %s" % (
            r["assay"], r.get("n_samples"), r["n_in_release"], r.get("n_features"), r.get("units", "") or ""))
    m = read_manifest(a.data)
    keep = [c for c in ("label", "platform", "matrix", "acquisition_date", "value_type", "units", "notes",
                        "n_samples", "n_features", "n_in_release") if c in assays.columns]
    m["metabolomics"] = {"source": os.path.abspath(a.metab), "samples": len(mids), "samples_in_release": len(both),
                         "assays": {r["assay"]: {k: (None if pd.isna(r[k]) else r[k]) for k in keep}
                                    for _, r in assays.iterrows()},
                         "tables": {k: int(len(df)) for k, df in out.items()}}
    with open(os.path.join(a.data, "manifest.json"), "w") as fh:
        json.dump(m, fh, indent=2, default=str)
    say("done: metab_* tables and manifest['metabolomics'] in %s" % a.data)
    return 0


if __name__ == "__main__":
    sys.exit(main())
