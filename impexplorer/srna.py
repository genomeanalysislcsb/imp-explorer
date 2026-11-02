"""imp-explorer srna: add the results of an sRNA-IMP cohort analysis to a release (optional layer).

Reads the Cohort folder written by sRNA-IMP (srna_cohort.py: tables/ and stats/) and writes srna_* tables into the
release. Nothing is recomputed: the numbers are those of the sRNA-IMP cohort report.

  per sample   srna_samples (fraction, group, block, QC), srna_qc (library metrics, long),
               srna_features (count matrices: human miRNAs, taxa from rRNA reads, RNA classes, Rfam, tRNA, functions)
  cohort       srna_tests (group tests per feature and metric), srna_models (+ _overview), srna_permanova,
               srna_dispersion, srna_ordination (NMDS), srna_covariates

Layers are '<fraction>.<set>', e.g. srna_ext.genus, srna_int.mirna, srna_ext.metrics (fraction = the sRNA-IMP type,
lower case, '-' -> '_').
"""
import argparse
import glob
import json
import os
import re
import sys

import numpy as np
import pandas as pd

from .release import read_manifest, write_table

SET_LABEL = {
    "mirna": "human miRNAs", "genus": "genera (rRNA reads)", "family": "families (rRNA reads)",
    "species": "species (rRNA reads)", "rfam_type": "Rfam types", "rfam_family": "Rfam families",
    "rna_class": "RNA classes", "trna_aa": "tRNA amino acids", "fn_kegg_ko": "KEGG KO (reads on CDS)",
    "fn_kegg_module": "KEGG modules (reads on CDS)", "fn_kegg_pathway": "KEGG pathways (reads on CDS)",
    "fn_pfam": "Pfam (reads on CDS)", "fn_cog": "COG (reads on CDS)", "metrics": "library metrics",
    "asidx_kegg_ko": "antisense index KEGG KO", "asidx_kegg_pathway": "antisense index KEGG pathways",
}
SKIP_SETS = re.compile(r"^(fnas_|antisense_index)")


def say(msg):
    sys.stderr.write("[srna] %s\n" % msg)


def fraction_key(t):
    return re.sub(r"[^a-z0-9]+", "_", str(t).lower()).strip("_")


def label(layer, types):
    frac, _, s = layer.partition(".")
    return "%s %s" % (types.get(frac, frac), SET_LABEL.get(s, s.replace("_", " ")))


def read_tsv(path):
    return pd.read_csv(path, sep="\t", dtype={"sample": str, "feature": str}, keep_default_na=True,
                       na_values=["NA", "nan", ""], low_memory=False)


def tidy_comparison(lab):
    """sRNA-IMP labels -> labels as in imp-explorer cohort: 'PD vs HC (within site: Kassel)' -> 'PD vs HC (Kassel)',
    'HC Berlin vs HC Kassel (site effect within HC)' -> 'HC: Berlin vs Kassel'"""
    lab = str(lab)
    m = re.match(r"^(\S+) (.+?) vs \1 (.+?) \(.*effect within \1\)$", lab)
    if m:
        return "%s: %s vs %s" % (m.group(1), m.group(2), m.group(3)), "block"
    lab2 = re.sub(r"\(within [^:]+: ([^)]+)\)", r"(\1)", lab)
    return lab2, ("pooled" if "pooled" in lab else "within")


class SrnaCohort:
    def __init__(self, cohort_dir, ids, block, sig_q):
        self.dir = cohort_dir
        self.tdir, self.sdir = os.path.join(cohort_dir, "tables"), os.path.join(cohort_dir, "stats")
        if not (os.path.isdir(self.tdir) and os.path.isdir(self.sdir)):
            sys.exit("[srna] %s has no tables/ and stats/ (the Cohort folder of sRNA-IMP)" % cohort_dir)
        self.ids, self.block, self.sig_q = ids, block, sig_q
        self.types = {}                                    # fraction key -> sRNA-IMP type
        self.out = {}

    def sid(self, s):
        return self.ids(str(s))

    # ------------------------------------------------------------------ per sample
    def samples(self):
        st = read_tsv(os.path.join(self.tdir, "cohort.sample_stats.tsv"))
        meta_p = os.path.join(self.tdir, "cohort.metadata.tsv")
        meta_cols = list(read_tsv(meta_p).columns) if os.path.exists(meta_p) else []
        self.types = {fraction_key(t): t for t in st["type"].dropna().unique()}
        st["sample_id"] = st["sample"].map(self.sid)
        st["omic"] = st["type"].map(fraction_key)
        block = self.block if self.block in st.columns else ("site" if "site" in st.columns else "")
        self.block = block
        smp = pd.DataFrame({"sample_id": st["sample_id"], "omic": st["omic"], "type": st["type"],
                            "group": st.get("group", pd.Series("all", index=st.index)).astype(str),
                            "block": st[block].astype(str) if block else "all",
                            "qc": st.get("qc", ""), "in_statistics": st.get("in_statistics", "")})
        flags = os.path.join(self.tdir, "cohort.qc_flags.tsv")
        if os.path.exists(flags):
            f = read_tsv(flags)
            f["sample_id"] = f["sample"].map(self.sid)
            f["omic"] = f["type"].map(fraction_key)
            smp = smp.merge(f[["sample_id", "omic", "rules", "details"]], on=["sample_id", "omic"], how="left")
        self.out["srna_samples"] = smp
        # library metrics: numeric columns that are not metadata
        skip = set(meta_cols) | {"sample", "type", "group", "qc", "in_statistics", "sample_id", "omic"}
        num = [c for c in st.columns if c not in skip and pd.api.types.is_numeric_dtype(st[c])]
        q = st[["sample_id", "omic"] + num].melt(id_vars=["sample_id", "omic"], var_name="metric", value_name="value")
        self.out["srna_qc"] = q.dropna(subset=["value"])
        self.metrics = num
        say("%d libraries (%s), %d metrics" % (len(smp), ", ".join(self.types.values()), len(num)))

    def features(self):
        parts = []
        for p in sorted(glob.glob(os.path.join(self.tdir, "cohort.*.counts.tsv"))):
            m = re.match(r"cohort\.([^.]+)\.([A-Za-z0-9_]+)\.counts\.tsv$", os.path.basename(p))
            if not m or SKIP_SETS.match(m.group(2)):
                continue
            frac, s = fraction_key(m.group(1)), m.group(2)
            d = read_tsv(p).set_index("feature").copy()
            d.columns = [self.sid(c) for c in d.columns]
            long = d.reset_index().melt(id_vars="feature", var_name="sample_id", value_name="reads")
            long = long[long["reads"] > 0]
            tot = long.groupby("sample_id")["reads"].transform("sum")
            long["rel"] = long["reads"] / tot
            long.insert(0, "set", s)
            long.insert(0, "omic", frac)
            parts.append(long[["sample_id", "omic", "set", "feature", "reads", "rel"]])
        self.out["srna_features"] = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(
            columns=["sample_id", "omic", "set", "feature", "reads", "rel"])
        say("features: %d rows in %d sets" % (len(self.out["srna_features"]), len(parts)))

    # ------------------------------------------------------------------ cohort statistics
    def stat_files(self, pattern):
        rx = re.compile(pattern)
        for p in sorted(glob.glob(os.path.join(self.sdir, "cohort.*.tsv"))):
            m = rx.match(os.path.basename(p))
            if m:
                yield p, m

    def comparisons(self):
        """file stem of matched tests -> (label, kind, n_a, n_b), from the matched PERMANOVA files"""
        out = {}
        for p, m in self.stat_files(r"cohort\.([^.]+)\.matched\.([^.]+)\.permanova\.tsv$"):
            d = read_tsv(p)
            if len(d):
                lab, kind = tidy_comparison(d["comparison"].iloc[0])
                out[(fraction_key(m.group(1)), m.group(2))] = (lab, kind, d["n_a"].iloc[0], d["n_b"].iloc[0])
        return out

    def tests(self):
        rows = []
        comps = self.comparisons()
        groups = sorted(self.out["srna_samples"]["group"].unique())
        # Kruskal-Wallis across groups per feature (by_group) and per metric (metrics_by_group)
        for p, m in self.stat_files(r"cohort\.([^.]+)\.(?:([A-Za-z0-9_]+)\.by_group|(metrics)_by_group)\.tsv$"):
            frac, s = fraction_key(m.group(1)), (m.group(2) or "metrics")
            if SKIP_SETS.match(s):
                continue
            d = read_tsv(p)
            key = "feature" if "feature" in d.columns else d.columns[0]
            mean = [c for c in d.columns if c.startswith(("mean_rel_", "median_"))]
            summ = d[mean].apply(lambda r: "; ".join("%s %.4g" % (c.split("_")[-1], v) for c, v in r.items()
                                                       if pd.notna(v)), axis=1) if mean else ""
            rows.append(pd.DataFrame({"layer": "%s.%s" % (frac, s), "feature": d[key].astype(str),
                                      "comparison": "all groups", "kind": "kruskal",
                                      "n_a": d["n"] if "n" in d.columns else np.nan, "effect": d.get("H"),
                                      "p": d["p"], "q": d["q"], "prevalence": d.get("prevalence"),
                                      "means" if s != "metrics" else "medians": summ}))
        # matched comparisons (within blocks, block contrast in the reference)
        for p, m in self.stat_files(r"cohort\.([^.]+)\.matched\.([^.]+)\.([A-Za-z0-9_]+)\.tsv$"):
            frac, stem, s = fraction_key(m.group(1)), m.group(2), m.group(3)
            if s in ("permanova",) or SKIP_SETS.match(s) or stem.endswith(".antisense_index"):
                continue
            lab, kind, na, nb = comps.get((frac, stem), (stem.replace("_", " "), "within", np.nan, np.nan))
            d = read_tsv(p)
            key = "feature" if "feature" in d.columns else d.columns[0]
            eff = "P(a>b)" if "P(a>b)" in d.columns else None
            r = pd.DataFrame({"layer": "%s.%s" % (frac, s), "feature": d[key].astype(str), "comparison": lab,
                              "kind": kind, "effect": d[eff] if eff else np.nan, "p": d["p"], "q": d["q"]})
            r["n_a"] = d["n_a"] if "n_a" in d.columns else na
            r["n_b"] = d["n_b"] if "n_b" in d.columns else nb
            for a, b in (("median_a", "median_a"), ("median_b", "median_b"), ("mean_rel_a", "mean_a"),
                         ("mean_rel_b", "mean_b"), ("prevalence", "prevalence")):
                if a in d.columns:
                    r[b] = d[a]
            rows.append(r)
        t = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
        if len(t):
            t["name"] = ""
        self.out["srna_tests"] = t
        say("tests: %d rows, comparisons: %s" % (len(t), ", ".join(sorted(t["comparison"].unique())) if len(t) else "-"))
        self.groups = groups

    def models(self):
        mt = read_tsv(os.path.join(self.sdir, "cohort.models.tsv")) if os.path.exists(
            os.path.join(self.sdir, "cohort.models.tsv")) else pd.DataFrame()
        forms, notes = {}, {}
        for r in mt.itertuples():
            forms[(r.model, fraction_key(r.type))] = r.formula
            notes[r.model] = "" if str(getattr(r, "notes", "-")) in ("-", "nan") else str(r.notes)
        rows, over = [], []
        for p, m in self.stat_files(r"cohort\.([^.]+)\.model_(M\d+)\.([A-Za-z0-9_]+)\.tsv$"):
            frac, model, s = fraction_key(m.group(1)), m.group(2), m.group(3)
            if s == "design":
                continue
            d = read_tsv(p)
            if not len(d):
                continue
            form = forms.get((model, frac), "")
            first = re.split(r"[+|]", form)[0].strip() if form else ""
            focus = d["term"].astype(str).map(lambda x: bool(first) and (x == first or x.startswith(first + ":") or
                                                                         x.startswith(first + " (per SD")))
            layer = "%s.%s" % (frac, s)
            r = pd.DataFrame({"model": model, "layer": layer, "feature": d["feature"].astype(str),
                              "n": d["n"], "term": d["term"].astype(str), "focus": focus, "estimate": d["estimate"],
                              "se": d["se"], "t": d.get("t"), "p": d["p"], "q": d["q"],
                              "name": d.get("name", pd.Series("", index=d.index)).fillna("")})
            rows.append(r)
            over.append({"model": model, "formula": form, "note": notes.get(model, ""), "layer": layer,
                         "n": int(d["n"].max()), "features": int(d["feature"].nunique()),
                         "features_q_focus": int(r.loc[r["focus"] & (r["q"] < self.sig_q), "feature"].nunique()),
                         "info": ""})
        self.out["srna_models"] = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
        self.out["srna_models_overview"] = pd.DataFrame(over)
        say("models: %d rows, %d model x layer fits" % (len(self.out["srna_models"]), len(over)))
        self.model_specs = [{"model": k, "formula": forms.get((k, f), ""), "note": notes.get(k, "")}
                            for k, f in sorted({(k, f) for k, f in forms}, key=lambda x: (int(x[0][1:]), x[1]))]

    def beta(self):
        pm, od = [], []
        for p, m in self.stat_files(r"cohort\.([^.]+)\.([A-Za-z0-9_]+)\.permanova\.tsv$"):
            frac, s = fraction_key(m.group(1)), m.group(2)
            d = read_tsv(p)
            for r in d.itertuples():
                pm.append({"layer": "%s.%s" % (frac, s), "comparison": "%s (all samples)" % r.variable,
                           "kind": "groups", "n": r.n, "F": r.pseudo_F, "R2": r.R2, "p": r.p})
        sets = {}
        for p, m in self.stat_files(r"cohort\.([^.]+)\.matched\.([^.]+)\.permanova\.tsv$"):
            frac = fraction_key(m.group(1))
            d = read_tsv(p)
            for r in d.itertuples():
                lab, kind = tidy_comparison(r.comparison)
                # the matched PERMANOVA of sRNA-IMP is on the genus profile
                pm.append({"layer": "%s.genus" % frac, "comparison": lab, "kind": kind, "n": r.n_a + r.n_b,
                           "F": r.pseudo_F, "R2": r.R2, "p": r.p})
        smp = self.out["srna_samples"]
        blk = {(a, b): c for a, b, c in zip(smp["sample_id"], smp["omic"], smp["block"])}
        for p, m in self.stat_files(r"cohort\.([^.]+)\.([A-Za-z0-9_]+)\.nmds\.tsv$"):
            frac, s = fraction_key(m.group(1)), m.group(2)
            d = read_tsv(p)
            sids = d["sample"].map(self.sid)
            od.append(pd.DataFrame({"layer": "%s.%s" % (frac, s), "sample_id": sids, "axis1": d["NMDS1"],
                                    "axis2": d["NMDS2"], "method": "NMDS", "group": d["group"].astype(str),
                                    "block": [blk.get((x, frac), "all") for x in sids],
                                    "explained1": d["stress"], "explained2": np.nan}))
        self.out["srna_permanova"] = pd.DataFrame(pm)
        self.out["srna_ordination"] = pd.concat(od, ignore_index=True) if od else pd.DataFrame()
        dp = os.path.join(self.sdir, "cohort.heterogeneity.dispersion.tsv")
        if os.path.exists(dp):
            d = read_tsv(dp)
            self.out["srna_dispersion"] = pd.DataFrame({
                "layer": d["type"].map(fraction_key) + "." + d["set"].astype(str),
                "comparison": d["comparison"], "kind": "within", "n_a": d["n_cases"], "n_b": d["n_controls"],
                "median_dist_a": d["median_dist_cases"], "median_dist_b": d["median_dist_controls"],
                "ratio": d["ratio"], "p": d["p"], "q": d["q"]})
        cv = os.path.join(self.sdir, "cohort.covariates_by_group.tsv")
        if os.path.exists(cv):
            self.out["srna_covariates"] = read_tsv(cv)

    def run(self):
        self.samples()
        self.features()
        self.tests()
        self.models()
        self.beta()
        return self.out


def main(argv=None):
    ap = argparse.ArgumentParser(prog="imp-explorer srna", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True, help="release directory (imp-explorer export)")
    ap.add_argument("--srna-cohort", required=True, help="Cohort folder of sRNA-IMP (with tables/ and stats/)")
    ap.add_argument("--block", default="site", help="block column of the sRNA metadata  [site]")
    ap.add_argument("--sig-q", type=float, default=0.1, help="FDR threshold for the hit counts  [0.1]")
    ap.add_argument("--id-map", help="TSV/CSV: sRNA sample name, release sample id")
    ap.add_argument("--id-regex", help="regular expression; its first group (or match) is the sample id")
    ap.add_argument("--id-case", choices=("keep", "upper", "lower"), default="keep")
    ap.add_argument("--format", choices=("auto", "parquet", "tsv"), default="auto")
    a = ap.parse_args(argv)
    from .export import derive_ids, read_id_map
    id_map = read_id_map(a.id_map) if a.id_map else None
    cache = {}

    def ids(name):
        if name not in cache:
            cache[name] = derive_ids([name], id_map, a.id_regex, a.id_case)[0]
        return cache[name]
    C = SrnaCohort(os.path.abspath(a.srna_cohort), ids, a.block, a.sig_q)
    out = C.run()
    rel = read_table_ids(a.data)
    srna_ids = set(out["srna_samples"]["sample_id"])
    both = srna_ids & rel
    say("%d sRNA samples, %d in the release%s" % (len(srna_ids), len(both),
                                                   "" if both else " - check --id-case / --id-regex / --id-map"))
    for name, df in out.items():
        if isinstance(df, pd.DataFrame) and len(df):
            d = df.copy()
            for c in d.columns:
                if d[c].dtype == object:
                    d[c] = d[c].astype(str).replace({"nan": ""})
            write_table(d, a.data, name, a.format)
    layers = sorted({*out["srna_tests"].get("layer", pd.Series(dtype=str)),
                     *out["srna_models_overview"].get("layer", pd.Series(dtype=str)),
                     *out["srna_ordination"].get("layer", pd.Series(dtype=str))})
    m = read_manifest(a.data)
    m["srna"] = {"source": os.path.abspath(a.srna_cohort), "fractions": C.types, "block": C.block,
                 "groups": C.groups, "sig_q": a.sig_q, "metrics": C.metrics, "models": C.model_specs,
                 "layers": {k: {"label": label(k, C.types)} for k in layers},
                 "samples": len(srna_ids), "samples_in_release": len(both),
                 "tables": {k: int(len(v)) for k, v in out.items() if isinstance(v, pd.DataFrame)}}
    with open(os.path.join(a.data, "manifest.json"), "w") as fh:
        json.dump(m, fh, indent=2, default=str)
    say("done: srna_* tables and manifest['srna'] in %s" % a.data)
    return 0


def read_table_ids(release):
    from .release import read_table
    s = read_table(release, "samples")
    return set(s["sample_id"].astype(str)) if len(s) else set()


if __name__ == "__main__":
    sys.exit(main())
