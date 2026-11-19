"""imp-explorer cohort: cohort analysis of an IMP3 release (QC and covariates, taxa and diversity, group comparisons
within blocks, covariate-adjusted models, MAGs, functions, metatranscriptome / metagenome ratios).

Writes <out>/tables/*.tsv, <out>/stats/*.tsv, <out>/cohort.report.pdf + .html and cohort_* tables into the release
for the app's Cohort page. Same statistics and settings as the sRNA-IMP cohort analysis; a study settings file with --config (also the
COHORT_* block of an sRNA-IMP config).
"""
import json
import os
import sys
import time

import pandas as pd

from .. import __version__
from ..release import read_manifest, table_path, write_table
from . import settings
from .data import Cohort, METRIC_LABEL
from .engine import Analysis
from .stats import load_term_names

RELEASE_TABLES = {   # cohort result -> release table read by the app
    "models": "cohort_models", "models_overview": "cohort_models_overview", "feature_tests": "cohort_tests",
    "metrics_tests": "cohort_metric_tests", "permanova": "cohort_permanova", "dispersion": "cohort_dispersion",
    "ordination": "cohort_ordination", "sample_metrics": "cohort_samples", "covariates": "cohort_covariates",
    "qc_flags": "cohort_qc", "mag_prevalence": "cohort_mag_prevalence",
    "paper_anova": "cohort_paper_anova", "paper_permanova": "cohort_paper_permanova",
    "paper_variance": "cohort_paper_variance", "concord_presence": "cohort_concord_presence",
    "concord_alpha": "cohort_concord_alpha", "concord_mantel": "cohort_concord_mantel",
    "concord_features": "cohort_concord_features", "concord_within": "cohort_concord_within",
    "fraction_diff": "cohort_fraction_diff", "fraction_summary": "cohort_fraction_summary"}


def write_outputs(C, A, res, a):
    tdir, sdir = os.path.join(a.out, "tables"), os.path.join(a.out, "stats")
    os.makedirs(tdir, exist_ok=True)
    os.makedirs(sdir, exist_ok=True)
    meta = C.samples.copy()
    meta.drop(columns=[c for c in ("imp_dir",) if c in meta.columns]).to_csv(
        os.path.join(tdir, "cohort.metadata.tsv"), sep="\t")
    for key, L in C.layers.items():
        L.M.T.to_csv(os.path.join(tdir, "cohort.%s.%s.tsv" % (key, {"counts": "counts", "metab": "values", "auc": "auc"}.get(L.kind, "log2ratio"))),
                     sep="\t", float_format="%.6g")
    for name, df in res.items():
        if isinstance(df, pd.DataFrame):
            df.to_csv(os.path.join(sdir if name not in ("sample_metrics", "qc_flags") else tdir,
                                   "cohort.%s.tsv" % name), sep="\t", index=False, float_format="%.6g")
    with open(os.path.join(sdir, "cohort.findings.tsv"), "w") as fh:
        fh.write("area\tfinding\n")
        for area, txt in res["findings"]:
            fh.write("%s\t%s\n" % (area, txt.replace("\t", " ")))


def settings_dict(C, A, a):
    return {"name": a.name or C.study, "created": time.strftime("%Y-%m-%dT%H:%M:%S"), "impexplorer": __version__,
            "group_col": C.group_col, "groups": C.groups, "reference": C.ref, "block": C.block,
            "covariates": C.covariates, "models": [{k: (sorted(v[1]) if k == "subset" and v else v)
                                                    for k, v in s.items() if k != "terms"}
                                                   for s in A.res["model_specs"]],
            "layers": {k: {"label": L.label, "omic": L.omic, "kind": L.kind, "db": L.db, "samples": int(L.M.shape[0]),
                           "features": int(L.M.shape[1])} for k, L in C.layers.items()},
            "metrics": {m: METRIC_LABEL.get(m, m) for m in C.metric_names},
            "permutations": a.permutations, "sig_q": a.sig_q, "prevalence": a.prevalence, "qc_drop": a.qc_drop,
            "nmds": a.nmds, "nmds_starts": getattr(a, "nmds_starts", None),
            "paper": {"model": a.paper_model, "transform": a.paper_transform, "qc": a.paper_qc}
            if a.paper_model else {},
            "concordance": a.concordance, "fraction_pairs": a.fraction_pairs, "fraction_min_r": a.fraction_min_r,
            "contig_min_share": a.contig_min_share,
            "samples": len(C.samples), "kept": len(C.kept), "excluded": C.excluded,
            "config": os.path.abspath(a.config) if a.config else "", "from_config": a.from_config, "notes": C.log}


def main(argv=None):
    a = settings.parse(argv)
    if a.write_config:
        sys.stderr.write("[cohort] settings template written to %s\n" % a.write_config)
        return 0
    t0 = time.time()
    if a.term_names:
        load_term_names(a.term_names)
    C = Cohort(a)
    C.note("%d samples, %d kept; groups %s (reference %s); block %s" % (
        len(C.samples), len(C.kept), ",".join(C.groups) or "–", C.ref or "–", C.block or "–"))
    default = a.layers == settings.DEFAULT_LAYERS
    if default and table_path(a.data, "metab_values"):
        a.layers += ",metab"                  # metabolomics of the release (imp-explorer metabolomics): every assay
    fdbs = read_manifest(a.data).get("sources", {}).get("function_dbs") or []
    if default and read_manifest(a.data).get("sources", {}).get("contig_taxonomy"):
        a.layers += ",contig.genus.mg,contig.genus.mt,contig.species.mg"     # contig taxonomy (export --contig-taxonomy)
    if default:                               # Pfam / COG of mg and mt, when exported (--function-dbs)
        a.layers += "".join(",%s.mg,%s.mt" % (d, d) for d in ("pfam", "cog") if d in fdbs)
    fr = read_manifest(a.data).get("fractions", {})
    if default and fr:                        # fractions added with 'imp-explorer fraction'
        for o in fr:
            if o.endswith("_remap"):          # reads on the main assembly: functions only
                a.layers += "".join(",%s.%s" % (d, o) for d in ("kegg_module", "kegg_pathway", "cazy", "kegg_ko") +
                                    tuple(x for x in ("pfam", "cog") if x in fdbs))
            else:
                tb = fr[o].get("tables") or ("taxonomy", "bin_abundance", "functions")
                a.layers += "".join(",%s.%s" % (k, o) for k in
                                    (("kraken.genus", "kraken.species") if "taxonomy" in tb else ()) +
                                    (("contig.genus", "contig.species") if fr[o].get("contig_taxonomy") else ()) +
                                    (("mag.genus",) if "bin_abundance" in tb else ()) +
                                    (("kegg_module", "kegg_pathway", "cazy", "kegg_ko") +
                                     tuple(x for x in ("pfam", "cog") if x in fdbs) if "functions" in tb else ()))
    if default and table_path(a.data, "prot_samples"):
        a.layers += ",prot"                   # metaproteomics: genus and KEGG layers of every fraction
    if a.concordance:                         # the layers the concordance sets need
        from .concord import needed_layers, parse
        have = [k.strip() for k in C.expand_metab(a.layers).split(",")]
        extra = [k for k in needed_layers(parse(a.concordance)) if k not in have]
        if extra:
            a.layers += "," + ",".join(extra)
    if a.fraction_pairs:                      # the layers the intra / extra pairs need
        from .fracdiff import needed_layers as fr_needed, parse as fr_parse
        have = [k.strip() for k in C.expand_metab(a.layers).split(",")]
        extra = [k for k in fr_needed(fr_parse(a.fraction_pairs)) if k not in have]
        if extra:
            a.layers += "," + ",".join(dict.fromkeys(extra))
    C.load_layers(a.layers)
    C.note("layers: " + ", ".join("%s (%d x %d)" % (k, L.M.shape[0], L.M.shape[1]) for k, L in C.layers.items()))
    if not C.layers:
        why = C.qc_flags.loc[C.excluded, "reasons"].value_counts().head(3) if C.excluded else None
        C.note("no layer to analyse: %d of %d samples kept by the sample QC (--qc-drop %s)%s; check the release "
               "(imp-explorer export) and the metadata" % (
                   len(C.kept), len(C.samples), a.qc_drop,
                   "; most frequent reasons: " + "; ".join("%s (%d)" % (r, n) for r, n in why.items())
                   if why is not None and len(why) else ""))
    A = Analysis(C)
    res = A.run()
    if a.paper_model:
        from .paper import run as run_paper
        res.update(run_paper(C, a))
    if a.concordance:
        from .concord import run as run_concord
        res.update(run_concord(C, a))
    if a.fraction_pairs:
        from .fracdiff import run as run_fracdiff
        res.update(run_fracdiff(C, a))
    write_outputs(C, A, res, a)
    st = settings_dict(C, A, a)
    with open(os.path.join(a.out, "cohort.settings.json"), "w") as fh:
        json.dump(st, fh, indent=2, default=str)
    if not a.no_release_tables:
        for name, tab in RELEASE_TABLES.items():
            df = res.get(name)
            if isinstance(df, pd.DataFrame) and len(df):
                d = df.copy()
                for c in d.columns:
                    if d[c].dtype == object:
                        d[c] = d[c].astype(str)
                write_table(d, a.data, tab)
        m = read_manifest(a.data)
        m["cohort"] = st
        with open(os.path.join(a.data, "manifest.json"), "w") as fh:
            json.dump(m, fh, indent=2, default=str)
    if not a.no_report:
        from .report import write_report
        write_report(C, A, res, st, a.out)
    C.note("done in %.0f s -> %s" % (time.time() - t0, a.out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
