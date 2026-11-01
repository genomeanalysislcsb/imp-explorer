"""imp-explorer cohort: cohort analysis of an IMP3 release (QC and covariates, taxa and diversity, group comparisons
within blocks, covariate-adjusted models, MAGs, functions, metatranscriptome / metagenome ratios).

Writes <out>/tables/*.tsv, <out>/stats/*.tsv, <out>/cohort.report.pdf + .html and cohort_* tables into the release
for the app's Cohort page. Same statistics and settings as the sRNA-IMP cohort analysis (--srna-config).
"""
import json
import os
import sys
import time

import pandas as pd

from .. import __version__
from ..release import read_manifest, write_table
from . import settings
from .data import Cohort, METRIC_LABEL
from .engine import Analysis
from .stats import load_term_names

RELEASE_TABLES = {   # cohort result -> release table read by the app
    "models": "cohort_models", "models_overview": "cohort_models_overview", "feature_tests": "cohort_tests",
    "metrics_tests": "cohort_metric_tests", "permanova": "cohort_permanova", "dispersion": "cohort_dispersion",
    "ordination": "cohort_ordination", "sample_metrics": "cohort_samples", "covariates": "cohort_covariates",
    "qc_flags": "cohort_qc", "mag_prevalence": "cohort_mag_prevalence"}


def write_outputs(C, A, res, a):
    tdir, sdir = os.path.join(a.out, "tables"), os.path.join(a.out, "stats")
    os.makedirs(tdir, exist_ok=True)
    os.makedirs(sdir, exist_ok=True)
    meta = C.samples.copy()
    meta.drop(columns=[c for c in ("imp_dir",) if c in meta.columns]).to_csv(
        os.path.join(tdir, "cohort.metadata.tsv"), sep="\t")
    for key, L in C.layers.items():
        L.M.T.to_csv(os.path.join(tdir, "cohort.%s.%s.tsv" % (key, "counts" if L.kind == "counts" else "log2ratio")),
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
            "samples": len(C.samples), "kept": len(C.kept), "excluded": C.excluded,
            "from_srna_config": a.from_srna, "notes": C.log}


def main(argv=None):
    a = settings.parse(argv)
    t0 = time.time()
    if a.term_names:
        load_term_names(a.term_names)
    C = Cohort(a)
    C.note("%d samples, %d kept; groups %s (reference %s); block %s" % (
        len(C.samples), len(C.kept), ",".join(C.groups) or "–", C.ref or "–", C.block or "–"))
    C.load_layers(a.layers)
    C.note("layers: " + ", ".join("%s (%d x %d)" % (k, L.M.shape[0], L.M.shape[1]) for k, L in C.layers.items()))
    A = Analysis(C)
    res = A.run()
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
