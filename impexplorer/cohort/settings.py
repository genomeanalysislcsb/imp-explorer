"""Settings of the cohort analysis: command line, optionally read from the cohort block of an sRNA-IMP config
(srna_pipeline.config.sh: COHORT_GROUP_COL, COHORT_MODEL, ...), so that IMP3 and sRNA use the same design."""
import argparse
import os
import re
import shlex

# sRNA-IMP config key -> option of this command
SRNA_KEYS = {
    "COHORT_METADATA": "metadata", "COHORT_ID_COL": "id_col", "COHORT_GROUP_COL": "group_col",
    "COHORT_GROUP_ORDER": "group_order", "COHORT_COVARIATES": "covariates", "COHORT_PERMUTATIONS": "permutations",
    "COHORT_BLOCK": "block", "COHORT_MATCH_COL": "block", "COHORT_REFERENCE": "reference",
    "COHORT_OUTLIER_Z": "outlier_z", "COHORT_SIG_Q": "sig_q", "COHORT_TERM_NAMES": "term_names",
    "COHORT_NAME": "name", "COHORT_MODEL": "model", "COHORT_MODEL_NOTES": "model_notes",
    "COHORT_QC_DROP": "qc_drop",
}

DEFAULT_LAYERS = ("kraken.genus.mg,kraken.genus.mt,kraken.species.mg,kraken.species.mt,"
                  "motus.species.mg,motus.species.mt,mag.genus.mg,mag.genus.mt,"
                  "kegg_module.mg,kegg_module.mt,kegg_pathway.mg,kegg_pathway.mt,cazy.mg,cazy.mt,"
                  "kegg_ko.mg,kegg_ko.mt,ratio.mag.genus,ratio.kegg_module,ratio.kegg_pathway")


def read_shell_config(path):
    """KEY="value" / KEY="${KEY-default}" / KEY="${KEY:-default}" assignments of a bash config -> dict.
    Values that refer to other variables (${OTHER}) are expanded from the keys read so far; nothing is executed."""
    out = {}
    for ln in open(path, encoding="utf-8", errors="replace"):
        m = re.match(r"\s*(?:export\s+)?([A-Z_][A-Z0-9_]*)=(.*)$", ln)
        if not m:
            continue
        key, rest = m.group(1), m.group(2).strip()
        try:
            toks = shlex.split(rest, comments=True, posix=True)
        except ValueError:
            continue
        val = toks[0] if toks else ""

        def expand(v):
            def rep(mm):
                name, op, default = mm.group(1), mm.group(2), mm.group(3)
                cur = out.get(name, os.environ.get(name))
                if op in (":-", "-") and (cur is None or (op == ":-" and cur == "")):
                    return expand(default or "")
                return cur or ""
            v = re.sub(r"\$\{([A-Za-z_][A-Za-z0-9_]*)(:?-)?([^{}]*(?:\{[^{}]*\}[^{}]*)*)\}", rep, v)
            return re.sub(r"\$([A-Za-z_][A-Za-z0-9_]*)", lambda mm: out.get(mm.group(1), ""), v)
        out[key] = expand(val)
    return out


def parser():
    ap = argparse.ArgumentParser(
        prog="imp-explorer cohort",
        description="Cohort analysis of an IMP3 release: QC and covariates, taxa and diversity, group comparisons "
                    "within blocks (e.g. site), covariate-adjusted models (OLS on CLR, HC3), MAGs, functions and "
                    "metatranscriptome / metagenome ratios. Writes TSVs, a PDF + HTML report and cohort tables "
                    "for the app.")
    ap.add_argument("--data", required=True, help="release directory (imp-explorer export)")
    ap.add_argument("--out", default="", help="output directory  [<release>/cohort]")
    ap.add_argument("--srna-config", default="", help="read the COHORT_* settings of an sRNA-IMP config "
                    "(srna_pipeline.config.sh); options given here override them")
    ap.add_argument("--metadata", default=None, help="extra metadata TSV/CSV merged on the sample id "
                    "(columns override those of the release)")
    ap.add_argument("--id-col", default=None, help="sample id column of --metadata  [first column]")
    ap.add_argument("--group-col", default=None, help="group column  [group]")
    ap.add_argument("--group-order", default=None, help="comma-separated group order, e.g. control,case_a,case_b")
    ap.add_argument("--reference", default=None, help="reference group  [first group]")
    ap.add_argument("--block", default=None, help="block column, e.g. site: groups are compared within blocks, "
                    "PERMANOVA permutes within blocks")
    ap.add_argument("--covariates", default=None, help="comma-separated covariate columns  [all metadata columns]")
    ap.add_argument("--model", default=None, help="adjusted models, comma-separated; terms joined by '+'; first term "
                    "= term of interest; 'terms | column=v1/v2' fits on a subset, e.g. "
                    "'group+site, group+site+sex+age, dose+sex+age | group=case_a'  [group(+block)]")
    ap.add_argument("--model-notes", default=None, help="question per model, ';'-separated: 'M3=dose response'")
    ap.add_argument("--layers", default=DEFAULT_LAYERS, help="feature sets to analyse [%(default)s]")
    ap.add_argument("--permutations", type=int, default=None, help="PERMANOVA permutations  [999]")
    ap.add_argument("--sig-q", type=float, default=None, help="FDR threshold for highlights  [0.1]")
    ap.add_argument("--outlier-z", type=float, default=None, help="robust z of a QC metric that gives WARN  [4]")
    ap.add_argument("--min-mg-pairs", type=float, default=1e6, help="FAIL below this many trimmed mg read pairs")
    ap.add_argument("--min-mt-pairs", type=float, default=1e6, help="mt layers: sample left out below this many "
                    "trimmed mt read pairs (WARN)")
    ap.add_argument("--qc-drop", default=None, choices=("FAIL", "WARN", "none"), help="samples left out  [FAIL]")
    ap.add_argument("--prevalence", type=float, default=0.1, help="features with >= 10 reads in at least this share "
                    "of samples are tested  [0.1]")
    ap.add_argument("--nmds", default="kraken.genus.mg,kegg_pathway.mg", help="layers with an NMDS ('' = none)")
    ap.add_argument("--term-names", default=None, help="TSV(s) term<TAB>name or KEGG flat files for readable labels")
    ap.add_argument("--name", default=None, help="title of the report  [study name]")
    ap.add_argument("--no-release-tables", action="store_true", help="do not write cohort_* tables into the release")
    ap.add_argument("--no-report", action="store_true", help="TSVs only")
    ap.add_argument("--threads", type=int, default=1, help="(reserved)")
    return ap


DEFAULTS = {"metadata": "", "id_col": "", "group_col": "group", "group_order": "", "reference": "", "block": "",
            "covariates": "", "model": "", "model_notes": "", "permutations": 999, "sig_q": 0.1, "outlier_z": 4.0,
            "qc_drop": "FAIL", "term_names": "", "name": ""}
TYPES = {"permutations": int, "sig_q": float, "outlier_z": float}


def parse(argv=None):
    a = parser().parse_args(argv)
    cfg = {}
    if a.srna_config:
        raw = read_shell_config(a.srna_config)
        for k, opt in SRNA_KEYS.items():
            if k in raw and raw[k] != "" and opt not in cfg:
                cfg[opt] = raw[k]
    for k, d in DEFAULTS.items():
        v = getattr(a, k)
        if v is None:
            v = cfg.get(k, d)
            v = TYPES[k](v) if k in TYPES and v != "" else v
        setattr(a, k, v)
    a.from_srna = sorted(cfg) if a.srna_config else []
    a.out = a.out or os.path.join(a.data, "cohort")
    return a
