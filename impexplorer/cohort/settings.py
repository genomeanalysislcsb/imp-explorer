"""Settings of the cohort analysis: command line and an optional settings file (--config).

The settings file is study-wide, not tied to sRNA: 'key: value' lines (e.g. written by --write-config) or bash
assignments KEY="value" / KEY="${KEY-default}". Keys are the option names (group_col, block, model, ...) with or
without the prefix COHORT_, in any case, so the COHORT_* block of an sRNA-IMP config works unchanged and one file can
hold the design of a study for all its data types."""
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
    "COHORT_QC_DROP": "qc_drop", "COHORT_DICTIONARY": "dictionary",
}

ALIASES = {"match_col": "block", "group": "group_col", "groups": "group_order", "title": "name"}
CONFIG_KEYS = ("metadata", "dictionary", "id_col", "group_col", "group_order", "reference", "block", "covariates", "model",
               "model_notes", "permutations", "sig_q", "outlier_z", "qc_drop", "term_names", "name", "layers",
               "prevalence", "min_mg_pairs", "min_mt_pairs", "nmds", "nmds_starts")
IGNORED = {"dir", "min_mirna_reads"}       # sRNA-specific COHORT_* keys without meaning here

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


def norm_key(k):
    k = k.strip().lower().replace("-", "_")
    k = k[len("cohort_"):] if k.startswith("cohort_") else k
    return ALIASES.get(k, k)


def read_config(path):
    """settings file -> {option: value} for the known options; unknown keys are reported, sRNA-only keys ignored"""
    raw = read_shell_config(path)
    for ln in open(path, encoding="utf-8", errors="replace"):
        m = re.match(r"\s*([A-Za-z_][\w-]*)\s*:\s*(.*?)\s*$", ln)
        if not m or ln.lstrip().startswith("#"):
            continue
        val = m.group(2)
        if val.startswith(("'", '"')):
            q = val[0]
            val = val[1:val.index(q, 1)] if q in val[1:] else val[1:]
        else:
            val = re.sub(r"\s+#.*$", "", val)
        raw[m.group(1)] = val
    out, unknown = {}, []
    for k, v in raw.items():
        nk = norm_key(k)
        if nk in CONFIG_KEYS:
            if v != "" and nk not in out:
                out[nk] = v
        elif k.upper().startswith("COHORT_") or ":" in k or k.islower():
            if nk not in IGNORED:
                unknown.append(k)
    return out, unknown


TEMPLATE = """# imp-explorer cohort settings (imp-explorer cohort --config this_file)
# key: value; empty = default. Keys also work as COHORT_<KEY>="..." (bash), e.g. in an sRNA-IMP config.
name:                     # report title  [study name of the release]
metadata:                 # extra metadata TSV/CSV merged on the sample id
dictionary:               # metadata.dictionary.tsv: column types (numeric / categorical), missing codes
id_col:                   # id column of the metadata  [first column]
group_col: group          # group column
group_order:              # e.g. control,case_a,case_b
reference:                # reference group  [first group]
block:                    # block column, e.g. site: groups compared within blocks
covariates:               # covariates checked for balance  [all metadata columns]
# models: comma-separated, terms joined by +, first term = term of interest; 'terms | col=v1/v2' = subset
#   e.g. group+site, group+site+sex+age, dose+sex+age | group=case_a    [empty: group(+block)]
model: {model}
model_notes:              # question per model, e.g. M2=adjusted for sex and age; M3=dose response
permutations: 999
sig_q: 0.1
outlier_z: 4
qc_drop: FAIL             # FAIL | WARN | none
prevalence: 0.1
min_mg_pairs: 1000000
min_mt_pairs: 1000000
term_names:               # term<TAB>name TSVs or KEGG flat files for readable labels
nmds: kraken.genus.mg,kegg_pathway.mg
nmds_starts: 5
layers: {layers}
"""


def parser():
    ap = argparse.ArgumentParser(
        prog="imp-explorer cohort",
        description="Cohort analysis of an IMP3 release: QC and covariates, taxa and diversity, group comparisons "
                    "within blocks (e.g. site), covariate-adjusted models (OLS on CLR, HC3), MAGs, functions and "
                    "metatranscriptome / metagenome ratios. Writes TSVs, a PDF + HTML report and cohort tables "
                    "for the app.")
    ap.add_argument("--data", default="", help="release directory (imp-explorer export)")
    ap.add_argument("--out", default="", help="output directory  [<release>/cohort]")
    ap.add_argument("--config", "--srna-config", dest="config", default="",
                    help="settings file of the study: 'key: value' lines (see --write-config) or bash COHORT_* "
                         "assignments, e.g. the cohort block of an sRNA-IMP config; command-line options override it")
    ap.add_argument("--write-config", default="", metavar="FILE",
                    help="write a commented settings template (with the options given) and exit")
    ap.add_argument("--metadata", default=None, help="extra metadata TSV/CSV merged on the sample id "
                    "(columns override those of the release)")
    ap.add_argument("--dictionary", default=None, help="metadata dictionary (imp-explorer metadata template): its "
                    "types decide numeric / categorical columns (e.g. 0/1 or 1-7 codes), its missing codes become NA")
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
    ap.add_argument("--layers", default=None, help="feature sets to analyse [%s]" % DEFAULT_LAYERS)
    ap.add_argument("--permutations", type=int, default=None, help="PERMANOVA permutations  [999]")
    ap.add_argument("--sig-q", type=float, default=None, help="FDR threshold for highlights  [0.1]")
    ap.add_argument("--outlier-z", type=float, default=None, help="robust z of a QC metric that gives WARN  [4]")
    ap.add_argument("--min-mg-pairs", type=float, default=None, help="FAIL below this many trimmed mg read pairs")
    ap.add_argument("--min-mt-pairs", type=float, default=None, help="mt layers: sample left out below this many "
                    "trimmed mt read pairs (WARN)")
    ap.add_argument("--qc-drop", default=None, choices=("FAIL", "WARN", "none"), help="samples left out  [FAIL]")
    ap.add_argument("--prevalence", type=float, default=None, help="features with >= 10 reads in at least this share "
                    "of samples are tested  [0.1]")
    ap.add_argument("--nmds", default=None, help="layers with an NMDS ('all' = every layer, '' = none); sRNA-IMP names (genus, "
                    "fn_kegg_pathway) mean the metagenome layers (kraken.genus.mg, kegg_pathway.mg)")
    ap.add_argument("--nmds-starts", type=int, default=None, help="random starts per NMDS (+ 1 PCoA start)  [5]")
    ap.add_argument("--term-names", default=None, help="TSV(s) term<TAB>name or KEGG flat files for readable labels")
    ap.add_argument("--name", default=None, help="title of the report  [study name]")
    ap.add_argument("--no-release-tables", action="store_true", help="do not write cohort_* tables into the release")
    ap.add_argument("--no-report", action="store_true", help="TSVs only")
    ap.add_argument("--threads", type=int, default=1, help="processes for the beta-diversity step (PCoA, NMDS, "
                    "PERMANOVA per layer); same results as 1  [1]")
    return ap


DEFAULTS = {"metadata": "", "dictionary": "", "id_col": "", "group_col": "group", "group_order": "", "reference": "", "block": "",
            "covariates": "", "model": "", "model_notes": "", "permutations": 999, "sig_q": 0.1, "outlier_z": 4.0,
            "qc_drop": "FAIL", "term_names": "", "name": "", "layers": DEFAULT_LAYERS, "prevalence": 0.1,
            "min_mg_pairs": 1e6, "min_mt_pairs": 1e6, "nmds": "kraken.genus.mg,kegg_pathway.mg", "nmds_starts": 5}
TYPES = {"permutations": int, "nmds_starts": int, "sig_q": float, "outlier_z": float, "prevalence": float, "min_mg_pairs": float,
         "min_mt_pairs": float}


RANKS = ("phylum", "class", "order", "family", "genus", "species")


def nmds_layers(value):
    """NMDS layers; names of the sRNA-IMP config (COHORT_NMDS: genus, fn_kegg_pathway, mirna, ...) are translated to
    the metagenome layers here (kraken.genus.mg, kegg_pathway.mg); sRNA-only sets (mirna, rfam_*, ...) are left out"""
    out = []
    for k in [x.strip() for x in str(value or "").split(",") if x.strip()]:
        if k.lower() in ("all", "*"):
            out.append("all")
        elif "." in k:
            out.append(k)
        elif k in RANKS:
            out.append("kraken.%s.mg" % k)
        elif k.startswith("fn_"):
            out.append("%s.mg" % k[3:])
    return ",".join(dict.fromkeys(out))


def parse(argv=None):
    ap = parser()
    a = ap.parse_args(argv)
    cfg, unknown = read_config(a.config) if a.config else ({}, [])
    if unknown:
        import sys
        sys.stderr.write("[cohort] %s: unknown settings ignored: %s\n" % (a.config, ", ".join(sorted(unknown))))
    for k, d in DEFAULTS.items():
        v = getattr(a, k)
        if v is None:
            v = cfg.get(k, d)
            v = TYPES[k](v) if k in TYPES and v != "" else v
        setattr(a, k, v)
    a.nmds = nmds_layers(a.nmds)
    if a.qc_drop not in ("FAIL", "WARN", "none"):
        ap.error("qc_drop must be FAIL, WARN or none (got %r)" % a.qc_drop)
    a.from_config = sorted(cfg)
    if a.write_config:
        write_template(a, a.write_config)
        return a
    if not a.data:
        ap.error("--data is required")
    a.out = a.out or os.path.join(a.data, "cohort")
    return a


def write_template(a, path):
    """settings template with the current values (defaults, --config, command line)"""
    txt = TEMPLATE.format(model=a.model or "", layers=a.layers)
    lines = []
    for ln in txt.splitlines():
        m = re.match(r"([a-z_]+):(\s*)(.*)$", ln)
        if m and m.group(1) in CONFIG_KEYS and m.group(1) not in ("model", "layers"):
            v = getattr(a, m.group(1))
            v = "" if v is None else (("%d" % v if float(v).is_integer() else "%g" % v) if isinstance(v, float)
                                       else str(v))
            comment = m.group(3).split("#", 1)[1] if "#" in m.group(3) else ""
            ln = "%s: %s" % (m.group(1), v)
            ln = ln.ljust(26) + ("# " + comment.strip() if comment else "")
        lines.append(ln.rstrip())
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
