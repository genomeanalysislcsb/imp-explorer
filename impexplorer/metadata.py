"""imp-explorer metadata: a template for the sample metadata of a study, and a check of a metadata table before it
goes into a release or a cohort analysis.

  imp-explorer metadata template --out DIR [--data RELEASE]
      metadata.dictionary.tsv  one row per column: type, unit, allowed values or range, missing codes, description
      metadata.template.tsv    the columns of the dictionary, one row per sample of the release (if --data)
  imp-explorer metadata check --metadata TABLE [--dictionary DICT] [--data RELEASE] [--write-clean OUT]
      checks ids (duplicates, matching the release, with --id-map / --id-regex / --id-case as for export),
      columns against the dictionary (types, units, allowed values, ranges), missing values, group sizes;
      --write-clean writes the table with the ids translated, missing codes emptied, decimal commas and the
      spelling of categories fixed. Exit code 1 when there are errors.

The dictionary is also read by 'imp-explorer cohort --dictionary': its types decide whether a column is numeric or
categorical in the models (a 0/1 or 1-7 coded column is otherwise numeric) and its missing codes become NA.

Dictionary columns: column, type (id | categorical | binary | integer | numeric | date | text), unit, allowed
(categories separated by |, or a range min..max), missing (extra missing codes separated by |, e.g. -9|999),
category (sample | clinical | medication | diet | lifestyle | technical | ...), description."""
import argparse
import os
import re
import sys

import pandas as pd

MISSING = {"", "na", "n/a", "nan", "none", "null", "-", "?", "nd", "n.d.", "not available", "unknown"}
TYPES = ("id", "categorical", "binary", "integer", "numeric", "date", "text")
DICT_COLS = ("column", "type", "unit", "allowed", "missing", "category", "description")

# example dictionary: generic columns of a human gut microbiome cohort; edit, delete and add rows for a study
TEMPLATE_ROWS = [
    ("sample_id", "id", "", "", "", "sample", "sample id as in the release (or translated with --id-map / --id-regex)"),
    ("group", "categorical", "", "control|case", "", "sample", "study group; the first is the reference"),
    ("site", "categorical", "", "", "", "sample", "recruitment site / centre (block of the cohort analysis)"),
    ("sex", "categorical", "", "F|M", "", "clinical", "sex"),
    ("age", "numeric", "years", "18..100", "", "clinical", "age at sampling"),
    ("bmi", "numeric", "kg/m2", "12..70", "", "clinical", "body mass index"),
    ("age_at_onset", "numeric", "years", "0..100", "", "clinical", "age at disease onset (cases)"),
    ("disease_duration", "numeric", "years", "0..80", "", "clinical", "years since diagnosis (cases)"),
    ("constipation", "binary", "", "yes|no", "", "clinical", "constipation (e.g. Rome IV criteria)"),
    ("bristol_stool_scale", "integer", "", "1..7", "", "clinical", "Bristol stool form of the sample"),
    ("antibiotics_3m", "binary", "", "yes|no", "", "medication", "antibiotics in the 3 months before sampling"),
    ("ppi", "binary", "", "yes|no", "", "medication", "proton pump inhibitors"),
    ("levodopa", "binary", "", "yes|no", "", "medication", "levodopa treatment"),
    ("ledd", "numeric", "mg/day", "0..5000", "", "medication", "levodopa equivalent daily dose"),
    ("methotrexate", "binary", "", "yes|no", "", "medication", "methotrexate treatment"),
    ("smoking", "categorical", "", "never|former|current", "", "lifestyle", "smoking status"),
    ("alcohol_units_week", "numeric", "units/week", "0..200", "", "lifestyle", "alcohol consumption"),
    ("diet_pattern", "categorical", "", "omnivore|pescatarian|vegetarian|vegan", "", "diet", "self-reported diet"),
    ("energy_kcal_day", "numeric", "kcal/day", "300..8000", "", "diet", "energy intake (food frequency questionnaire)"),
    ("fibre_g_day", "numeric", "g/day", "0..200", "", "diet", "dietary fibre intake"),
    ("protein_pct_energy", "numeric", "% energy", "0..100", "", "diet", "protein share of energy intake"),
    ("fat_pct_energy", "numeric", "% energy", "0..100", "", "diet", "fat share of energy intake"),
    ("carbohydrate_pct_energy", "numeric", "% energy", "0..100", "", "diet", "carbohydrate share of energy intake"),
    ("coffee_cups_day", "numeric", "cups/day", "0..30", "", "diet", "coffee consumption"),
    ("sampling_date", "date", "", "", "", "technical", "date of stool sampling (YYYY-MM-DD)"),
]


def say(msg):
    sys.stderr.write("[metadata] %s\n" % msg)


# ------------------------------------------------------------------------------------------------- reading
def read_table_any(path):
    """TSV or CSV (by the first line), all as text, comment rows (#...) dropped"""
    head = open(path, encoding="utf-8-sig", errors="replace").read(4096)
    sep = "\t" if head.count("\t") >= head.count(",") else ","
    d = pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False, encoding="utf-8-sig")
    d.columns = [str(c).strip() for c in d.columns]
    if len(d.columns):
        d = d[~d.iloc[:, 0].astype(str).str.startswith("#")]
    return d.reset_index(drop=True)


def read_dictionary(path):
    d = read_table_any(path)
    miss = [c for c in ("column", "type") if c not in d.columns]
    if miss:
        raise ValueError("dictionary %s lacks the columns %s" % (path, ", ".join(miss)))
    for c in DICT_COLS:
        if c not in d.columns:
            d[c] = ""
    d["type"] = d["type"].str.strip().str.lower()
    bad = sorted(set(d["type"]) - set(TYPES))
    if bad:
        raise ValueError("dictionary %s: unknown types %s (use %s)" % (path, ", ".join(bad), ", ".join(TYPES)))
    return {r["column"].strip(): {k: str(r[k]).strip() for k in DICT_COLS} for _, r in d.iterrows()
            if r["column"].strip()}


def parse_range(s):
    m = re.match(r"^\s*(-?[\d.]+)?\s*\.\.\s*(-?[\d.]+)?\s*$", s or "")
    if not m:
        return None
    return (float(m.group(1)) if m.group(1) else None, float(m.group(2)) if m.group(2) else None)


def to_number(v):
    """value -> float or None; accepts a decimal comma (12,5) and surrounding spaces"""
    t = str(v).strip().replace("\u2212", "-")
    if re.match(r"^-?\d+,\d+$", t):
        t = t.replace(",", ".")
    try:
        return float(t)
    except ValueError:
        return None


# ------------------------------------------------------------------------------------------------- normalising
def normalize(df, dictionary=None, id_col=None):
    """apply the dictionary (or, without one, only the missing codes) to a metadata table (all text).

    Returns (clean table, issues, kinds): issues are (level, column, message, examples); kinds maps a column to
    numeric / categorical / text for the cohort analysis (only columns of the dictionary)."""
    dictionary = dictionary or {}
    out = df.copy()
    issues, kinds = [], {}
    for col in out.columns:
        spec = dictionary.get(col)
        codes = set(MISSING)
        if spec and spec.get("missing"):
            codes |= {x.strip().lower() for x in spec["missing"].split("|") if x.strip()}
        vals = out[col].astype(str).map(lambda v: v.strip())
        is_na = vals.str.lower().isin(codes)
        vals = vals.where(~is_na, "")
        if not spec or col == id_col:
            out[col] = vals
            continue
        t = spec["type"]
        present = vals[vals != ""]
        if t in ("numeric", "integer"):
            nums = pd.Series([to_number(v) for v in present], index=present.index, dtype=float)
            bad = present[nums.isna()]
            if len(bad):
                issues.append(("ERROR", col, "%d values are not numbers (set to missing in the clean table)"
                               % len(bad), sorted(set(bad))[:5]))
            comma = present[present.str.match(r"^-?\d+,\d+$")]
            if len(comma):
                issues.append(("WARN", col, "%d values with a decimal comma (read as a dot)" % len(comma),
                               sorted(set(comma))[:3]))
            ok = nums.dropna()
            if t == "integer":
                frac = ok[(ok % 1) != 0]
                if len(frac):
                    issues.append(("ERROR", col, "%d values are not whole numbers" % len(frac),
                                   [repr(x) for x in sorted(set(frac))[:5]]))
            rng = parse_range(spec.get("allowed"))
            if rng:
                lo, hi = rng
                outside = ok[((ok < lo) if lo is not None else False) | ((ok > hi) if hi is not None else False)]
                if len(outside):
                    issues.append(("WARN", col, "%d values outside %s %s (wrong unit?)"
                                   % (len(outside), spec["allowed"], spec.get("unit", "")),
                                   [repr(x) for x in sorted(set(outside))[:5]]))
            num = vals.map(lambda v: to_number(v) if v != "" else None)
            out[col] = num.map(lambda x: "" if x is None or pd.isna(x) else
                               (str(int(x)) if float(x).is_integer() else repr(float(x))))
            kinds[col] = "numeric"
        elif t in ("categorical", "binary"):
            allowed = [x.strip() for x in (spec.get("allowed") or "").split("|") if x.strip()]
            if allowed:
                canon = {a.lower(): a for a in allowed}
                fixed = present.map(lambda v: canon.get(v.lower()))
                case = present[fixed.notna() & (fixed != present)]
                if len(case):
                    issues.append(("WARN", col, "%d values differ from the dictionary only in case (fixed)"
                                   % len(case), sorted(set(case))[:5]))
                unknown = present[fixed.isna()]
                if len(unknown):
                    issues.append(("ERROR", col, "%d values not in %s" % (len(unknown), "|".join(allowed)),
                                   sorted(set(unknown))[:5]))
                out[col] = vals.map(lambda v: canon.get(v.lower(), v) if v else "")
            else:
                out[col] = vals
            levels = set(out[col][out[col] != ""])
            if t == "binary" and len(levels) > 2 and not allowed:
                issues.append(("ERROR", col, "binary column with %d different values" % len(levels),
                               sorted(levels)[:5]))
            kinds[col] = "categorical"
        elif t == "date":
            dt = pd.to_datetime(present, errors="coerce", dayfirst=False)
            bad = present[dt.isna()]
            if len(bad):
                issues.append(("ERROR", col, "%d values are not dates (YYYY-MM-DD)" % len(bad), sorted(set(bad))[:5]))
            out[col] = vals.map(lambda v: "" if v == "" else
                                (lambda x: "" if pd.isna(x) else x.strftime("%Y-%m-%d"))(pd.to_datetime(v, errors="coerce")))
            kinds[col] = "text"
        else:
            out[col] = vals
            kinds[col] = "text"
    return out, issues, kinds


# ------------------------------------------------------------------------------------------------- template
def cmd_template(a):
    os.makedirs(a.out, exist_ok=True)
    rows = list(TEMPLATE_ROWS)
    samples = None
    if a.data:
        from .release import read_table
        samples = read_table(a.data, "samples")
        groups = sorted({str(x) for x in samples.get("group", pd.Series(dtype=str)).dropna() if str(x).strip()})
        if groups:
            rows[1] = ("group",) + rows[1][1:2] + ("", "|".join(groups)) + rows[1][4:]
    dct = pd.DataFrame(rows, columns=DICT_COLS)
    p_dict = os.path.join(a.out, "metadata.dictionary.tsv")
    dct.to_csv(p_dict, sep="\t", index=False)
    cols = list(dct["column"])
    tmpl = pd.DataFrame(columns=cols)
    if samples is not None and len(samples):
        tmpl = pd.DataFrame({c: "" for c in cols}, index=range(len(samples)))
        tmpl["sample_id"] = samples["sample_id"].astype(str).values
        for c in cols[1:]:
            if c in samples.columns:
                tmpl[c] = samples[c].fillna("").astype(str).values
    p_tmpl = os.path.join(a.out, "metadata.template.tsv")
    tmpl.to_csv(p_tmpl, sep="\t", index=False)
    say("wrote %s (%d columns: edit types, units, allowed values) and %s (%d samples)"
        % (p_dict, len(dct), p_tmpl, len(tmpl)))
    return 0


# ------------------------------------------------------------------------------------------------- check
def cmd_check(a):
    md = read_table_any(a.metadata)
    dictionary = read_dictionary(a.dictionary) if a.dictionary else {}
    issues = []
    idc = a.id_col if a.id_col else ("sample_id" if "sample_id" in md.columns else md.columns[0])
    if idc not in md.columns:
        sys.exit("metadata: id column %s not in %s" % (idc, a.metadata))
    ids = md[idc].astype(str).str.strip()
    if (ids == "").any():
        issues.append(("ERROR", idc, "%d rows without a sample id" % int((ids == "").sum()), []))
    if a.id_map or a.id_regex or a.id_case != "keep":
        from .export import derive_ids, read_id_map
        try:
            ids = pd.Series(derive_ids(list(ids), read_id_map(a.id_map) if a.id_map else None, a.id_regex,
                                       a.id_case), index=md.index)
        except ValueError as e:
            issues.append(("ERROR", idc, str(e), []))
    md[idc] = ids
    dup = ids[ids.duplicated(keep=False) & (ids != "")]
    if len(dup):
        issues.append(("ERROR", idc, "%d rows with a duplicated sample id" % len(dup), sorted(set(dup))[:5]))
    # columns vs dictionary
    if dictionary:
        missing_cols = [c for c in dictionary if c not in md.columns and dictionary[c]["type"] != "id"]
        if missing_cols:
            issues.append(("WARN", "-", "%d dictionary columns not in the table" % len(missing_cols), missing_cols[:8]))
        undocumented = [c for c in md.columns if c not in dictionary and c != idc]
        if undocumented:
            issues.append(("INFO", "-", "%d columns not in the dictionary (type guessed in the cohort analysis)"
                           % len(undocumented), undocumented[:8]))
    clean, more, kinds = normalize(md, dictionary, idc)
    issues += more
    # without a dictionary: suspicious columns
    empty = [c for c in clean.columns if c != idc and (clean[c] == "").all()]
    if empty:
        issues.append(("INFO", "-", "%d columns without any value (not filled in yet?)" % len(empty), empty[:8]))
    for col in clean.columns:
        if col == idc or col in empty:
            continue
        v = clean[col][clean[col] != ""]
        n_missing = int((clean[col] == "").sum())
        if len(clean) and n_missing / len(clean) > a.max_missing:
            issues.append(("WARN", col, "%.0f %% missing" % (100.0 * n_missing / len(clean)), []))
        if len(v) and v.nunique() == 1:
            issues.append(("INFO", col, "constant (%s)" % v.iloc[0], []))
        if col in dictionary or not len(v):
            continue
        nums = pd.Series([to_number(x) for x in v], index=v.index, dtype=float)
        if 0 < nums.isna().sum() <= max(2, 0.1 * len(v)):
            issues.append(("WARN", col, "mostly numbers, but %d values are not (e.g. '<5', 'n.d.')"
                           % int(nums.isna().sum()), sorted(set(v[nums.isna()]))[:5]))
        if nums.isna().all():
            low = v.str.lower()
            if low.nunique() < v.nunique():
                issues.append(("WARN", col, "categories that differ only in case", sorted(set(v))[:6]))
        elif nums.notna().all() and v.nunique() <= 7:
            issues.append(("INFO", col, "numeric codes (%s): numeric in the models unless the dictionary says "
                           "categorical" % ", ".join(sorted(set(v))[:7]), []))
    # release
    if a.data:
        from .release import read_table
        from .cohort.data import norm_id
        rel = read_table(a.data, "samples")
        rel_ids = {norm_id(x) for x in rel["sample_id"].astype(str)}
        mine = {norm_id(x) for x in ids if x}
        no_rel = sorted(x for x in ids if x and norm_id(x) not in rel_ids)
        no_md = sorted(x for x in rel["sample_id"].astype(str) if norm_id(x) not in mine)
        say("release %s: %d of %d release samples have metadata" % (a.data, len(rel_ids & mine), len(rel_ids)))
        if no_rel:
            issues.append(("WARN", idc, "%d metadata ids not in the release (check --id-map / --id-regex / --id-case)"
                           % len(no_rel), no_rel[:8]))
        if no_md:
            issues.append(("WARN", idc, "%d release samples without metadata" % len(no_md), no_md[:8]))
    # groups
    gc = a.group_col
    if gc in clean.columns:
        gt = clean[clean[gc] != ""].groupby([gc] + ([a.block] if a.block in clean.columns else [])).size()
        say("samples per group%s:\n%s" % (" and %s" % a.block if a.block in clean.columns else "",
                                           gt.to_string()))
        small = gt[gt < a.min_group]
        if len(small):
            issues.append(("WARN", gc, "groups with fewer than %d samples" % a.min_group,
                           ["%s: %d" % (k if isinstance(k, str) else "/".join(map(str, k)), n)
                            for k, n in small.items()][:6]))
    # summary per column
    summ = []
    for col in clean.columns:
        v = clean[col][clean[col] != ""]
        kind = kinds.get(col) or ("id" if col == idc else
                                  ("numeric" if len(v) and v.map(to_number).notna().all() and v.nunique() > 2
                                   else "categorical"))
        if kind == "numeric" and len(v):
            x = v.map(to_number).astype(float)
            desc = "median %.4g [%.4g–%.4g], range %.4g–%.4g" % (x.median(), x.quantile(.25), x.quantile(.75),
                                                                  x.min(), x.max())
        elif kind in ("categorical",) and len(v):
            desc = ", ".join("%s %d" % (k, n) for k, n in v.value_counts().head(8).items())
        else:
            desc = "%d distinct" % v.nunique()
        unit = (dictionary.get(col) or {}).get("unit", "")
        summ.append((col, kind + (" (%s)" % unit if unit else ""), len(v), len(clean) - len(v), desc))
    print("%-28s %-26s %6s %8s  %s" % ("column", "type", "values", "missing", "summary"))
    for r in summ:
        print("%-28s %-26s %6d %8d  %s" % (r[0][:28], r[1][:26], r[2], r[3], r[4][:90]))
    print()
    order = {"ERROR": 0, "WARN": 1, "INFO": 2}
    for lvl, col, msg, ex in sorted(issues, key=lambda x: order[x[0]]):
        print("%-5s %-24s %s%s" % (lvl, col, msg, ("  e.g. " + ", ".join(map(str, ex))) if ex else ""))
    n_err = sum(1 for x in issues if x[0] == "ERROR")
    n_warn = sum(1 for x in issues if x[0] == "WARN")
    print("\n%d errors, %d warnings, %d rows, %d columns" % (n_err, n_warn, len(clean), len(clean.columns)))
    if a.report:
        pd.DataFrame([(l, c, m, "; ".join(map(str, e))) for l, c, m, e in issues],
                     columns=["level", "column", "message", "examples"]).to_csv(a.report, sep="\t", index=False)
    if a.write_clean:
        clean.rename(columns={idc: "sample_id"}).to_csv(a.write_clean, sep="\t", index=False)
        say("clean table: %s (ids%s, missing codes empty, numbers with a dot, categories as in the dictionary)"
            % (a.write_clean, " translated" if (a.id_map or a.id_regex or a.id_case != "keep") else ""))
    return 1 if n_err else 0


def main(argv=None):
    ap = argparse.ArgumentParser(prog="imp-explorer metadata", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    t = sub.add_parser("template", help="write a metadata dictionary and an empty table")
    t.add_argument("--out", required=True, help="output folder")
    t.add_argument("--data", default="", help="release: one template row per sample, its groups as categories")
    c = sub.add_parser("check", help="check a metadata table")
    c.add_argument("--metadata", required=True, help="TSV/CSV, one row per sample")
    c.add_argument("--dictionary", default="", help="metadata.dictionary.tsv (types, units, allowed values)")
    c.add_argument("--data", default="", help="release: check that the ids match its samples")
    c.add_argument("--id-col", default="", help="id column  [sample_id or the first column]")
    c.add_argument("--id-map", help="TSV/CSV: id in the table, release sample id")
    c.add_argument("--id-regex", help="regular expression; its first group (or match) is the sample id")
    c.add_argument("--id-case", choices=("keep", "upper", "lower"), default="keep")
    c.add_argument("--group-col", default="group", help="group column for the group sizes  [group]")
    c.add_argument("--block", default="site", help="block column for the group sizes  [site]")
    c.add_argument("--min-group", type=int, default=5, help="warn for smaller groups (per block)  [5]")
    c.add_argument("--max-missing", type=float, default=0.3, help="warn above this share of missing values  [0.3]")
    c.add_argument("--report", default="", help="write the issues as TSV")
    c.add_argument("--write-clean", default="", help="write the normalised table (TSV)")
    a = ap.parse_args(argv)
    return cmd_template(a) if a.cmd == "template" else cmd_check(a)


if __name__ == "__main__":
    sys.exit(main())
