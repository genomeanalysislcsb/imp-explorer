"""imp-explorer export: IMP3 sample directories -> one release directory of tables for the explorer.

Samples come from a sheet (--samples: sample_id, imp_dir, any metadata columns) or from a root directory
(--imp-root: every sub-directory with IMP3 output is a sample, named after the directory). A metadata table
(--metadata, TSV/CSV with sample_id or the first column as id) is joined to the samples table.

Sample ids of --imp-root folders: the folder name, optionally translated (in this order) by
  --id-map FILE     two columns: folder name, sample id (folders missing from the map keep their name)
  --id-regex REGEX  the first capture group (or the whole match) of REGEX in the folder name
  --id-case         keep | upper | lower
"""
import argparse
import concurrent.futures as cf
import os
import re
import sys

import pandas as pd

from .readers import imp3
from .release import write_manifest, write_table

TABLE_ORDER = ("qc", "taxonomy", "bins", "bin_abundance", "functions")


def discover(root):
    out = []
    for d in sorted(os.listdir(root)):
        p = os.path.join(root, d)
        if os.path.isdir(p) and (os.path.isdir(os.path.join(p, "Stats")) or os.path.isdir(os.path.join(p, "status"))):
            out.append((d, p))
    return out


def derive_ids(names, id_map=None, id_regex=None, id_case="keep"):
    """folder names -> sample ids (see the module doc); raises ValueError on duplicate or empty ids"""
    out = []
    rx = re.compile(id_regex) if id_regex else None
    for n in names:
        sid = (id_map or {}).get(n, n)
        if rx is not None and n not in (id_map or {}):
            m = rx.search(sid)
            if m is None:
                raise ValueError("--id-regex %r does not match folder %r" % (id_regex, n))
            sid = m.group(1) if m.groups() else m.group(0)
        if id_case == "upper":
            sid = sid.upper()
        elif id_case == "lower":
            sid = sid.lower()
        out.append(sid)
    empty = [n for n, s in zip(names, out) if not s]
    if empty:
        raise ValueError("empty sample id for folder(s): %s" % ", ".join(empty))
    dup = sorted({s for s in out if out.count(s) > 1})
    if dup:
        raise ValueError("duplicate sample ids after translation: %s" % ", ".join(dup))
    return out


def read_id_map(path):
    df = pd.read_csv(path, sep="," if path.endswith(".csv") else "\t", dtype=str, header=None,
                     keep_default_na=False, comment="#")
    if df.shape[1] < 2:
        raise ValueError("--id-map needs two columns: folder name, sample id")
    df = df.iloc[:, :2]
    if df.iloc[0, 0].lower() in ("folder", "dir", "directory", "imp_dir", "name") or df.iloc[0, 1].lower() == "sample_id":
        df = df.iloc[1:]
    return dict(zip(df.iloc[:, 0].str.strip(), df.iloc[:, 1].str.strip()))


def store_paths(paths, mode, root=None):
    """how imp_dir is stored in the release: absolute, relative (to --imp-root / the sheet) or none"""
    if mode == "absolute":
        return [os.path.abspath(p) for p in paths]
    if mode == "none":
        return [""] * len(paths)
    base = os.path.abspath(root) if root else None
    return [os.path.relpath(os.path.abspath(p), base) if base else os.path.basename(os.path.normpath(p))
            for p in paths]


def read_sheet(path):
    sep = "," if path.endswith(".csv") else "\t"
    df = pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False, comment="#")
    if "sample_id" not in df.columns:
        df = df.rename(columns={df.columns[0]: "sample_id"})
    return df


def _one(args):
    sid, base, dbs = args
    try:
        return sid, base, imp3.read_sample(base, sid, dbs), None
    except Exception as e:                       # one broken sample must not stop the export
        return sid, base, None, "%s: %s" % (type(e).__name__, e)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="imp-explorer export", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--samples", help="TSV/CSV with sample_id, imp_dir (+ metadata columns)")
    src.add_argument("--imp-root", help="directory whose sub-directories are IMP3 sample outputs")
    ap.add_argument("--metadata", help="TSV/CSV of sample metadata (sample_id or first column = id)")
    ap.add_argument("--id-map", help="TSV/CSV: folder name, sample id (with --imp-root)")
    ap.add_argument("--id-regex", help="regular expression; its first group (or match) is the sample id (with --imp-root)")
    ap.add_argument("--id-case", choices=("keep", "upper", "lower"), default="keep",
                    help="change the case of folder-derived sample ids [keep]")
    ap.add_argument("--paths", choices=("relative", "absolute", "none"), default="relative",
                    help="how the IMP3 folders are stored in the release (imp_dir): relative to --imp-root or the "
                         "sample sheet, absolute, or not at all [relative]")
    ap.add_argument("--out", required=True, help="release directory to write (e.g. release/v0.1)")
    ap.add_argument("--study", default="study", help="study name shown in the explorer")
    ap.add_argument("--function-dbs", default=",".join(imp3.DEFAULT_FUNCTION_DBS),
                    help="Mantis databases to export [%(default)s]")
    ap.add_argument("--format", choices=("auto", "parquet", "tsv"), default="auto")
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args(argv)

    if a.samples:
        sheet = read_sheet(a.samples)
        if "imp_dir" not in sheet.columns:
            sys.exit("--samples needs a column imp_dir")
        sheet_dir = os.path.dirname(os.path.abspath(a.samples))
        sheet["imp_dir"] = [p if os.path.isabs(p) else os.path.join(sheet_dir, p) for p in sheet["imp_dir"]]
        root = sheet_dir
    else:
        found = discover(a.imp_root)
        try:
            ids = derive_ids([n for n, _ in found], read_id_map(a.id_map) if a.id_map else None, a.id_regex,
                             a.id_case)
        except ValueError as e:
            sys.exit("[export] %s" % e)
        sheet = pd.DataFrame({"sample_id": ids, "imp_dir": [p for _, p in found]})
        root = a.imp_root
    if sheet.empty:
        sys.exit("no samples found")
    unmatched_md = []
    if a.metadata:
        md = read_sheet(a.metadata)
        unmatched_md = sorted(set(md["sample_id"]) - set(sheet["sample_id"]))
        no_md = sorted(set(sheet["sample_id"]) - set(md["sample_id"]))
        print("[export] metadata: %d of %d samples matched%s" % (
            len(sheet) - len(no_md), len(sheet), ("; no metadata for: %s" % ", ".join(no_md[:10]) +
                                                  (" ..." if len(no_md) > 10 else "")) if no_md else ""),
              file=sys.stderr)
        if unmatched_md:
            print("[export] %d metadata ids without IMP3 output (check --id-map / --id-regex / --id-case): %s"
                  % (len(unmatched_md), ", ".join(unmatched_md[:10]) + (" ..." if len(unmatched_md) > 10 else "")),
                  file=sys.stderr)
        sheet = sheet.merge(md, on="sample_id", how="left", suffixes=("", "_metadata"))
    os.makedirs(a.out, exist_ok=True)
    dbs = tuple(x for x in a.function_dbs.split(",") if x)

    parts = {k: [] for k in TABLE_ORDER}
    omics, errors = {}, {}
    jobs = list(zip(sheet["sample_id"], sheet["imp_dir"], [dbs] * len(sheet)))
    with cf.ProcessPoolExecutor(max_workers=max(1, a.threads)) as ex:
        for i, (sid, base, res, err) in enumerate(ex.map(_one, jobs), 1):
            if err:
                errors[sid] = err
                print("[export] %s: FAILED %s" % (sid, err), file=sys.stderr)
                continue
            for k in TABLE_ORDER:
                if len(res[k]):
                    parts[k].append(res[k])
            omics[sid] = res["_omics"]
            print("[export] %d/%d %s: %s" % (i, len(jobs), sid, ", ".join(
                "%s %d" % (k, len(res[k])) for k in TABLE_ORDER)), file=sys.stderr)

    sheet["imp_dir"] = store_paths(sheet["imp_dir"], a.paths, root)
    sheet["has_mg"] = sheet["sample_id"].map(lambda s: "mg" in omics.get(s, []))
    sheet["has_mt"] = sheet["sample_id"].map(lambda s: "mt" in omics.get(s, []))
    first = ["sample_id", "imp_dir", "has_mg", "has_mt"]
    sheet = sheet[first + [c for c in sheet.columns if c not in first]]
    written = {"samples": (sheet, write_table(sheet, a.out, "samples", a.format))}
    from .model import TABLES
    for k in TABLE_ORDER:
        df = pd.concat(parts[k], ignore_index=True) if parts[k] else pd.DataFrame(columns=TABLES[k]["columns"])
        written[k] = (df, write_table(df, a.out, k, a.format))
    m = write_manifest(a.out, written, a.study, {"pipeline": "IMP3", "function_dbs": list(dbs)},
                       {"failed_samples": errors, "metadata_without_output": unmatched_md})
    print("[export] %s: %s" % (a.out, ", ".join("%s %d" % (k, v["rows"]) for k, v in m["tables"].items())),
          file=sys.stderr)
    return 0 if not errors else 2


if __name__ == "__main__":
    sys.exit(main())
