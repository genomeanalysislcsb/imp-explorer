"""imp-explorer export: IMP3 sample directories -> one release directory of tables for the explorer.

Samples come from a sheet (--samples: sample_id, imp_dir, any metadata columns) or from a root directory
(--imp-root: every sub-directory with IMP3 output is a sample, named after the directory). A metadata table
(--metadata, TSV/CSV with sample_id or the first column as id) is joined to the samples table.
"""
import argparse
import concurrent.futures as cf
import os
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
    ap.add_argument("--id-case", choices=("keep", "upper", "lower"), default="keep",
                    help="change the case of directory-derived sample ids (e.g. bhc01 -> BHC01) [keep]")
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
    else:
        sheet = pd.DataFrame(discover(a.imp_root), columns=["sample_id", "imp_dir"])
        if a.id_case != "keep":
            sheet["sample_id"] = getattr(sheet["sample_id"].str, a.id_case)()
    if sheet.empty:
        sys.exit("no samples found")
    if a.metadata:
        md = read_sheet(a.metadata)
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
                       {"failed_samples": errors})
    print("[export] %s: %s" % (a.out, ", ".join("%s %d" % (k, v["rows"]) for k, v in m["tables"].items())),
          file=sys.stderr)
    return 0 if not errors else 2


if __name__ == "__main__":
    sys.exit(main())
