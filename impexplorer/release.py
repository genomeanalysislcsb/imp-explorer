"""Release directory: one file per table (Parquet if pyarrow is installed, else TSV.gz) + manifest.json."""
import datetime
import json
import os

import pandas as pd

from . import __version__
from .model import TABLES


def parquet_ok():
    try:
        import pyarrow  # noqa: F401
        return True
    except ImportError:
        return False


def write_table(df, outdir, name, fmt="auto"):
    fmt = ("parquet" if parquet_ok() else "tsv") if fmt == "auto" else fmt
    path = os.path.join(outdir, name + (".parquet" if fmt == "parquet" else ".tsv.gz"))
    if fmt == "parquet":
        df.to_parquet(path, index=False)
    else:
        df.to_csv(path, sep="\t", index=False, compression="gzip")
    return path


def read_table(release, name):
    for ext in (".parquet", ".tsv.gz", ".tsv"):
        p = os.path.join(release, name + ext)
        if os.path.exists(p):
            return pd.read_parquet(p) if ext == ".parquet" else pd.read_csv(p, sep="\t")
    return pd.DataFrame(columns=TABLES.get(name, {}).get("columns", []))


def table_path(release, name):
    for ext in (".parquet", ".tsv.gz", ".tsv"):
        p = os.path.join(release, name + ext)
        if os.path.exists(p):
            return p
    return None


def write_manifest(outdir, tables, study, sources, extra=None):
    m = {"release": os.path.basename(os.path.abspath(outdir)), "study": study,
         "created": datetime.datetime.now().isoformat(timespec="seconds"), "impexplorer": __version__,
         "tables": {k: {"rows": int(len(v)), "file": os.path.basename(p), "doc": TABLES.get(k, {}).get("doc", "")}
                    for k, (v, p) in tables.items()},
         "sources": sources}
    m.update(extra or {})
    with open(os.path.join(outdir, "manifest.json"), "w") as fh:
        json.dump(m, fh, indent=2)
    return m


def read_manifest(release):
    p = os.path.join(release, "manifest.json")
    return json.load(open(p)) if os.path.exists(p) else {}
