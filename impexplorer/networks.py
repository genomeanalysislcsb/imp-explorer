"""imp-explorer networks: add the metabolic networks of IMP-metabolic (github.com/genomeanalysislcsb/IMP-metabolic)
to a release.

Input: the output folders of 'imp-metabolic map' (--mapping) and 'imp-metabolic networks' (--networks), optionally
the IMP-metabolic reference database (--refdb) for readable reaction equations. Writes into the release:

  met_topology            one row per reaction set (specific / all) and network (community, MAG, pan-reactome of a
                          species / genus / family): size, components, degree, clustering, path length, scope
  met_hubs                hub metabolites (betweenness) of the community and pan-reactome networks
  met_mag_coverage        share of each community's reactions found in its MAGs
  met_pan_summary         per taxon and rank: MAGs, pan, core, accessory and rare reactions
  met_community_reactions sample x reaction group: confidence, genes, metaG / metaT reads
  met_reactions           the reactions used: equation with metabolite names, ECs (with --refdb)

Sample ids of IMP-metabolic (IMP3 folder names) are translated with --id-map / --id-regex / --id-case as for the
other layers; MAG ids (sample:bin) follow their sample.
"""
import argparse
import json
import os
import re
import sys

import pandas as pd

from .release import read_manifest, read_table, write_table

TABLES = {"met_topology": ("networks", "topology"), "met_hubs": ("networks", "hubs"),
          "met_mag_coverage": ("networks", "mag_coverage"), "met_pan_summary": ("networks", "pan_summary"),
          "met_community_reactions": ("mapping", "community_reactions")}


def say(msg):
    sys.stderr.write("[networks] %s\n" % msg)


def read_any(folder, name, columns=None):
    for ext in (".parquet", ".tsv.gz", ".tsv"):
        p = os.path.join(folder, name + ext)
        if os.path.exists(p):
            return pd.read_parquet(p, columns=columns) if ext == ".parquet" else \
                pd.read_csv(p, sep="\t", usecols=columns, low_memory=False)
    return None


def equations(refdb, ids):
    """mnx_id, equation (metabolite names), ec_list for the reactions in ids"""
    r = read_any(refdb, "reactions", ["mnx_id", "equation", "ec_list"])
    m = read_any(refdb, "metabolites", ["mnx_id", "name"])
    if r is None:
        return None
    r = r[r["mnx_id"].isin(ids)].copy()
    names = m.drop_duplicates("mnx_id").set_index("mnx_id")["name"] if m is not None else pd.Series(dtype=str)
    rx = re.compile(r"(\S+)@(\S+)")

    def nice(eq):
        return rx.sub(lambda x: names.get(x.group(1), x.group(1)) + ("[e]" if x.group(2) == "MNXD2" else ""),
                      str(eq)).replace(" = ", " <=> ")
    r["equation"] = r["equation"].map(nice)
    r["ec_list"] = r["ec_list"].fillna("").str.replace(";", " ")
    return r


def main(argv=None):
    ap = argparse.ArgumentParser(prog="imp-explorer networks", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True, help="release directory (imp-explorer export)")
    ap.add_argument("--networks", required=True, help="output of imp-metabolic networks")
    ap.add_argument("--mapping", required=True, help="output of imp-metabolic map")
    ap.add_argument("--refdb", help="IMP-metabolic reference database (reaction equations with names)")
    ap.add_argument("--id-map", help="TSV/CSV: IMP-metabolic sample id, release sample id")
    ap.add_argument("--id-regex", help="regular expression; its first group (or match) is the sample id")
    ap.add_argument("--id-case", choices=("keep", "upper", "lower"), default="keep")
    ap.add_argument("--format", choices=("auto", "parquet", "tsv"), default="auto")
    a = ap.parse_args(argv)
    from .export import derive_ids, read_id_map
    t = {}
    for name, (src, f) in TABLES.items():
        df = read_any(a.networks if src == "networks" else a.mapping, f)
        if df is None:
            sys.exit("networks: %s/%s not found" % (a.networks if src == "networks" else a.mapping, f))
        t[name] = df
    names = sorted(set(t["met_community_reactions"]["sample_id"].astype(str)))
    ids = dict(zip(names, derive_ids(names, read_id_map(a.id_map) if a.id_map else None, a.id_regex, a.id_case)))

    def tr_sample(s):
        return s.astype(str).map(lambda x: ids.get(x, x))

    def tr_bin(s):
        return s.astype(str).map(lambda x: ids.get(x.split(":")[0], x.split(":")[0]) + ":" + x.split(":", 1)[1]
                                 if ":" in x else x)
    for name, df in t.items():
        if "sample_id" in df.columns:
            df["sample_id"] = tr_sample(df["sample_id"].fillna(""))
        if name in ("met_topology", "met_hubs"):
            lv = df["level"]
            df.loc[lv == "community", "id"] = tr_sample(df.loc[lv == "community", "id"])
            df.loc[lv == "bin", "id"] = tr_bin(df.loc[lv == "bin", "id"])
    rel = read_table(a.data, "samples")
    rel_ids = set(rel["sample_id"].astype(str))
    mine = set(ids.values())
    say("%d samples, %d in the release%s" % (len(mine), len(mine & rel_ids),
                                             "" if mine & rel_ids else " - check --id-case / --id-regex / --id-map"))
    if a.refdb:
        eq = equations(a.refdb, set(t["met_community_reactions"]["mnx_id"]))
        if eq is not None:
            t["met_reactions"] = eq
    rows = {}
    for name, df in t.items():
        p = write_table(df, a.data, name, a.format)
        rows[name] = {"rows": int(len(df)), "file": os.path.basename(p)}
    m = read_manifest(a.data)
    m.setdefault("tables", {}).update({k: dict(v, doc="IMP-metabolic (imp-explorer networks)") for k, v in rows.items()})
    nm = read_manifest(a.networks) if os.path.exists(os.path.join(a.networks, "manifest.json")) else {}
    m["networks"] = {"networks": os.path.abspath(a.networks), "mapping": os.path.abspath(a.mapping),
                     "impmetabolic": nm.get("impmetabolic", ""), "sets": nm.get("sets", []),
                     "ranks": nm.get("ranks", []), "refdb": os.path.abspath(a.refdb) if a.refdb else ""}
    with open(os.path.join(a.data, "manifest.json"), "w") as fh:
        json.dump(m, fh, indent=2)
    say("%s: %s" % (a.data, ", ".join("%s %d" % (k, v["rows"]) for k, v in rows.items())))
    return 0


if __name__ == "__main__":
    sys.exit(main())
