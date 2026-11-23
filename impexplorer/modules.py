"""KEGG module completeness of the community (Mantis, per sample) and of every MAG (from the KOs of its genes).

Mantis (IMP3: -km) writes Analysis/annotation/intermediary/mantis/kegg_modules.tsv.gz with the completeness of
each module in the sample's whole gene catalogue, and sample_kos.tsv.gz with the KOs of every protein. With KEGG
module definitions (imp-explorer modules --download FILE, needs access to rest.kegg.jp), the same score is computed
for the genes of each MAG (gene -> contig from Analysis/annotation/annotation_CDS_RNA_hmms.gff, contig -> bin from
Binning/MetaWrap/bins.contigs).

Score (as Mantis): the path through the module definition with the highest share of its KOs present; completeness =
present / (present + missing) KOs of that path. Optional components (-K...) are ignored, a sub-module (M...) counts
as one unit, present when it is complete.

  imp-explorer modules --download kegg_module_definitions.tsv          # once, where rest.kegg.jp is reachable
  imp-explorer modules --release REL [--imp-root ROOT] [--definitions kegg_module_definitions.tsv] [--threads 8]
"""
import argparse
import concurrent.futures as cf
import multiprocessing as mp
import os
import re
import sys
import time
import urllib.request

import pandas as pd

from .readers import imp3

MANTIS_DIR = "Analysis/annotation/intermediary/mantis"
GFF = "Analysis/annotation/annotation_CDS_RNA_hmms.gff"
KEGG_REST = "https://rest.kegg.jp"
COLUMNS = ["sample_id", "scope", "bin_id", "module", "completeness", "kos", "missing"]
INFO_COLUMNS = ["module", "name", "class", "definition"]


# ----------------------------------------------------------------------------- definitions
_TOKEN = re.compile(r"K\d{5}|M\d{5}|--|[()+,\- ]")


def parse(definition):
    """KEGG module definition -> tree of ('ko', id) | ('mod', id) | ('none',) | ('and', [...]) | ('or', [...]);
    ' ' = next step (and), ',' = alternatives (or), '+' = complex subunits (and), '-' = optional subunit (dropped)"""
    toks = [t for t in _TOKEN.findall(re.sub(r"\s+", " ", definition.strip()))]
    pos = [0]

    def peek():
        return toks[pos[0]] if pos[0] < len(toks) else None

    def take():
        pos[0] += 1
        return toks[pos[0] - 1]

    def item():
        t = peek()
        if t == "(":
            take()
            node = alt()
            if peek() == ")":
                take()
            return node
        if t is None or t in (")", ",", " "):
            return ("none",)
        take()
        if t == "--":
            return ("none",)
        return ("ko", t) if t.startswith("K") else ("mod", t) if t.startswith("M") else ("none",)

    def cplx():
        parts = []
        optional = False
        if peek() == "-":                       # leading optional subunit
            take()
            optional = True
        while True:
            node = item()
            if not optional:
                parts.append(node)
            t = peek()
            if t == "+":
                take()
                optional = False
            elif t == "-":
                take()
                optional = True
            else:
                break
        parts = [p for p in parts if p != ("none",)] or [("none",)]
        return parts[0] if len(parts) == 1 else ("and", parts)

    def seq():
        parts = [cplx()]
        while peek() == " ":
            take()
            if peek() in (")", ",", None):
                break
            parts.append(cplx())
        return parts[0] if len(parts) == 1 else ("and", parts)

    def alt():
        parts = [seq()]
        while peek() == ",":
            take()
            parts.append(seq())
        return parts[0] if len(parts) == 1 else ("or", parts)

    return alt()


def _best(node, kos, lam, done):
    """(value, present, missing) of the path maximising present - lam * (present + missing)"""
    kind = node[0]
    if kind == "none":
        return 0.0, frozenset(), frozenset()
    if kind in ("ko", "mod"):
        have = node[1] in kos if kind == "ko" else node[1] in done
        return ((1.0 - lam, frozenset([node[1]]), frozenset()) if have
                else (-lam, frozenset(), frozenset([node[1]])))
    parts = [_best(c, kos, lam, done) for c in node[1]]
    if kind == "and":
        return (sum(p[0] for p in parts), frozenset().union(*[p[1] for p in parts]),
                frozenset().union(*[p[2] for p in parts]))
    return max(parts, key=lambda p: (p[0], len(p[1])))


def score(tree, kos, done=frozenset()):
    """completeness, present ids, missing ids of the best path (maximum share present; Dinkelbach iteration)"""
    lam, best = 0.0, None
    for _ in range(50):
        _, pres, miss = _best(tree, kos, lam, done)
        n = len(pres) + len(miss)
        r = len(pres) / n if n else 0.0
        if best is not None and r <= lam + 1e-12:
            break
        best, lam = (r, pres, miss), r
    return best


def score_all(trees, kos):
    """{module: (completeness, present, missing)} for one KO set; sub-modules resolved first (complete = present)"""
    out, done = {}, set()
    for _ in range(3):                              # nesting depth of KEGG modules is small
        for m, t in trees.items():
            out[m] = score(t, kos, frozenset(done))
        done = {m for m, r in out.items() if r[0] >= 1.0}
    return out


def read_definitions(path):
    d = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    return d[[c for c in INFO_COLUMNS if c in d.columns]]


def download(path, sleep=0.4):
    """KEGG REST: all modules with name, class and definition -> TSV (KEGG terms of use apply)"""
    def get(url):
        for i in range(4):
            try:
                with urllib.request.urlopen(url, timeout=60) as r:
                    return r.read().decode()
            except Exception:
                if i == 3:
                    raise
                time.sleep(2 + 3 * i)
    ids = [ln.split("\t")[0].replace("md:", "") for ln in get(KEGG_REST + "/list/module").splitlines() if ln]
    rows = []
    for i in range(0, len(ids), 10):
        for entry in get(KEGG_REST + "/get/" + "+".join(ids[i:i + 10])).split("///"):
            rec, key = {}, None
            for ln in entry.splitlines():
                if ln[:12].strip():
                    key = ln[:12].strip()
                if key:
                    rec.setdefault(key, []).append(ln[12:].strip())
            if "ENTRY" in rec:
                rows.append({"module": rec["ENTRY"][0].split()[0], "name": " ".join(rec.get("NAME", [""])),
                             "class": " ".join(rec.get("CLASS", [""])),
                             "definition": " ".join(rec.get("DEFINITION", [""]))})
        print("[modules] %d/%d definitions" % (min(i + 10, len(ids)), len(ids)), file=sys.stderr, end="\r")
        time.sleep(sleep)
    print(file=sys.stderr)
    df = pd.DataFrame(rows, columns=INFO_COLUMNS)
    df.to_csv(path, sep="\t", index=False)
    return df


# ----------------------------------------------------------------------------- one sample
def mantis_modules(base):
    """Mantis' community completeness: DataFrame module, name, completeness, kos, missing"""
    t = imp3.table(imp3.find(base, MANTIS_DIR + "/kegg_modules.tsv"))
    if t.empty or "Module_ID" not in t.columns:
        return pd.DataFrame(columns=["module", "name", "completeness", "kos", "missing"])
    na = lambda s: s.replace({"NA": "", "nan": ""})                      # noqa: E731
    return pd.DataFrame({"module": t["Module_ID"], "name": t.get("Module_name", ""),
                         "completeness": imp3.num(t["Completeness_score_mantis"]),
                         "kos": na(t.get("KOs_mantis", "")).str.replace(",", " "),
                         "missing": na(t.get("Missing_KOs_mantis", "")).str.replace(",", " ")})


def gene_kos(base):
    """{gene id: set of KOs} from Mantis' sample_kos.tsv(.gz) (Sample = mantis@<gene id>)"""
    t = imp3.table(imp3.find(base, MANTIS_DIR + "/sample_kos.tsv"))
    if t.empty:
        return {}
    t = t[t.iloc[:, 1].ne("NA") & t.iloc[:, 1].ne("")]
    ids = t.iloc[:, 0].str.replace(r"^mantis@", "", regex=True)
    return {g: set(re.findall(r"K\d{5}", k)) for g, k in zip(ids, t.iloc[:, 1])}


def gene_contigs(base, genes=None):
    """{gene id: contig} from the CDS of the IMP3 annotation GFF (ID attribute)"""
    p = imp3.find(base, GFF)
    if p is None:
        return {}
    out = {}
    rx = re.compile(r"(?:^|;)ID=([^;]+)")
    for chunk in pd.read_csv(p, sep="\t", header=None, comment="#", usecols=[0, 2, 8], dtype=str,
                             chunksize=500000, quoting=3, on_bad_lines="skip"):
        chunk = chunk[chunk[2] == "CDS"]
        ids = chunk[8].str.extract(rx, expand=False)
        for g, c in zip(ids, chunk[0]):
            if isinstance(g, str) and (genes is None or g in genes):
                out[g] = c
    return out


def sample_rows(base, sample, trees=None):
    """-> (rows in COLUMNS, module info rows, notes)"""
    notes = []
    mm = mantis_modules(base)
    info = mm[["module", "name"]].copy()
    rows = []
    if len(mm):
        x = mm[mm["completeness"] > 0]
        rows.append(pd.DataFrame({"sample_id": sample, "scope": "community", "bin_id": "", "module": x["module"],
                                  "completeness": x["completeness"], "kos": x["kos"], "missing": x["missing"]}))
    else:
        notes.append("modules: no Mantis kegg_modules.tsv")
    if trees:
        gk = gene_kos(base)
        if not gk:
            notes.append("modules: no Mantis sample_kos.tsv, no MAG modules")
        else:
            g2c = gene_contigs(base, set(gk))
            c2b = imp3.contig2bin(base)
            per_bin = {}
            for g, ks in gk.items():
                b = c2b.get(g2c.get(g))
                if b is not None:
                    per_bin.setdefault(b, set()).update(ks)
            if len(mm):                                 # the same score on the whole catalogue, against Mantis
                allk = set().union(*gk.values())
                own = score_all(trees, allk)
                cmp = mm.set_index("module")["completeness"]
                both = [m for m in cmp.index if m in own]
                agree = sum(abs(own[m][0] - cmp[m]) < 0.01 for m in both)
                notes.append("modules: community score = Mantis for %d of %d modules" % (agree, len(both)))
            for b, ks in sorted(per_bin.items()):
                res = score_all(trees, ks)
                rows.append(pd.DataFrame([(sample, "bin", sample + ":" + b, m, round(r[0], 3),
                                           " ".join(sorted(r[1])), " ".join(sorted(r[2])))
                                          for m, r in res.items() if r[0] > 0], columns=COLUMNS))
            notes.append("modules: %d genes with KOs, %d on binned contigs, %d bins"
                         % (len(gk), sum(c2b.get(g2c.get(g)) is not None for g in gk), len(per_bin)))
    df = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=COLUMNS)
    return df[COLUMNS], info, notes


def _one(args):
    sid, base, defs_path = args
    try:
        trees = None
        if defs_path:
            d = read_definitions(defs_path)
            trees = {m: parse(x) for m, x in zip(d["module"], d["definition"]) if x}
        return sid, sample_rows(base, sid, trees), None
    except Exception as e:
        return sid, None, "%s: %s" % (type(e).__name__, e)


def main(argv=None):
    from .release import read_manifest, read_table, write_table
    import json
    ap = argparse.ArgumentParser(prog="imp-explorer modules", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--download", metavar="FILE", help="download the KEGG module definitions to FILE and stop")
    ap.add_argument("--release", help="release directory (written by 'imp-explorer export')")
    ap.add_argument("--imp-root", help="base of the relative imp_dir paths of the release [the release's samples "
                    "sheet location is unknown: give it when the release was exported with --paths relative]")
    ap.add_argument("--definitions", help="KEGG module definitions (--download): also score every MAG")
    ap.add_argument("--format", choices=("auto", "parquet", "tsv"), default="auto")
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args(argv)
    if a.download:
        d = download(a.download)
        print("[modules] %d module definitions -> %s" % (len(d), a.download), file=sys.stderr)
        return 0
    if not a.release:
        ap.error("--release or --download is needed")
    s = read_table(a.release, "samples")
    dirs = [p if os.path.isabs(str(p)) else os.path.join(a.imp_root or "", str(p)) for p in s["imp_dir"]]
    jobs = [(str(sid), d, a.definitions) for sid, d in zip(s["sample_id"], dirs) if str(d) and os.path.isdir(d)]
    if len(jobs) < len(s):
        print("[modules] %d of %d samples without an IMP3 folder (use --imp-root)" % (len(s) - len(jobs), len(s)),
              file=sys.stderr)
    parts, info, errors = [], [], {}
    with cf.ProcessPoolExecutor(max_workers=max(1, a.threads), mp_context=mp.get_context("spawn")) as ex:
        for i, (sid, res, err) in enumerate(ex.map(_one, jobs), 1):
            if err:
                errors[sid] = err
                print("[modules] %s: FAILED %s" % (sid, err), file=sys.stderr)
                continue
            df, inf, notes = res
            parts.append(df)
            info.append(inf)
            print("[modules] %d/%d %s: %d rows; %s" % (i, len(jobs), sid, len(df), "; ".join(notes)),
                  file=sys.stderr)
    df = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=COLUMNS)
    inf = pd.concat(info, ignore_index=True).drop_duplicates("module") if info else pd.DataFrame(columns=["module"])
    if a.definitions:
        d = read_definitions(a.definitions)
        inf = d.merge(inf, on="module", how="outer", suffixes=("", "_mantis"))
        inf["name"] = inf["name"].where(inf["name"].fillna("") != "", inf.get("name_mantis"))
        inf = inf.drop(columns=[c for c in inf.columns if c.endswith("_mantis")])
    for c in INFO_COLUMNS:
        if c not in inf.columns:
            inf[c] = ""
    write_table(df, a.release, "modules", a.format)
    write_table(inf[INFO_COLUMNS], a.release, "module_info", a.format)
    m = read_manifest(a.release)
    from .model import TABLES
    for k, v in (("modules", df), ("module_info", inf)):
        m.setdefault("tables", {})[k] = {"rows": int(len(v)), "file": os.path.basename(
            [p for p in (os.path.join(a.release, k + e) for e in (".parquet", ".tsv.gz")) if os.path.exists(p)][0]),
            "doc": TABLES[k]["doc"]}
    m.setdefault("sources", {})["kegg_module_definitions"] = os.path.abspath(a.definitions) if a.definitions else ""
    if errors:
        m["modules_failed_samples"] = errors
    with open(os.path.join(a.release, "manifest.json"), "w") as fh:
        json.dump(m, fh, indent=2)
    print("[modules] %s: modules %d rows (%d community, %d MAG), module_info %d"
          % (a.release, len(df), (df["scope"] == "community").sum(), (df["scope"] == "bin").sum(), len(inf)),
          file=sys.stderr)
    return 0 if not errors else 2


if __name__ == "__main__":
    sys.exit(main())
