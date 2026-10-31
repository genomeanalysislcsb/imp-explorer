"""Readers for one IMP3 sample output directory (IMP3 3.0, integrative mg + mt runs; single-omic runs work too).

Every function takes the sample directory and the sample id and returns a pandas DataFrame in the release model
(impexplorer.model.TABLES); missing files give an empty frame, never an error, so partly finished samples still load.
"""
import glob
import gzip
import io
import os
import re

import pandas as pd

from ..model import GTDB_RANKS, KRAKEN_RANKS, TABLES, mimag_tier

OMICS = ("mg", "mt")
DEFAULT_FUNCTION_DBS = ("kegg_ko", "kegg_module", "kegg_pathway", "cazy", "enzyme_ec")
FLAGSTAT = ("total", "primary", "secondary", "supplementary", "duplicates", "primary_duplicates", "mapped",
            "primary_mapped")


# ----------------------------------------------------------------------------- helpers
def find(base, *patterns):
    """first existing file among the patterns (each also tried with / without .gz)"""
    for p in patterns:
        for cand in (p, p + ".gz") if not p.endswith(".gz") else (p, p[:-3]):
            hits = sorted(glob.glob(os.path.join(base, cand)))
            hits = [h for h in hits if os.path.isfile(h) and os.path.getsize(h) > 0]
            if hits:
                return hits[0]
    return None


def text(path):
    if path is None:
        return ""
    with open(path, "rb") as fh:
        raw = fh.read()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    return raw.decode("utf-8", errors="replace")


def table(path, header="infer", comment=None, names=None):
    """TSV (plain or gzipped; cells padded with spaces as in Stats/all_bin_stats.tsv are stripped)"""
    t = text(path)
    if not t.strip():
        return pd.DataFrame()
    df = pd.read_csv(io.StringIO(t), sep="\t", header=0 if header == "infer" else header, comment=comment,
                     names=names, dtype=str, keep_default_na=False)
    df.columns = [str(c).strip() for c in df.columns]
    for c in df.columns:
        df[c] = df[c].str.strip()
    return df


def num(s):
    return pd.to_numeric(s, errors="coerce")


def empty(name):
    return pd.DataFrame(columns=TABLES[name]["columns"])


def assembly_prefix(base):
    """'mgmt' for integrative runs, 'mg' / 'mt' for single-omic runs"""
    for p in sorted(glob.glob(os.path.join(base, "Stats", "*", "*.assembly.contig_read_count.txt*"))):
        return os.path.basename(p).split(".assembly.")[0]
    for p in sorted(glob.glob(os.path.join(base, "Assembly", "*.assembly.merged.fa*"))):
        return os.path.basename(p).split(".assembly.")[0]
    return "mgmt"


def omics_present(base):
    return [o for o in OMICS if find(base, "Stats/%s/%s.read_stats.tsv" % (o, o), "Stats/%s/%s.stats" % (o, o),
                                       "Preprocessing/%s.r1.preprocessed.fq.gz" % o)]


# ----------------------------------------------------------------------------- QC
def qc(base, sample):
    rows = []

    def add(omic, metric, value):
        v = num(pd.Series([value])).iloc[0]
        if pd.notna(v):
            rows.append((sample, omic, metric, float(v)))

    # IMP3 summary statistics (Stats/all_stats.tsv: type, step, number; reads per step, mapping, annotation)
    st = table(find(base, "Stats/all_stats.tsv"))
    if {"type", "step", "number"} <= set(st.columns):
        for _, r in st.iterrows():
            add(r["type"], "imp:" + r["step"], r["number"])
    # read quality of raw and trimmed R1 (seqkit stats)
    for o in OMICS:
        rs = table(find(base, "Stats/%s/%s.read_stats.tsv" % (o, o)))
        if "file" in rs.columns:
            for _, r in rs.iterrows():
                tag = {"%s.r1.fq" % o: "raw", "%s.r1.trimmed.fq" % o: "trimmed"}.get(r["file"])
                if tag:
                    for col, m in (("num_seqs", "pairs"), ("avg_len", "avg_len"), ("Q30(%)", "q30_pct"),
                                   ("GC(%)", "gc_pct")):
                        add(o, "%s:%s" % (tag, m), r.get(col))
        # Nonpareil diversity (mean of R1 / R2)
        vals = []
        for rr in ("r1", "r2"):
            t = table(find(base, "Stats/%s/%s.%s_nonpareil.diversity.txt" % (o, o, rr)), header=None,
                      names=["k", "v"])
            if len(t):
                vals += list(num(t.loc[t["k"].str.contains("diversity"), "v"]))
        if vals:
            add(o, "nonpareil_diversity", sum(vals) / len(vals))
    # mapping to the assembly (samtools flagstat, first column)
    ass = assembly_prefix(base)
    for o in OMICS:
        t = text(find(base, "Stats/%s/%s.assembly.contig_flagstat.txt" % (o, ass))).split()
        vals = dict(zip(FLAGSTAT, t)) if len(t) >= len(FLAGSTAT) else {}
        for k in ("total", "primary", "mapped", "primary_mapped"):
            add(o, "flagstat:" + k, vals.get(k))
        if vals:
            p, pm = num(pd.Series([vals["primary"], vals["primary_mapped"]]))
            if p and pm == pm:
                add(o, "mapped_fraction", pm / p)
    # assembly
    ln = table(find(base, "Stats/%s/%s.assembly.length.txt" % (ass, ass)), header=None, names=["contig", "len"])
    if len(ln):
        L = num(ln["len"]).dropna().sort_values(ascending=False)
        cs = L.cumsum()
        add(ass, "contigs", len(L))
        add(ass, "assembly_bp", L.sum())
        add(ass, "contigs_ge_1kb", int((L >= 1000).sum()))
        add(ass, "n50", L[cs >= L.sum() / 2].iloc[0] if len(L) else None)
        add(ass, "longest_contig", L.iloc[0] if len(L) else None)
    return pd.DataFrame(rows, columns=TABLES["qc"]["columns"]) if rows else empty("qc")


# ----------------------------------------------------------------------------- taxonomy
def _kreport(path, sample, omic, rank_code):
    """Kraken-style report (also Bracken's kraken_bracken_<rank>.report): pct, clade reads, direct reads, rank,
    taxid, name"""
    t = text(path)
    rows = []
    for line in t.splitlines():
        f = line.split("\t")
        if len(f) < 6 or f[3].strip() != rank_code:
            continue
        rows.append((sample, omic, "kraken2_bracken", KRAKEN_RANKS[rank_code], f[5].strip(), f[4].strip(),
                     float(f[1])))
    return rows


def _bracken(path, sample, omic, rank):
    t = table(path)
    if "name" not in t.columns:
        return []
    return [(sample, omic, "kraken2_bracken", rank, r["name"], r.get("taxonomy_id", ""), float(r["new_est_reads"]))
            for _, r in t.iterrows()]


MOTU_RANKS = {"k": "domain", "p": "phylum", "c": "class", "o": "order", "f": "family", "g": "genus", "s": "species"}


def _motus(path, sample, omic):
    """mOTUs profile: '#'-comment lines (last: '#consensus_taxonomy\t<sample>'); plain (one row per mOTU) or
    hierarchical (-p -A / -k all: 'k__Bacteria|p__...|...' rows for every rank). Ranked rows keep their rank
    (domain ... species, 'motu' for single mOTUs); 'unassigned' rows are dropped."""
    rows = []
    for line in text(path).splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 2:
            continue
        v = num(pd.Series([f[-1]])).iloc[0]
        if pd.isna(v) or v <= 0:
            continue
        leaf = f[0].strip().split("|")[-1].strip()
        if leaf.lower().startswith("unassigned") or leaf == "-1":
            continue
        m = re.search(r"\[((?:ref|meta|ext)_mOTU_v\d+_\d+)\]", leaf)
        pm = re.match(r"([kpcofgst])__(.*)$", leaf)
        if m or (pm and pm.group(1) == "t"):
            rank, name = "motu", re.sub(r"^[a-z]__", "", leaf)
        elif pm:
            rank, name = MOTU_RANKS[pm.group(1)], pm.group(2)
        else:
            rank, name = "motu", leaf
        rows.append((sample, omic, "motus", rank, re.sub(r"\s*\[[^\]]+\]\s*$", "", name).strip(),
                     m.group(1) if m else "", float(v)))
    return rows


def taxonomy(base, sample, ranks=("D", "P", "C", "O", "F", "G", "S")):
    rows = []
    level_name = {"D": "domains", "P": "phylums", "C": "classes", "O": "orders", "F": "families", "G": "genuses",
                  "S": "species"}
    for o in OMICS:
        for rc in ranks:
            # IMP3 >= 3.0: Analysis/taxonomy/bracken/<omic>/<rank>.bracken; older: kraken/<omic>.kraken_bracken_*.report
            b = find(base, "Analysis/taxonomy/bracken/%s/%s.bracken" % (o, rc))
            if b:
                rows += _bracken(b, sample, o, KRAKEN_RANKS[rc])
                continue
            k = find(base, "Analysis/taxonomy/kraken/%s.kraken_bracken_%s.report" % (o, level_name[rc]),
                     "Analysis/taxonomy/bracken/%s/kraken_bracken_%s.report" % (o, level_name[rc]))
            if k:
                rows += _kreport(k, sample, o, rc)
        m = find(base, "Analysis/taxonomy/mOTUs/%s.mOTU.counts.tsv" % o)
        if m:
            rows += _motus(m, sample, o)
    if not rows:
        return empty("taxonomy")
    df = pd.DataFrame(rows, columns=["sample_id", "omic", "source", "rank", "taxon", "taxid", "reads"])
    df = df[df["reads"] > 0]
    tot = df.groupby(["omic", "source", "rank"])["reads"].transform("sum")
    df["rel_abundance"] = df["reads"] / tot
    return df[TABLES["taxonomy"]["columns"]].reset_index(drop=True)


# ----------------------------------------------------------------------------- bins
_BIN_COLS = {"bin": "bin", "completeness": "completeness", "contamination": "contamination",
             "contig_n50": "contig_n50", "genome_size": "genome_size", "gc_content": "gc",
             "classification": "classification", "closest_genome_reference": "closest_genome_reference",
             "closest_genome_ani": "closest_genome_ani", "grid": "grid", "coverage": "coverage"}


def bins(base, sample):
    ass = assembly_prefix(base)
    t = table(find(base, "Stats/all_bin_stats.tsv", "Stats/%s/%s.bins.tsv" % (ass, ass)))
    if t.empty or "Bin" not in t.columns:
        return empty("bins")
    t.columns = [re.sub(r"[^a-z0-9]+", "_", c.lower()).strip("_") for c in t.columns]
    df = pd.DataFrame({new: t[old] if old in t.columns else "" for old, new in _BIN_COLS.items()})
    for c in ("completeness", "contamination", "contig_n50", "genome_size", "gc", "closest_genome_ani", "grid",
              "coverage"):
        df[c] = num(df[c])
    df.insert(0, "sample_id", sample)
    df.insert(0, "bin_id", sample + ":" + df["bin"])
    df["quality"] = [mimag_tier(c, x) for c, x in zip(df["completeness"], df["contamination"])]
    ranks = df["classification"].fillna("").str.split(";")
    for i, r in enumerate(GTDB_RANKS):
        df[r] = [(x[i][3:] if len(x) > i and len(x[i]) > 3 else "") for x in ranks]
    return df[TABLES["bins"]["columns"]]


def contig2bin(base):
    t = table(find(base, "Binning/MetaWrap/bins.contigs"), header=None, names=["contig", "bin"])
    return dict(zip(t["contig"], t["bin"])) if len(t) else {}


def bin_abundance(base, sample):
    c2b = contig2bin(base)
    if not c2b:
        return empty("bin_abundance")
    ass = assembly_prefix(base)
    out = []
    for o in OMICS:
        rc = table(find(base, "Stats/%s/%s.assembly.contig_read_count.txt" % (o, ass)), header=None,
                   names=["contig", "start", "end", "reads"])
        if rc.empty:
            continue
        rc["reads"] = num(rc["reads"]).fillna(0)
        rc["length_bp"] = num(rc["end"]) - num(rc["start"])
        dp = table(find(base, "Stats/%s/%s.assembly.contig_depth.txt" % (o, ass)), header=None,
                   names=["contig", "depth"])
        depth = dict(zip(dp["contig"], num(dp["depth"]))) if len(dp) else {}
        rc["bin"] = rc["contig"].map(c2b)
        rc["depth_bp"] = rc["contig"].map(depth).fillna(0) * rc["length_bp"]
        mapped = rc["reads"].sum()
        g = rc.dropna(subset=["bin"]).groupby("bin").agg(reads=("reads", "sum"), length_bp=("length_bp", "sum"),
                                                          depth_bp=("depth_bp", "sum")).reset_index()
        g["mean_depth"] = g["depth_bp"] / g["length_bp"].where(g["length_bp"] > 0)
        g["share_mapped"] = g["reads"] / mapped if mapped else float("nan")
        g["rpkm"] = g["reads"] / (g["length_bp"] / 1e3) / (mapped / 1e6) if mapped else float("nan")
        g["omic"] = o
        g["sample_id"] = sample
        g["bin_id"] = sample + ":" + g["bin"]
        out.append(g)
    if not out:
        return empty("bin_abundance")
    return pd.concat(out, ignore_index=True)[TABLES["bin_abundance"]["columns"]]


# ----------------------------------------------------------------------------- functions
def _featurecounts(path):
    """featureCounts table -> {feature: reads} (sum of the BAM columns: paired + single reads)"""
    t = table(path, comment="#")
    if t.empty or "Geneid" not in t.columns:
        return {}
    cnt = t[t.columns[6:]].apply(num).fillna(0).sum(axis=1)
    s = pd.Series(cnt.values, index=t["Geneid"])
    s = s.groupby(level=0).sum()
    return s[s > 0].to_dict()


def functions(base, sample, dbs=DEFAULT_FUNCTION_DBS):
    rows = []
    for o in OMICS:
        for db in dbs:
            for f, v in _featurecounts(find(base, "Analysis/annotation/%s.mantis.%s.counts.tsv" % (o, db))).items():
                rows.append((sample, o, db, f, v))
    return pd.DataFrame(rows, columns=TABLES["functions"]["columns"]) if rows else empty("functions")


# ----------------------------------------------------------------------------- one sample
def read_sample(base, sample, function_dbs=DEFAULT_FUNCTION_DBS):
    return {"qc": qc(base, sample), "taxonomy": taxonomy(base, sample), "bins": bins(base, sample),
            "bin_abundance": bin_abundance(base, sample), "functions": functions(base, sample, function_dbs),
            "_omics": omics_present(base)}
