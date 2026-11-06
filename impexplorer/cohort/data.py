"""Cohort input: samples + metadata, QC flags, per-sample metrics and feature layers (count matrices) of a release."""
import math
import os
import re
import sys
from collections import OrderedDict

import numpy as np
import pandas as pd

from ..release import read_manifest, read_table
from .stats import clr, column_kind, is_missing

GLOBAL_MAPS = re.compile(r"^(ko|map)01[12]\d\d$")       # KEGG global / overview maps: not tested
LAYER_LABEL = {"kraken": "Kraken2/Bracken", "motus": "mOTUs", "mag": "MAGs (GTDB)", "kegg_module": "KEGG module",
               "kegg_pathway": "KEGG pathway", "kegg_ko": "KEGG KO", "cazy": "CAZy", "enzyme_ec": "EC", "pfam": "Pfam", "cog": "COG"}
OMIC = {"mg": "metagenome", "mt": "metatranscriptome"}

# per-sample metrics: name -> (label, unit/scale note); computed in Cohort.metrics_table
METRIC_LABEL = OrderedDict([
    ("mg_pairs_M", "mg read pairs after trimming (M)"), ("mt_pairs_M", "mt read pairs after trimming (M)"),
    ("mg_removed_fraction", "mg reads removed after trimming (PhiX, host)"),
    ("mt_rrna_fraction", "mt rRNA reads (share of trimmed)"),
    ("mt_removed_fraction", "mt reads removed after trimming (rRNA, PhiX, host)"),
    ("mg_mapped_fraction", "mg reads mapped to the assembly"), ("mt_mapped_fraction", "mt reads mapped to the assembly"),
    ("mg_nonpareil", "Nonpareil diversity (mg)"), ("assembly_Mb", "assembly length (Mb)"),
    ("assembly_n50_kb", "assembly N50 (kb)"), ("mags_hq", "high-quality MAGs"), ("mags_mq", "medium-quality MAGs"),
    ("mg_binned_fraction", "mg reads on MAGs (share of mapped)"),
    ("mt_binned_fraction", "mt reads on MAGs (share of mapped)"),
])
QC_OUTLIER = {"mg_mapped_fraction": 0.05, "mt_mapped_fraction": 0.05, "mg_removed_fraction": 0.05,
              "mt_rrna_fraction": 0.10, "mg_nonpareil": 0.5, "assembly_n50_kb": 0.5}


def norm_id(x):
    return str(x).strip().upper().split(".")[0]


class Layer:
    """one feature set: counts (samples x features), for ratio layers values (CLR mt - CLR mg), for metabolomics
    (metab) the measured values (concentrations / normalised intensities; >= 0, 0 = not detected)"""

    def __init__(self, key, label, omic, kind, M, db=None):
        self.key, self.label, self.omic, self.kind, self.M, self.db = key, label, omic, kind, M, db

    def rel(self):
        if self.kind not in ("counts", "auc"):
            return self.M
        return self.M.div(self.M.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)

    def tested(self, prevalence, min_reads=10):
        """features with >= min_reads in >= prevalence of the samples (counts), detected (> 0) in >= prevalence of
        the samples and not constant (metab), all features (ratio)"""
        if self.kind == "auc":                                   # proteins: detected in >= prevalence
            return list(self.M.columns[(self.M > 0).mean(axis=0) >= prevalence])
        if self.kind == "metab":
            ok = ((self.M > 0).mean(axis=0) >= prevalence) & (self.M.std(axis=0) > 0)
            return list(self.M.columns[ok])
        if self.kind != "counts":
            return list(self.M.columns)
        return list(self.M.columns[(self.M >= min_reads).mean(axis=0) >= prevalence])


class Cohort:
    def __init__(self, a):
        self.a = a
        self.release = a.data
        self.manifest = read_manifest(a.data)
        self.study = self.manifest.get("study") or os.path.basename(os.path.abspath(a.data))
        self.log = []
        s = read_table(a.data, "samples")
        s["sample_id"] = s["sample_id"].astype(str)
        self.unmatched_meta = []
        if a.metadata:
            s = self._merge_metadata(s, a.metadata, a.id_col)
        dict_kinds = {}
        if getattr(a, "dictionary", ""):
            s, dict_kinds = self._apply_dictionary(s, a.dictionary)
        self.samples = s.set_index("sample_id")
        self._qct = read_table(a.data, "qc")
        self._metrics = self.metrics_table()
        self.metric_names = list(self._metrics.columns)
        # groups and blocks
        gc = a.group_col if a.group_col in self.samples.columns else ""
        if a.group_col and not gc:
            self.note("group column '%s' not in the metadata: no group comparisons" % a.group_col)
        self.group_col = gc
        self.block = a.block if a.block in self.samples.columns else ""
        if a.block and not self.block:
            self.note("block column '%s' not in the metadata: ignored" % a.block)
        vals = [v for v in self.samples[gc].map(self._s).unique() if v != "NA"] if gc else []
        order = [g.strip() for g in (a.group_order or "").split(",") if g.strip()]
        self.groups = [g for g in order if g in vals] + sorted(v for v in vals if v not in order)
        self.ref = a.reference if a.reference in self.groups else (self.groups[0] if self.groups else "")
        # metadata columns and kinds
        hidden = {"imp_dir", "has_mg", "has_mt"}
        self.meta_cols = [c for c in self.samples.columns if c not in hidden]
        self.kinds = {c: column_kind([self._s(v) for v in self.samples[c]]) for c in self.meta_cols}
        for c, k in dict_kinds.items():          # the dictionary decides numeric / categorical; text / dates are no
            if c in self.kinds and self.kinds[c] != "empty":   # covariates
                self.kinds[c] = k
        cov = [c.strip() for c in (a.covariates or "").split(",") if c.strip()]
        self.covariates = [c for c in (cov or self.meta_cols) if c in self.meta_cols and c != gc
                           and self.kinds.get(c) not in ("empty", "text")]
        self.qc_flags = self.flag_samples()
        drop = {"FAIL": ("FAIL",), "WARN": ("FAIL", "WARN"), "none": ()}[a.qc_drop]
        self.kept = [x for x in self.samples.index if self.qc_flags.loc[x, "level"] not in drop]
        self.excluded = [x for x in self.samples.index if x not in self.kept]
        self.layers = OrderedDict()

    # ------------------------------------------------------------------ metadata
    @staticmethod
    def _s(v):
        return "NA" if v is None or (isinstance(v, float) and math.isnan(v)) or is_missing(v) else str(v).strip()

    def _merge_metadata(self, s, path, id_col):
        sep = "\t" if open(path, encoding="utf-8-sig").read(4096).count("\t") else ","
        m = pd.read_csv(path, sep=sep, dtype=str, comment=None, encoding="utf-8-sig")
        m = m[~m.iloc[:, 0].astype(str).str.startswith("#")]
        idc = id_col if id_col in m.columns else ("sample_id" if "sample_id" in m.columns else m.columns[0])
        m["_key"] = m[idc].map(norm_id)
        s["_key"] = s["sample_id"].map(norm_id)
        known = set(s["_key"])
        self.unmatched_meta = sorted(m.loc[~m["_key"].isin(known), idc].astype(str))
        m = m.drop_duplicates("_key").drop(columns=[idc])
        over = [c for c in m.columns if c in s.columns and c != "_key"]
        s = s.drop(columns=over).merge(m, on="_key", how="left").drop(columns=["_key"])
        n = int(s[[c for c in m.columns if c != "_key"]].notna().any(axis=1).sum()) if len(m.columns) > 1 else 0
        self.note("metadata %s: %d of %d samples matched, %d metadata ids without IMP3 output"
                  % (os.path.basename(path), n, len(s), len(self.unmatched_meta)))
        return s

    def _apply_dictionary(self, s, path):
        """metadata dictionary: missing codes -> NA, numbers / categories normalised, kinds per column"""
        from ..metadata import normalize, read_dictionary
        dct = read_dictionary(path)
        cols = [c for c in s.columns if c in dct and c != "sample_id"]
        sub = s[cols].fillna("").astype(str)
        clean, issues, kinds = normalize(sub, dct)
        for c in cols:
            s[c] = clean[c].replace("", np.nan)
        for lvl, col, msg, ex in issues:
            if lvl == "ERROR":
                self.note("dictionary: %s: %s (set to NA)%s" % (col, msg, (" e.g. " + ", ".join(map(str, ex))) if ex else ""))
        self.note("dictionary %s: %d columns typed (%d numeric, %d categorical)" % (
            os.path.basename(path), len(kinds), sum(1 for k in kinds.values() if k == "numeric"),
            sum(1 for k in kinds.values() if k == "categorical")))
        return s, kinds

    def note(self, msg):
        self.log.append(msg)
        sys.stderr.write("[cohort] %s\n" % msg)

    def meta(self, sample, col):
        if col == "__group__":
            col = self.group_col
        if col in self.samples.columns:
            return self._s(self.samples.at[sample, col])
        if col in self._metrics.columns:
            v = self._metrics.at[sample, col]
            return "NA" if pd.isna(v) else repr(float(v))
        return "NA"

    def group(self, sample):
        return self.meta(sample, self.group_col) if self.group_col else "all"

    def block_of(self, sample):
        return self.meta(sample, self.block) if self.block else "all"

    # ------------------------------------------------------------------ metrics
    def _qc(self, omic, metric):
        q = self._qct
        d = q[(q["omic"] == omic) & (q["metric"] == metric)]
        return d.groupby("sample_id")["value"].first()

    def _qc_re(self, omic, pattern):
        q = self._qct
        d = q[(q["omic"] == omic) & q["metric"].str.contains(pattern, regex=True)]
        return d.groupby(["sample_id", "metric"])["value"].first().unstack()

    def metrics_table(self):
        idx = self.samples.index
        m = pd.DataFrame(index=idx)
        m["mg_pairs_M"] = self._qc("mg", "trimmed:pairs") / 1e6
        m["mt_pairs_M"] = self._qc("mt", "trimmed:pairs") / 1e6
        for o in ("mg", "mt"):                     # read numbers per IMP3 filtering step (R1)
            st = self._qc_re(o, r"\b%s\.r1\.trimmed" % o)
            if st.empty:
                continue
            trimmed = [c for c in st.columns if re.search(r"\.r1\.trimmed\.fq$", c)]
            steps = [c for c in st.columns if "__" in c]
            if trimmed and steps:
                last = max(steps, key=lambda c: c.count("__"))
                m["%s_removed_fraction" % o] = 1 - st[last] / st[trimmed[0]]
            rr = [c for c in st.columns if re.search(r"\.r1\.trimmed\.rrna\.fq$", c)]
            if o == "mt" and rr and trimmed:
                m["mt_rrna_fraction"] = st[rr[0]] / st[trimmed[0]]
        m["mg_mapped_fraction"] = self._qc("mg", "mapped_fraction")
        m["mt_mapped_fraction"] = self._qc("mt", "mapped_fraction")
        m["mg_nonpareil"] = self._qc("mg", "nonpareil_diversity")
        m["assembly_Mb"] = self._qc("mgmt", "assembly_bp") / 1e6
        m["assembly_n50_kb"] = self._qc("mgmt", "n50") / 1e3
        bins = read_table(self.release, "bins")
        if len(bins):
            c = bins.groupby(["sample_id", "quality"]).size().unstack(fill_value=0)
            m["mags_hq"] = c.get("high", pd.Series(dtype=float)).reindex(idx).fillna(0)
            m["mags_mq"] = c.get("medium", pd.Series(dtype=float)).reindex(idx).fillna(0)
        ba = read_table(self.release, "bin_abundance")
        if len(ba):
            sh = ba.groupby(["sample_id", "omic"])["share_mapped"].sum().unstack()
            for o in ("mg", "mt"):
                if o in sh.columns:
                    m["%s_binned_fraction" % o] = sh[o]
        m = m.reindex(idx)
        return m[[c for c in METRIC_LABEL if c in m.columns and m[c].notna().any()]]

    def add_alpha(self):
        """Shannon and richness of the main taxonomic layers -> metrics"""
        for key in ("kraken.genus.mg", "kraken.genus.mt", "kraken.species.mg", "motus.species.mg",
                    "motus.species.mt"):
            L = self.layers.get(key)
            if L is None:
                continue
            R = L.rel()
            with np.errstate(divide="ignore", invalid="ignore"):
                H = -(R * np.log(R.where(R > 0, 1))).sum(axis=1)
            tag = key.replace(".", "_")
            self._metrics.loc[R.index, "shannon_%s" % tag] = H
            self._metrics.loc[R.index, "richness_%s" % tag] = (L.M > 0).sum(axis=1)
            METRIC_LABEL["shannon_%s" % tag] = "Shannon diversity, %s %s %s" % (
                LAYER_LABEL[key.split(".")[0]], key.split(".")[1], key.split(".")[2])
            METRIC_LABEL["richness_%s" % tag] = "richness (taxa > 0), %s %s %s" % (
                LAYER_LABEL[key.split(".")[0]], key.split(".")[1], key.split(".")[2])
        self.metric_names = list(self._metrics.columns)

    @property
    def metrics(self):
        return self._metrics

    # ------------------------------------------------------------------ QC
    def flag_samples(self):
        m = self._metrics
        rows = OrderedDict()
        for sid in self.samples.index:
            hits = []
            mg = m.at[sid, "mg_pairs_M"] if "mg_pairs_M" in m.columns else np.nan
            if pd.isna(mg):
                hits.append(("FAIL", "no metagenome reads"))
            elif mg * 1e6 < self.a.min_mg_pairs:
                hits.append(("FAIL", "mg read pairs %.2g < %.2g" % (mg * 1e6, self.a.min_mg_pairs)))
            mt = m.at[sid, "mt_pairs_M"] if "mt_pairs_M" in m.columns else np.nan
            if pd.isna(mt):
                hits.append(("WARN", "no metatranscriptome reads (mt layers without this sample)"))
            elif mt * 1e6 < self.a.min_mt_pairs:
                hits.append(("WARN", "mt read pairs %.2g < %.2g (mt layers without this sample)" % (mt * 1e6, self.a.min_mt_pairs)))
            rows[sid] = hits
        z = self.a.outlier_z
        if z and z > 0:
            for met, min_abs in QC_OUTLIER.items():
                if met not in m.columns:
                    continue
                x = m[met].dropna()
                if len(x) < 10:
                    continue
                med = x.median()
                mad = 1.4826 * (x - med).abs().median()
                if not mad:
                    continue
                for sid, v in x.items():
                    rz = (v - med) / mad
                    if abs(rz) > z and abs(v - med) >= min_abs * (1 if "fraction" in met else abs(med)):
                        rows[sid].append(("WARN", "outlier %s = %.3g (robust z %.1f)" % (met, v, rz)))
        out = []
        for sid, hits in rows.items():
            lv = "FAIL" if any(h[0] == "FAIL" for h in hits) else "WARN" if hits else "PASS"
            out.append({"sample_id": sid, "level": lv, "reasons": "; ".join(h[1] for h in hits)})
        return pd.DataFrame(out).set_index("sample_id")

    def mt_ok(self, sid):
        mt = self._metrics.at[sid, "mt_pairs_M"] if "mt_pairs_M" in self._metrics.columns else np.nan
        return not pd.isna(mt) and mt * 1e6 >= self.a.min_mt_pairs

    # ------------------------------------------------------------------ layers
    def load_layers(self, keys):
        tax = fun = ba = bins = None
        keys = self.expand_metab(keys)
        for key in [k.strip() for k in keys.split(",") if k.strip()]:
            parts = key.split(".")
            try:
                if parts[0] == "metab":
                    self.metab_layer(key, parts[1])
                elif parts[0] == "prot":
                    self.prot_layer(key, parts[1], parts[2])
                elif parts[0] in ("kraken", "motus"):
                    if tax is None:
                        tax = read_table(self.release, "taxonomy")
                    src = "kraken2_bracken" if parts[0] == "kraken" else "motus"
                    d = tax[(tax["source"] == src) & (tax["rank"] == parts[1]) & (tax["omic"] == parts[2])]
                    if src == "motus" and parts[1] == "species" and d.empty:
                        d = self._old_motus(tax, parts[2])
                    M = d.pivot_table(index="sample_id", columns="taxon", values="reads", aggfunc="sum", fill_value=0)
                    lab = "%s %s, %s" % (LAYER_LABEL[parts[0]], parts[1], OMIC[parts[2]])
                    self._add(key, lab, parts[2], M)
                elif parts[0] == "mag":
                    if ba is None:
                        ba, bins = read_table(self.release, "bin_abundance"), read_table(self.release, "bins")
                    rank = parts[1]
                    b = bins[["bin_id", rank]].copy()
                    b[rank] = b[rank].fillna("").replace("", "unclassified")
                    d = ba[ba["omic"] == parts[2]].merge(b, on="bin_id")
                    M = d.pivot_table(index="sample_id", columns=rank, values="reads", aggfunc="sum", fill_value=0)
                    M = M.drop(columns=[c for c in M.columns if c == "unclassified"])
                    self._add(key, "MAG reads per GTDB %s, %s" % (rank, OMIC[parts[2]]), parts[2], M)
                elif parts[0] == "ratio":
                    continue                      # after the count layers
                else:
                    if fun is None:
                        fun = read_table(self.release, "functions")
                    db, omic = parts[0], parts[1]
                    d = fun[(fun["db"] == db) & (fun["omic"] == omic)]
                    if db == "kegg_pathway":
                        d = d[~d["feature"].str.match(GLOBAL_MAPS)]
                    M = d.pivot_table(index="sample_id", columns="feature", values="reads", aggfunc="sum", fill_value=0)
                    if db == "kegg_pathway":          # Mantis reports pathways as koNNNNN and mapNNNNN: keep ko
                        M = M.drop(columns=[c for c in M.columns if c.startswith("map") and "ko" + c[3:] in M.columns])
                    self._add(key, "%s, %s" % (LAYER_LABEL.get(db, db), OMIC[omic]), omic, M, db=db)
            except (KeyError, IndexError) as e:
                self.note("layer %s skipped (%s)" % (key, e))
        for key in [k.strip() for k in keys.split(",") if k.strip().startswith("ratio.")]:
            base = key[len("ratio."):]
            mg, mt = self.layers.get(base + ".mg"), self.layers.get(base + ".mt")
            if mg is None or mt is None:
                self.note("layer %s skipped: needs %s.mg and %s.mt" % (key, base, base))
                continue
            Y = self.ratio(mg, mt)
            if Y is not None:
                self.layers[key] = Layer(key, "mt/mg ratio, %s" % mg.label.rsplit(",", 1)[0], "mt/mg", "ratio", Y,
                                         db=mg.db)
        self.add_alpha()

    # ------------------------------------------------------------------ metabolomics (imp-explorer metabolomics)
    def metab_assays(self):
        if getattr(self, "_metab_assays", None) is None:
            a = read_table(self.release, "metab_assays")
            self._metab_assays = a.set_index("assay") if len(a) else pd.DataFrame()
        return self._metab_assays

    def expand_metab(self, keys):
        """'metab' = every assay of the release: metab.<assay>,...; 'prot' = the default protein layers
        (genus, KEGG KO / module / pathway of every fraction)"""
        out = []
        for k in [k.strip() for k in keys.split(",") if k.strip()]:
            if k in ("metab", "metab.all"):
                out += ["metab." + x for x in self.metab_assays().index]
            elif k in ("prot", "prot.all"):
                s = read_table(self.release, "prot_samples")
                fr = sorted({str(o)[3:] for o in s["omic"]}) if len(s) else []
                pdbs = read_manifest(self.release).get("proteomics", {}).get("dbs") or []
                out += ["prot.%s.%s" % (x, f) for x in ("genus", "kegg_module", "kegg_pathway", "kegg_ko") +
                        tuple(d for d in ("pfam", "cog") if d in pdbs) for f in fr]
            else:
                out.append(k)
        return ",".join(dict.fromkeys(out))

    def prot_layer(self, key, what, fraction):
        """metaproteomics (imp-explorer proteomics): AUC per taxon (rank) or KEGG term of one fraction"""
        omic = "mp_" + fraction
        if what in ("domain", "phylum", "class", "order", "family", "genus", "species"):
            if getattr(self, "_prot_tax", None) is None:
                self._prot_tax = read_table(self.release, "prot_taxonomy")
            d = self._prot_tax
            d = d[(d["rank"] == what) & (d["omic"] == omic)] if len(d) else d
            col, lab = "taxon", what
        else:
            if getattr(self, "_prot_fun", None) is None:
                self._prot_fun = read_table(self.release, "prot_functions")
            d = self._prot_fun
            d = d[(d["db"] == what) & (d["omic"] == omic)] if len(d) else d
            if what == "kegg_pathway" and len(d):
                d = d[~d["feature"].str.match(GLOBAL_MAPS)]
            col, lab = "feature", LAYER_LABEL.get(what, what)
        if d.empty:
            self.note("layer %s skipped: no %s of %s in the release (imp-explorer proteomics)" % (key, what, omic))
            return
        M = d.pivot_table(index="sample_id", columns=col, values="auc", aggfunc="sum", fill_value=0)
        M.index = M.index.astype(str)
        M = M.reindex([s for s in self.kept if s in M.index]).fillna(0)
        M = M.loc[M.sum(axis=1) > 0, M.sum(axis=0) > 0]
        if M.shape[0] < 6 or M.shape[1] < 2:
            self.note("layer %s skipped: %d samples, %d features" % (key, M.shape[0], M.shape[1]))
            return
        frac = {"extra": "extracellular", "intra": "intracellular"}.get(fraction, fraction)
        self.layers[key] = Layer(key, "Metaproteome %s, %s" % (frac, lab), omic, "auc", M.astype(float),
                                 db=None if col == "taxon" else what)

    def metab_layer(self, key, assay):
        if getattr(self, "_metab_values", None) is None:
            self._metab_values = read_table(self.release, "metab_values")
        d = self._metab_values
        d = d[d["assay"].astype(str) == assay]
        if d.empty:
            self.note("layer %s skipped: assay %s not in the release (imp-explorer metabolomics)" % (key, assay))
            return
        M = d.pivot_table(index="sample_id", columns="feature_id", values="value", aggfunc="first")
        M.index = M.index.astype(str)
        M = M.reindex([s for s in self.kept if s in M.index])
        # missing values: features missing in > 20 % of the samples, then samples with a missing value left out
        nf = M.shape[1]
        M = M.loc[:, M.isna().mean() <= 0.2]
        ns = M.shape[0]
        M = M.dropna()
        if M.shape[1] < nf or M.shape[0] < ns:
            self.note("layer %s: %d features (missing in > 20 %% of the samples) and %d samples (a missing value) "
                      "left out" % (key, nf - M.shape[1], ns - M.shape[0]))
        if (M < 0).any().any():
            self.note("layer %s: negative values (already log-transformed?) - analysed as they are" % key)
        if M.shape[0] < 6 or M.shape[1] < 2:
            self.note("layer %s skipped: %d samples, %d features" % (key, M.shape[0], M.shape[1]))
            return
        a = self.metab_assays()
        lab = a["label"].get(assay, assay) if "label" in a.columns else assay
        L = Layer(key, "Metabolomics: %s" % lab, "metab", "metab", M.astype(float), db="metab")
        f = read_table(self.release, "metab_features")
        f = f[f["assay"].astype(str) == assay] if len(f) else f
        L.names = dict(zip(f["feature_id"].astype(str), f["name"].astype(str))) if len(f) and "name" in f else {}
        self.layers[key] = L

    def _old_motus(self, tax, omic):
        """mOTUs of releases exported before imp-explorer 0.1.1: every row has rank 'motu' and the full lineage as
        name; the species-level mOTUs are the rows with an mOTU id -> their species name"""
        d = tax[(tax["source"] == "motus") & (tax["rank"] == "motu") & (tax["omic"] == omic)]
        d = d[d["taxid"].fillna("").astype(str).str.contains("mOTU")].copy()
        if len(d):
            d["taxon"] = (d["taxon"].astype(str).str.split("|").str[-1].str.replace(r"^[a-z]__", "", regex=True)
                          .str.replace(r"\s*\[[^\]]+\]\s*$", "", regex=True).str.strip())
            self.note("mOTUs %s: release from an older export (rank 'motu'); %d species-level rows used — "
                      "re-export for all mOTUs ranks in the app" % (omic, len(d)))
        return d

    def _add(self, key, label, omic, M, db=None):
        M.index = M.index.astype(str)
        keep = [s for s in self.kept if s in M.index and (omic != "mt" or self.mt_ok(s))]
        M = M.reindex(keep).fillna(0)
        M = M.loc[M.sum(axis=1) > 0, M.sum(axis=0) > 0]
        if M.shape[0] < 6 or M.shape[1] < 2:
            self.note("layer %s skipped: %d samples, %d features" % (key, M.shape[0], M.shape[1]))
            return
        self.layers[key] = Layer(key, label, omic, "counts", M.astype(float), db=db)

    def ratio(self, mg, mt):
        """expression ratio per feature: CLR(mt) - CLR(mg) of the features tested in both (log scale, relative to
        the average feature of the sample); > 0 = more transcribed than its gene copies suggest"""
        ss = [s for s in mg.M.index if s in mt.M.index]
        f = sorted(set(mg.tested(self.a.prevalence)) & set(mt.tested(self.a.prevalence)))
        if len(ss) < 6 or len(f) < 3:
            return None
        A = clr(mg.M.loc[ss, f].values)
        B = clr(mt.M.loc[ss, f].values)
        return pd.DataFrame((B - A) / np.log(2), index=ss, columns=f)        # log2 units
