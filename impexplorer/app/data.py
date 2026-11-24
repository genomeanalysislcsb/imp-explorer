"""Data access for the app: DuckDB views over the release tables (Parquet or TSV.gz), sample metadata in memory."""
import os
import threading

import duckdb
import pandas as pd

from ..release import read_manifest, table_path

CORE_TABLES = ("samples", "qc", "taxonomy", "bins", "bin_abundance", "functions")
COHORT_TABLES = ("cohort_models", "cohort_models_overview", "cohort_tests", "cohort_metric_tests", "cohort_permanova",
                 "cohort_dispersion", "cohort_ordination", "cohort_samples", "cohort_covariates", "cohort_qc",
                 "cohort_mag_prevalence", "cohort_paper_anova", "cohort_paper_permanova",
                 "cohort_paper_variance", "cohort_concord_presence", "cohort_concord_alpha",
                 "cohort_concord_mantel", "cohort_concord_features", "cohort_concord_within",
                 "cohort_fraction_diff", "cohort_fraction_summary")   # written by 'imp-explorer cohort'; optional
SRNA_TABLES = ("srna_samples", "srna_qc", "srna_features", "srna_tests", "srna_models", "srna_models_overview",
               "srna_permanova", "srna_dispersion", "srna_ordination", "srna_covariates")   # 'imp-explorer srna'
PROT_TABLES = ("prot_taxonomy", "prot_functions", "prot_samples", "prot_proteins", "prot_abundance")  # 'proteomics'
METAB_TABLES = ("metab_assays", "metab_features", "metab_values", "metab_samples")    # 'imp-explorer metabolomics'
MODULE_TABLES = ("modules", "module_info")                                             # 'imp-explorer modules'
NETWORK_TABLES = ("met_topology", "met_hubs", "met_mag_coverage", "met_pan_summary", "met_community_reactions",
                  "met_reactions", "met_keystones", "met_stats",
                  "met_taxon_activity")                                                       # 'imp-explorer networks'
QC_TABLES = ("sample_qc",)                                                            # 'imp-explorer qc'
# views over the cohort results of IMP3 and of sRNA-IMP (same columns where they mean the same)
UNION_VIEWS = {"all_tests": ("cohort_tests", "srna_tests"), "all_models": ("cohort_models", "srna_models"),
               "all_models_overview": ("cohort_models_overview", "srna_models_overview"),
               "all_permanova": ("cohort_permanova", "srna_permanova"),
               "all_dispersion": ("cohort_dispersion", "srna_dispersion"),
               "all_ordination": ("cohort_ordination", "srna_ordination")}
HIDDEN_COLS = {"imp_dir"}


class Data:
    def __init__(self, release, extra_tables=()):
        self.release = os.path.abspath(release)
        self.manifest = read_manifest(self.release)
        self.study = self.manifest.get("study") or os.path.basename(self.release)
        self._local = threading.local()
        self.tables = {}
        for name in (tuple(CORE_TABLES) + COHORT_TABLES + SRNA_TABLES + METAB_TABLES + PROT_TABLES + QC_TABLES + MODULE_TABLES + NETWORK_TABLES
                     + tuple(extra_tables)):
            p = table_path(self.release, name)
            if p:
                self.tables[name] = p
        con = duckdb.connect()
        for name, p in self.tables.items():
            con.execute("CREATE VIEW %s AS SELECT * FROM %s" % (name, self._reader(p)))
        # metaproteomics next to mg / mt: the protein profiles join the taxonomy and functions views as omics
        # mp_<fraction> (source 'metaproteomics'; reads = AUC, rel_abundance = share of the total AUC)
        for core, prot, cols in (("taxonomy", "prot_taxonomy",
                                  "sample_id, omic, 'metaproteomics' AS source, rank, taxon, NULL AS taxid, auc AS reads, "
                                  "rel AS rel_abundance"),
                                 ("functions", "prot_functions", "sample_id, omic, db, feature, auc AS reads")):
            if prot in self.tables:
                parts = ["SELECT %s FROM %s" % (cols, self._reader(self.tables[prot]))]
                if core in self.tables:
                    parts.insert(0, "SELECT * FROM %s" % self._reader(self.tables[core]))
                con.execute("CREATE OR REPLACE VIEW %s AS %s" % (core, " UNION ALL BY NAME ".join(parts)))
                self.tables.setdefault(core, "(view)")
        for view, parts in UNION_VIEWS.items():
            have = [t for t in parts if t in self.tables]
            if have:
                con.execute("CREATE VIEW %s AS %s" % (view, " UNION ALL BY NAME ".join(
                    "SELECT * FROM %s" % t for t in have)))
                self.tables[view] = "(view)"
        self._con = con
        self.samples = self.q("SELECT * FROM samples") if "samples" in self.tables else pd.DataFrame(
            columns=["sample_id"])
        self.samples["sample_id"] = self.samples["sample_id"].astype(str)
        self.meta_cols = [c for c in self.samples.columns if c not in HIDDEN_COLS | {"sample_id", "has_mg", "has_mt"}
                          and not str(c).startswith("has_")]
        # columns usable for grouping: few distinct values
        self.group_cols = [c for c in self.meta_cols if 1 < self.samples[c].nunique(dropna=True) <= 12]
        # filter: columns that can select a subset (identifiers such as run accessions, unique per sample, left out)
        self.filter_cols = [c for c in self.meta_cols
                            if 1 < self.samples[c].nunique(dropna=True) < max(2, self.samples[c].notna().sum())]
        # QC status per sample and layer ('imp-explorer qc'); {layer: {sample_id: (status, reasons)}}
        self.sample_qc = {}
        if "sample_qc" in self.tables:
            for r in self.q("SELECT sample_id, layer, status, reasons FROM sample_qc").itertuples(index=False):
                self.sample_qc.setdefault(str(r.layer), {})[str(r.sample_id)] = (str(r.status), r.reasons or "")

    @staticmethod
    def _reader(p):
        p = p.replace("'", "''")
        if p.endswith(".parquet"):
            return "read_parquet('%s')" % p
        return "read_csv_auto('%s', delim='\t', header=true, sample_size=-1)" % p

    def con(self):
        c = getattr(self._local, "con", None)
        if c is None:
            c = self._con.cursor()
            self._local.con = c
        return c

    def q(self, sql, params=None):
        return self.con().execute(sql, params or []).df()

    def has(self, table):
        if table not in self.tables:
            return False
        return bool(self.q("SELECT count(*) AS n FROM %s" % table)["n"].iloc[0])

    def distinct(self, table, col, where="", params=None):
        if table not in self.tables:
            return []
        w = (" WHERE " + where) if where else ""
        return [x for x in self.q("SELECT DISTINCT %s AS v FROM %s%s ORDER BY 1" % (col, table, w), params)["v"]
                if x is not None and x == x]

    # ------------------------------------------------------------------ samples / groups
    def subset(self, filter_col=None, filter_vals=None):
        s = self.samples
        if filter_col and filter_vals:
            s = s[s[filter_col].astype(str).isin([str(v) for v in filter_vals])]
        return s

    def qc_drop(self, layer, mode):
        """sample ids to leave out of `layer` for QC mode 'all' (none), 'nofail' (FAIL) or 'pass' (FAIL + WARN);
        samples without data in the layer ('missing') are not in its tables anyway"""
        if not mode or mode == "all" or not layer:
            return set()
        drop = ("FAIL",) if mode == "nofail" else ("FAIL", "WARN")
        return {s for s, (st, _) in self.sample_qc.get(str(layer), {}).items() if st in drop}

    def qc_layer(self, layer):
        """the QC layer of an omic / fraction / assay name as used in the tables (mg, mt, srna_int, mp_extra,
        metab.<assay>); None when the release has no QC for it"""
        if layer is None:
            return None
        layer = str(layer)
        return layer if layer in self.sample_qc else None

    def group_of(self, group_col, samples=None):
        s = self.samples if samples is None else samples
        if not group_col or group_col not in s.columns:
            return pd.Series("all", index=s["sample_id"])
        return pd.Series(s[group_col].fillna("NA").astype(str).values, index=s["sample_id"])
