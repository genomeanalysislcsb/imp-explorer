"""Data access for the app: DuckDB views over the release tables (Parquet or TSV.gz), sample metadata in memory."""
import os
import threading

import duckdb
import pandas as pd

from ..release import read_manifest, table_path

CORE_TABLES = ("samples", "qc", "taxonomy", "bins", "bin_abundance", "functions")
COHORT_TABLES = ("cohort_models", "cohort_models_overview", "cohort_tests", "cohort_metric_tests", "cohort_permanova",
                 "cohort_dispersion", "cohort_ordination", "cohort_samples", "cohort_covariates", "cohort_qc",
                 "cohort_mag_prevalence")     # written by 'imp-explorer cohort'; optional
HIDDEN_COLS = {"imp_dir"}


class Data:
    def __init__(self, release, extra_tables=()):
        self.release = os.path.abspath(release)
        self.manifest = read_manifest(self.release)
        self.study = self.manifest.get("study") or os.path.basename(self.release)
        self._local = threading.local()
        self.tables = {}
        for name in tuple(CORE_TABLES) + COHORT_TABLES + tuple(extra_tables):
            p = table_path(self.release, name)
            if p:
                self.tables[name] = p
        con = duckdb.connect()
        for name, p in self.tables.items():
            con.execute("CREATE VIEW %s AS SELECT * FROM %s" % (name, self._reader(p)))
        self._con = con
        self.samples = self.q("SELECT * FROM samples") if "samples" in self.tables else pd.DataFrame(
            columns=["sample_id"])
        self.samples["sample_id"] = self.samples["sample_id"].astype(str)
        self.meta_cols = [c for c in self.samples.columns if c not in HIDDEN_COLS | {"sample_id", "has_mg", "has_mt"}]
        # columns usable for grouping: few distinct values
        self.group_cols = [c for c in self.meta_cols if 1 < self.samples[c].nunique(dropna=True) <= 12]

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

    def group_of(self, group_col, samples=None):
        s = self.samples if samples is None else samples
        if not group_col or group_col not in s.columns:
            return pd.Series("all", index=s["sample_id"])
        return pd.Series(s[group_col].fillna("NA").astype(str).values, index=s["sample_id"])
