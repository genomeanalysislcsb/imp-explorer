"""Analyses of the cohort: every result is a DataFrame (written as TSV, used by the report and the app)."""
import time

import numpy as np
import pandas as pd

from .data import METRIC_LABEL
from .stats import (TERM_NAMES, bh, bray_curtis, chi2_test, clr, design_matrix, design_ok, euclidean_z, kruskal,
                    log_values, mannwhitney,
                    nmds, ols_hc3, pcoa, permanova, permanova_partial, van_elteren, isolated_samples)

MIN_N = 3          # samples per side for a comparison


class Analysis:
    def __init__(self, C):
        self.C, self.a = C, C.a
        self.res = {}
        self.timing = {}

    # ------------------------------------------------------------------ helpers
    def comparisons(self, samples):
        """[(label, kind, a_ids, b_ids)]: each group vs the reference within the blocks they share
        (kind 'within'; no block column: 'pooled'), and the reference split by block (kind 'block': site effect
        without disease)"""
        C = self.C
        out = []
        if len(C.groups) < 2:
            return out
        g_of = {s: C.group(s) for s in samples}
        b_of = {s: C.block_of(s) for s in samples}
        ref = C.ref
        rb = {b_of[s] for s in samples if g_of[s] == ref} - {"NA"}
        for g in C.groups:
            if g == ref:
                continue
            gb = {b_of[s] for s in samples if g_of[s] == g} - {"NA"}
            sh = sorted(gb & rb) if C.block else ["all"]
            if not sh:
                continue
            A = [s for s in samples if g_of[s] == g and b_of[s] in sh]
            B = [s for s in samples if g_of[s] == ref and b_of[s] in sh]
            if len(A) < MIN_N or len(B) < MIN_N:
                continue
            lab = "%s vs %s" % (g, ref) + (" (%s)" % ", ".join(sh) if C.block and len(sh) <= 3 else
                                           " (%d %s levels)" % (len(sh), C.block) if C.block else "")
            out.append((lab, "within" if C.block else "pooled", A, B))
        if C.block:
            lv = sorted(b for b in rb if sum(1 for s in samples if g_of[s] == ref and b_of[s] == b) >= MIN_N)
            if 2 <= len(lv) <= 6:
                for b2 in lv[1:]:
                    A = [s for s in samples if g_of[s] == ref and b_of[s] == b2]
                    B = [s for s in samples if g_of[s] == ref and b_of[s] == lv[0]]
                    out.append(("%s: %s vs %s" % (ref, b2, lv[0]), "block", A, B))
        return out

    def test_matrix(self, V, samples, layer):
        """per feature (columns of V, samples x features): Kruskal-Wallis across groups and every comparison;
        BH per layer and comparison"""
        C = self.C
        rows = []
        gs = [s for s in samples if C.group(s) in C.groups]
        if len(C.groups) >= 2 and len(gs) >= 6:
            labs = np.array([C.group(s) for s in gs])
            X = V.loc[gs].values
            for j, f in enumerate(V.columns):
                h, p = kruskal([X[labs == g, j] for g in C.groups])
                med = {g: float(np.median(X[labs == g, j])) for g in C.groups}
                rows.append([layer, f, "all groups", "kruskal", len(gs), np.nan,
                             np.nan, np.nan, h, p, "; ".join("%s %.4g" % (g, med[g]) for g in C.groups)])
        for lab, kind, A, B in self.comparisons(samples):
            A = [s for s in A if s in V.index]
            B = [s for s in B if s in V.index]
            if len(A) < MIN_N or len(B) < MIN_N:
                continue
            XA, XB = V.loc[A].values, V.loc[B].values
            ba = [C.block_of(s) for s in A]
            bb = [C.block_of(s) for s in B]
            for j, f in enumerate(V.columns):
                if kind == "within" and len(set(ba) | set(bb)) > 1:
                    auc, p, _ = van_elteren(XA[:, j], ba, XB[:, j], bb)
                else:
                    auc, p = mannwhitney(XA[:, j], XB[:, j])
                rows.append([layer, f, lab, kind, len(A), len(B), float(np.median(XA[:, j])),
                             float(np.median(XB[:, j])), auc, p, ""])
        df = pd.DataFrame(rows, columns=["layer", "feature", "comparison", "kind", "n_a", "n_b", "median_a",
                                         "median_b", "effect", "p", "medians"])
        df["q"] = np.nan
        for _, idx in df.groupby("comparison").groups.items():
            df.loc[idx, "q"] = bh(df.loc[idx, "p"].values)
        df["name"] = df["feature"].map(lambda f: TERM_NAMES.get(f, TERM_NAMES.get(str(f).replace("ko", "map", 1), "")))
        return df

    # ------------------------------------------------------------------ overview / covariates / QC
    def design(self):
        C = self.C
        s = C.samples.index
        d = pd.DataFrame({"group": [C.group(x) for x in s], "block": [C.block_of(x) for x in s],
                          "kept": [x in C.kept for x in s]})
        t = d.groupby(["group", "block"]).agg(samples=("kept", "size"), kept=("kept", "sum")).reset_index()
        self.res["design"] = t

    def covariates(self):
        C = self.C
        rows = []
        ks = C.kept
        groups = C.groups or ["all"]
        comps = self.comparisons(ks)
        for c in C.covariates:
            kind = C.kinds.get(c)
            vals = {s: C.meta(s, c) for s in ks}
            n = sum(1 for v in vals.values() if v != "NA")
            if not n:
                continue
            desc, extra = [], []
            if kind == "numeric":
                bys = [np.array([float(vals[s]) for s in ks if C.group(s) == g and vals[s] != "NA"]) for g in groups]
                _, p = kruskal([x for x in bys if len(x)]) if len(groups) > 1 else (np.nan, np.nan)
                for x in bys:
                    desc.append("%.3g [%.3g–%.3g] (n=%d)" % (np.median(x), np.percentile(x, 25), np.percentile(x, 75),
                                                             len(x)) if len(x) else "–")
                for lab, kd, A, B in comps:
                    xa = [float(vals[s]) for s in A if vals[s] != "NA"]
                    xb = [float(vals[s]) for s in B if vals[s] != "NA"]
                    extra.append(mannwhitney(xa, xb)[1] if len(xa) >= 2 and len(xb) >= 2 else np.nan)
                test = "Kruskal-Wallis"
            else:
                levels = sorted({v for v in vals.values() if v != "NA"})
                if len(levels) > 15:
                    continue
                tab = np.array([[sum(1 for s in ks if C.group(s) == g and vals[s] == l) for l in levels]
                                for g in groups], float)
                p = chi2_test(tab) if len(groups) > 1 else np.nan
                desc = [", ".join("%s %d" % (l, int(k)) for l, k in zip(levels, row) if k) or "–" for row in tab]
                for lab, kd, A, B in comps:
                    t2 = np.array([[sum(1 for s in A if vals[s] == l) for l in levels],
                                   [sum(1 for s in B if vals[s] == l) for l in levels]], float)
                    extra.append(chi2_test(t2))
                test = "chi-square"
            rows.append([c, kind, n, test, p] + desc + extra)
        cols = ["covariate", "kind", "n", "test", "p"] + list(groups) + ["p " + lab for lab, _, _, _ in comps]
        self.res["covariates"] = pd.DataFrame(rows, columns=cols)

    def qc(self):
        C = self.C
        f = C.qc_flags.copy()
        f["group"] = [C.group(s) for s in f.index]
        f["block"] = [C.block_of(s) for s in f.index]
        f["excluded"] = [s in C.excluded for s in f.index]
        self.res["qc_flags"] = f.reset_index()
        rows, tab = [], []
        for g in (C.groups or ["all"]):
            x = f[f["group"] == g]
            rows.append([g, len(x), int((x["level"] == "PASS").sum()), int((x["level"] == "WARN").sum()),
                         int((x["level"] == "FAIL").sum()), int(x["excluded"].sum())])
            tab.append([len(x) - int(x["excluded"].sum()), int(x["excluded"].sum())])
        p = chi2_test(np.array(tab, float)) if len(tab) > 1 else np.nan
        q = pd.DataFrame(rows, columns=["group", "samples", "PASS", "WARN", "FAIL", "excluded"])
        q.attrs["p_excluded"] = p
        self.res["qc_summary"] = q

    # ------------------------------------------------------------------ metrics, alpha diversity
    def metrics(self):
        C = self.C
        M = C.metrics.loc[C.kept]
        rows = []
        for m in M.columns:
            x = M[m].dropna()
            if x.nunique() < 3:
                continue
            t = self.test_matrix(x.to_frame(m), list(x.index), "metrics")
            rows.append(t)
        df = pd.concat(rows) if rows else pd.DataFrame()
        if len(df):                      # BH across metrics within each comparison
            for _, idx in df.groupby("comparison").groups.items():
                df.loc[idx, "q"] = bh(df.loc[idx, "p"].values)
            df["name"] = df["feature"].map(lambda f: METRIC_LABEL.get(f, f))
        self.res["metrics_tests"] = df
        out = M.copy()
        out.insert(0, "block", [C.block_of(s) for s in out.index])
        out.insert(0, "group", [C.group(s) for s in out.index])
        self.res["sample_metrics"] = out.reset_index()

    # ------------------------------------------------------------------ features
    def features(self):
        C = self.C
        parts = []
        for key, L in C.layers.items():
            t0 = time.time()
            f = L.tested(self.a.prevalence)
            if not f:
                continue
            V = L.rel()[f] if L.kind in ("counts", "auc") else L.M[f]
            t = self.test_matrix(V, list(V.index), key)
            t["mean"] = t["feature"].map(V.mean())
            t["prevalence"] = t["feature"].map((L.M[f] > 0).mean()) if L.kind in ("counts", "metab", "auc") else np.nan
            if L.kind == "metab":
                t["name"] = t["feature"].map(getattr(L, "names", {})).fillna("")
            parts.append(t)
            self.timing["tests " + key] = time.time() - t0
        self.res["feature_tests"] = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()

    # ------------------------------------------------------------------ beta diversity, dispersion, ordination
    def beta(self):
        C = self.C
        perm, disp, ords = [], [], []
        nmds_layers = {k.strip() for k in (self.a.nmds or "").split(",") if k.strip()}
        P = self.a.permutations
        for key, L in C.layers.items():
            if L.kind not in ("counts", "metab", "auc"):
                continue
            t0 = time.time()
            R = L.rel()
            ss = list(R.index)
            # counts: Bray-Curtis of relative abundances; metabolomics: Euclidean of autoscaled log values
            D = bray_curtis(R.values) if L.kind in ("counts", "auc") else euclidean_z(R[L.tested(self.a.prevalence) or list(R.columns)].values)
            pos = {s: i for i, s in enumerate(ss)}
            X, expl = pcoa(D, 2)
            o = pd.DataFrame({"layer": key, "sample_id": ss, "axis1": X[:, 0], "axis2": X[:, 1], "method": "PCoA",
                              "group": [C.group(s) for s in ss], "block": [C.block_of(s) for s in ss]})
            o.attrs["explained"] = list(map(float, expl))
            o["explained1"], o["explained2"] = (float(expl[0]) if len(expl) else np.nan,
                                                float(expl[1]) if len(expl) > 1 else np.nan)
            ords.append(o)
            if key in nmds_layers:
                iso = isolated_samples(D)
                keep = [i for i in range(len(ss)) if i not in iso]
                Xn, stress = nmds(D[np.ix_(keep, keep)], k=2, starts=5)
                if Xn is not None:
                    ords.append(pd.DataFrame({"layer": key, "sample_id": [ss[i] for i in keep], "axis1": Xn[:, 0],
                                              "axis2": Xn[:, 1], "method": "NMDS",
                                              "group": [C.group(ss[i]) for i in keep],
                                              "block": [C.block_of(ss[i]) for i in keep],
                                              "explained1": stress, "explained2": np.nan}))
            gs = [s for s in ss if C.group(s) in C.groups]
            if len(C.groups) >= 2 and len(gs) >= 6:
                idx = [pos[s] for s in gs]
                Dg = D[np.ix_(idx, idx)]
                lab = np.array([C.group(s) for s in gs], object)
                blk = np.array([C.block_of(s) for s in gs], object) if C.block else None
                F, r2, p = permanova(Dg, lab, P, strata=blk)
                perm.append([key, "all groups" + (" (permuted within %s)" % C.block if C.block else ""), "groups",
                             len(gs), F, r2, p])
                if C.block and len(set(blk)) > 1:
                    F, r2, p = permanova_partial(Dg, lab, blk, P)
                    perm.append([key, "groups | %s (partial)" % C.block, "partial", len(gs), F, r2, p])
                    F, r2, p = permanova_partial(Dg, blk, lab, P)
                    perm.append([key, "%s | groups (partial)" % C.block, "partial", len(gs), F, r2, p])
            for lab_, kind, A, B in self.comparisons(ss):
                ids = A + B
                idx = [pos[s] for s in ids]
                Dc = D[np.ix_(idx, idx)]
                lab = np.array(["a"] * len(A) + ["b"] * len(B), object)
                blk = np.array([C.block_of(s) for s in ids], object)
                F, r2, p = permanova(Dc, lab, P, strata=blk if kind == "within" and len(set(blk)) > 1 else None)
                perm.append([key, lab_, kind, len(ids), F, r2, p])
                # dispersion: distance to the own centroid in the full PCoA space of the comparison
                Xf, _ = pcoa(Dc, k=len(ids) - 1)
                d = np.zeros(len(ids))
                for g in ("a", "b"):
                    i = lab == g
                    d[i] = np.sqrt(((Xf[i] - Xf[i].mean(axis=0)) ** 2).sum(axis=1))
                auc, p = mannwhitney(d[lab == "a"], d[lab == "b"])
                ma, mb = float(np.median(d[lab == "a"])), float(np.median(d[lab == "b"]))
                disp.append([key, lab_, kind, len(A), len(B), ma, mb, ma / mb if mb else np.nan, p])
            self.timing["beta " + key] = time.time() - t0
        pm = pd.DataFrame(perm, columns=["layer", "comparison", "kind", "n", "F", "R2", "p"])
        pm["q"] = bh(pm["p"].values) if len(pm) else []
        self.res["permanova"] = pm
        ds = pd.DataFrame(disp, columns=["layer", "comparison", "kind", "n_a", "n_b", "median_dist_a",
                                         "median_dist_b", "ratio", "p"])
        ds["q"] = bh(ds["p"].values) if len(ds) else []
        self.res["dispersion"] = ds
        self.res["ordination"] = pd.concat(ords, ignore_index=True) if ords else pd.DataFrame()

    # ------------------------------------------------------------------ adjusted models
    def model_specs(self):
        a, C = self.a, self.C
        forms = [m.strip() for m in (a.model or "").split(",") if m.strip()]
        if not forms and C.group_col:
            forms = [C.group_col + ("+" + C.block if C.block else "")]
            if C.covariates:
                extra = [c for c in C.covariates if c != C.block and C.kinds.get(c) in ("numeric", "categorical")]
                if extra and len(extra) <= 4:
                    forms.append(forms[0] + "+" + "+".join(extra))
        notes = {}
        for part in (a.model_notes or "").split(";"):
            if "=" in part:
                k, v = part.split("=", 1)
                notes[k.strip()] = v.strip()
        out = []
        for i, f in enumerate(forms):
            terms = [t.strip() for t in f.split("|", 1)[0].split("+") if t.strip()]
            subset = None
            if "|" in f and "=" in f.split("|", 1)[1]:
                col, vals = f.split("|", 1)[1].split("=", 1)
                subset = (col.strip(), {v.strip() for v in vals.split("/") if v.strip()})
            out.append({"model": "M%d" % (i + 1), "formula": f, "terms": terms, "subset": subset,
                        "focus": terms[0] if terms else None, "note": notes.get("M%d" % (i + 1), "")})
        return out

    def model_design(self, spec, samples):
        C = self.C
        kinds, use = {}, []
        for t in spec["terms"]:
            if t == C.group_col:
                kinds[t] = "categorical"
            elif t in C.meta_cols:
                kinds[t] = C.kinds.get(t, "categorical")
                if kinds[t] == "empty":
                    continue
            elif t in C.metric_names:
                kinds[t] = "numeric"
            else:
                continue
            use.append(t)
        rows, keys = [], []
        for s in samples:
            if spec["subset"]:
                col, vals = spec["subset"]
                if C.meta(s, col) not in vals:
                    continue
            r = {t: C.meta(s, t) for t in use}
            if all(v not in ("NA", "") for v in r.values()):
                rows.append(r)
                keys.append(s)
        if len(keys) < max(8, len(use) + 4) or not use:
            return None
        ref = {C.group_col: C.ref} if C.group_col else {}
        X, names, owner, info = design_matrix(rows, use, kinds, ref)
        if len(keys) - X.shape[1] < 3 or X.shape[1] < 2:
            return None
        missing = [t for t in spec["terms"] if t not in use]
        if missing:
            info.append("not in the metadata / metrics: %s" % ", ".join(missing))
        return keys, X, names, owner, info

    def fit_set(self, X, Y, feats, names, owner, spec, layer, n_rows):
        out = []
        if not np.isnan(Y).any():
            if not design_ok(X):
                return out
            B, SE, T, P, df = ols_hc3(X, Y)
            for c in range(1, X.shape[1]):
                for j, f in enumerate(feats):
                    out.append([spec["model"], layer, f, X.shape[0], names[c], owner[c] == spec["focus"],
                                B[c, j], SE[c, j], T[c, j], P[c, j]])
        else:                                       # per feature on its complete cases (metrics)
            for j, f in enumerate(feats):
                ok = ~np.isnan(Y[:, j])
                if ok.sum() < X.shape[1] + 3 or np.nanstd(Y[ok, j]) == 0 or len(np.unique(Y[ok, j])) < 5:
                    continue
                Xo = X[ok]
                if np.linalg.matrix_rank(Xo) < Xo.shape[1] or not design_ok(Xo):
                    continue
                B, SE, T, P, df = ols_hc3(Xo, Y[ok, j][:, None])
                for c in range(1, Xo.shape[1]):
                    out.append([spec["model"], layer, f, int(ok.sum()), names[c], owner[c] == spec["focus"],
                                B[c, 0], SE[c, 0], T[c, 0], P[c, 0]])
        return out

    def models(self):
        C = self.C
        specs = self.model_specs()
        rows, over, designs = [], [], []
        for spec in specs:
            t0 = time.time()
            sets = [("metrics", None)] + list(C.layers.items())
            for key, L in sets:
                samples = C.kept if L is None else list(L.M.index)
                md = self.model_design(spec, samples)
                if md is None:
                    over.append([spec["model"], spec["formula"], spec["note"], key, 0, 0, 0, "too few samples"])
                    continue
                keys, X, names, owner, info = md
                if L is None:
                    feats = [m for m in C.metric_names if m not in spec["terms"]]
                    Y = C.metrics.loc[keys, feats].values.astype(float)
                elif L.kind == "counts":
                    sub = L.M.loc[keys]
                    keepf = (sub >= 10).mean(axis=0) >= self.a.prevalence
                    feats = list(sub.columns[keepf])
                    if len(feats) < 2:
                        continue
                    Y = clr(sub[feats].values)
                elif L.kind == "auc":
                    # proteins: CLR of the AUC shares, pseudocount half the smallest positive share of the layer
                    sub = L.rel().loc[keys]
                    keepf = (sub > 0).mean(axis=0) >= self.a.prevalence
                    feats = list(sub.columns[keepf])
                    if len(feats) < 2:
                        continue
                    pos = sub[feats].values[sub[feats].values > 0]
                    Y = clr(sub[feats].values, pseudo=pos.min() / 2 if len(pos) else 1e-9)
                elif L.kind == "metab":
                    sub = L.M.loc[keys]
                    keepf = ((sub > 0).mean(axis=0) >= self.a.prevalence) & (sub.std(axis=0) > 0)
                    feats = list(sub.columns[keepf])
                    if len(feats) < 2:
                        continue
                    Y = log_values(sub[feats].values)          # log2 units
                else:
                    feats = list(L.M.columns)
                    Y = L.M.loc[keys, feats].values
                r = self.fit_set(X, Y, feats, names, owner, spec, key, len(keys))
                for term in sorted({x[4] for x in r}):
                    sel = [x for x in r if x[4] == term]
                    for x, q in zip(sel, bh([x[9] for x in sel])):
                        x.append(q)
                rows += r
                nsig = len({x[2] for x in r if x[5] and x[10] < self.a.sig_q})
                fitted = len({x[2] for x in r})
                if not fitted:
                    info = info + ["not fitted: fewer than 20 samples, < 3 samples in a level or a sample with "
                                   "leverage > 0.95"]
                over.append([spec["model"], spec["formula"], spec["note"], key, len(keys), fitted, nsig,
                             "; ".join(info)])
                if key == "metrics":
                    designs.append([spec["model"], spec["formula"], len(keys), ", ".join(names[1:]),
                                    "; ".join(info)])
            self.timing["model " + spec["model"]] = time.time() - t0
        df = pd.DataFrame(rows, columns=["model", "layer", "feature", "n", "term", "focus", "estimate", "se", "t",
                                         "p", "q"])
        df["name"] = [METRIC_LABEL.get(f, "") if l == "metrics" else
                      TERM_NAMES.get(f, TERM_NAMES.get(str(f).replace("ko", "map", 1), ""))
                      for l, f in zip(df["layer"], df["feature"])]
        for key, L in C.layers.items():
            if L.kind == "metab" and getattr(L, "names", None):
                i = df["layer"] == key
                df.loc[i, "name"] = df.loc[i, "feature"].map(L.names).fillna("")
        df = df.sort_values(["model", "layer", "term", "p"])
        self.res["models"] = df
        self.res["models_overview"] = pd.DataFrame(over, columns=["model", "formula", "note", "layer", "n",
                                                                  "features", "features_q_focus", "info"])
        self.res["model_designs"] = pd.DataFrame(designs, columns=["model", "formula", "n", "columns", "info"])
        self.res["model_specs"] = specs

    # ------------------------------------------------------------------ MAGs
    def mags(self):
        """per GTDB genus: share of samples with a medium/high-quality MAG, by group; each group vs the reference
        within blocks (van Elteren on presence)"""
        C = self.C
        from ..release import read_table
        bins = read_table(C.release, "bins")
        if not len(bins):
            self.res["mag_prevalence"] = pd.DataFrame()
            return
        bins = bins[bins["sample_id"].astype(str).isin(C.kept)]
        bins = bins[bins["genus"].fillna("") != ""]
        pres = bins.groupby(["genus", "sample_id"]).size().unstack(fill_value=0).T.reindex(C.kept).fillna(0) > 0
        phylum = bins.groupby("genus")["phylum"].first()
        keep = pres.columns[pres.mean() >= 0.05]
        rows = []
        comps = self.comparisons(C.kept)
        for gname in keep:
            x = pres[gname].astype(float)
            r = [gname, phylum.get(gname, ""), int(bins[bins["genus"] == gname].shape[0]), float(x.mean())]
            for g in (C.groups or ["all"]):
                sel = [s for s in C.kept if C.group(s) == g]
                r.append(float(x.loc[sel].mean()) if sel else np.nan)
            for lab, kind, A, B in comps:
                if kind == "within" and C.block:
                    auc, p, _ = van_elteren(x.loc[A].values, [C.block_of(s) for s in A], x.loc[B].values,
                                            [C.block_of(s) for s in B])
                else:
                    auc, p = mannwhitney(x.loc[A].values, x.loc[B].values)
                r.append(p)
            rows.append(r)
        cols = ["genus", "phylum", "mags", "prevalence"] + ["prev " + g for g in (C.groups or ["all"])] + \
               ["p " + lab for lab, _, _, _ in comps]
        df = pd.DataFrame(rows, columns=cols)
        for lab, _, _, _ in comps:
            df["q " + lab] = bh(df["p " + lab].values) if len(df) else []
        self.res["mag_prevalence"] = df.sort_values("prevalence", ascending=False)

    # ------------------------------------------------------------------ all
    def run(self):
        steps = [("design", self.design), ("covariates", self.covariates), ("qc", self.qc),
                 ("metrics", self.metrics), ("features", self.features), ("beta", self.beta),
                 ("models", self.models), ("mags", self.mags)]
        for name, fn in steps:
            t0 = time.time()
            fn()
            self.timing[name] = time.time() - t0
            import sys
            sys.stderr.write("[cohort] %-10s %.1f s\n" % (name, self.timing[name]))
        self.res["findings"] = self.findings()
        return self.res

    # ------------------------------------------------------------------ key findings
    def findings(self):
        """[(area, text)] of the results that pass the thresholds; a guide to the report, not a result list"""
        C, R, q = self.C, self.res, self.a.sig_q
        out = []
        des = R["design"]
        if C.block and C.groups:
            only = []
            for g in C.groups:
                bl = sorted(set(des.loc[(des["group"] == g) & (des["kept"] > 0), "block"]))
                if len(bl) == 1:
                    only.append("%s only at %s" % (g, bl[0]))
            if only:
                out.append(("design", "nested design: %s; group effects are estimated within the shared %s, and the "
                            "%s contrast within %s shows how large %s differences alone are" %
                            ("; ".join(only), C.block, C.block, C.ref, C.block)))
        qs = R["qc_summary"]
        out.append(("QC", "%d of %d samples used (%d excluded: %s); exclusion vs group χ² p = %s" % (
            len(C.kept), len(C.samples), len(C.excluded), ", ".join(C.excluded[:8]) or "none",
            "%.2g" % qs.attrs.get("p_excluded") if qs.attrs.get("p_excluded") == qs.attrs.get("p_excluded") else "–")))
        cv = R["covariates"]
        if len(cv):
            for lab in [c for c in cv.columns if c.startswith("p ")]:
                bad = cv[(cv[lab] < 0.05) & (cv["covariate"] != C.block)]
                if len(bad):
                    out.append(("covariates", "%s differ in %s (p < 0.05): adjust for them in the models" % (
                        ", ".join(bad["covariate"]), lab[2:])))
        mt = R["metrics_tests"]
        if len(mt):
            for comp, d in mt[(mt["kind"] != "kruskal") & (mt["q"] < q)].groupby("comparison"):
                out.append(("metrics", "%s: %s" % (comp, ", ".join(
                    "%s (%s)" % (r.feature, "higher" if r.effect > 0.5 else "lower") for r in d.itertuples()))))
        pm = R["permanova"]
        if len(pm):
            for r in pm[pm["p"] < 0.05].itertuples():
                if r.kind in ("within", "block", "groups"):
                    out.append(("beta diversity", "%s · %s: R² = %.3f, p = %.3g" % (C.layers[r.layer].label,
                                                                                   r.comparison, r.R2, r.p)))
        ft = R["feature_tests"]
        if len(ft):
            for (lay, comp), d in ft[(ft["kind"] != "kruskal") & (ft["q"] < q)].groupby(["layer", "comparison"]):
                top = d.sort_values("p").head(5)
                out.append(("features", "%s · %s: %d feature(s) q < %g, e.g. %s" % (
                    C.layers[lay].label, comp, len(d), q, ", ".join(
                        "%s (%s)" % (f, "↑" if e > 0.5 else "↓") for f, e in zip(top["feature"], top["effect"])))))
        mo = R["models"]
        if len(mo):
            for (m, lay), d in mo[mo["focus"] & (mo["q"] < q)].groupby(["model", "layer"]):
                top = d.sort_values("p").head(4)
                lab = "library metrics" if lay == "metrics" else C.layers[lay].label
                out.append(("models", "%s · %s: %d feature(s) with %s q < %g, e.g. %s" % (
                    m, lab, d["feature"].nunique(), "/".join(sorted(d["term"].unique())), q, ", ".join(
                        "%s (%+.2f)" % (f, e) for f, e in zip(top["feature"], top["estimate"])))))
        ds = R["dispersion"]
        if len(ds):
            for r in ds[ds["q"] < q].itertuples():
                out.append(("heterogeneity", "%s · %s: dispersion ratio %.2f (q = %.2g)" % (
                    C.layers[r.layer].label, r.comparison, r.ratio, r.q)))
        mp = R.get("mag_prevalence")
        if mp is not None and len(mp):
            for c in [c for c in mp.columns if c.startswith("q ")]:
                d = mp[mp[c] < q]
                if len(d):
                    out.append(("MAGs", "%s: MAG recovery differs for %s" % (c[2:], ", ".join(d["genus"].head(8)))))
        return out
