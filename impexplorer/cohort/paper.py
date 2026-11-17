"""Paper analysis of the cohort command (--paper-model): one multivariable linear model per layer, as in a paper's
Methods, next to the default analysis (group tests within blocks, OLS on CLR with HC3).

Per layer (the same model formula for every layer; complete cases of its terms):

  alpha     Shannon and Pielou evenness of the layer (counts and protein AUC layers) ~ model; ANOVA per term
  features  every tested feature ~ model; ANOVA per term, BH per layer and term
            response: relative abundance (sum scaling, --paper-transform tss; default for count layers) or CLR
            (clr), CLR for protein AUC and metabolomics (pseudocount: the smallest positive value of the layer)
  beta      PERMANOVA with all terms of the model, marginal sums of squares per term (like adonis2(by = "margin")):
            Bray-Curtis of relative abundances for count layers, Euclidean distance of CLR values (Aitchison) for
            protein AUC and metabolomics; samples permuted (--permutations)

ANOVA: sums of squares of each term given all others (type II; with main effects only the same as type III), F test
on the residual degrees of freedom; r2 = SS(term) / SS(total), the share of the variance explained by the term.
The term 'depth' is log10 of the total reads (or AUC) of the layer's sample, for count and protein layers. Samples
flagged FAIL (or WARN, --paper-qc) in the layer's sample_qc ('imp-explorer qc') are left out of that layer.

Tables: paper_anova (layer, analysis, feature, term, n, df, ss, r2, F, p, q), paper_permanova (layer, term, n,
distance, df, ss, r2, F, p, q), paper_variance (per layer and term: r2 of Shannon and of the PERMANOVA, median r2
over the features, share of features with q < sig_q; plus the residual share).
"""
import time

import numpy as np
import pandas as pd

from ..release import read_table, table_path
from .stats import bh, bray_curtis, clr, design_matrix

try:
    from scipy.special import fdtrc as _fdtrc

    def f_p(F, d1, d2):
        F = np.asarray(F, float)
        return np.where(np.isfinite(F), _fdtrc(d1, d2, np.nan_to_num(F)), np.nan)
except ImportError:                                     # pragma: no cover
    from .stats import _betai

    def f_p(F, d1, d2):
        return np.vectorize(lambda f: _betai(d2 / 2.0, d1 / 2.0, d2 / (d2 + d1 * f)) if np.isfinite(f) and f >= 0
                            else np.nan)(np.asarray(F, float))

MIN_N = 20


def terms_of(formula):
    return [t.strip() for t in str(formula or "").split("+") if t.strip()]


def qc_layer(L):
    """layer name of the sample_qc table for a cohort layer"""
    if L.kind == "metab":
        return L.key
    return L.omic if L.omic not in (None, "mt/mg") else None


def qc_drop(release, mode):
    """{qc layer: set(sample ids)} to leave out; mode FAIL | WARN | none"""
    if mode == "none" or not table_path(release, "sample_qc"):
        return {}
    t = read_table(release, "sample_qc")
    bad = ("FAIL",) if mode == "FAIL" else ("FAIL", "WARN")
    t = t[t["status"].isin(bad)]
    return {str(k): set(v["sample_id"].astype(str)) for k, v in t.groupby("layer")}


def rss(X, Y):
    B, *_ = np.linalg.lstsq(X, Y, rcond=None)
    E = Y - X @ B
    return (E ** 2).sum(axis=0)


def anova(X, owner, Y):
    """type II ANOVA of every column of Y on the design X (columns owned by terms): {term: (df, SS)},
    residual (df, SS), SS total (all arrays over the columns of Y)"""
    n = X.shape[0]
    Yc = Y - Y.mean(axis=0)
    sst = (Yc ** 2).sum(axis=0)
    full = rss(X, Y)
    df_res = n - np.linalg.matrix_rank(X)
    out = {}
    for t in [o for o in dict.fromkeys(owner) if o != "(intercept)"]:
        keep = [j for j, o in enumerate(owner) if o != t]
        Xr = X[:, keep]
        df = np.linalg.matrix_rank(X) - np.linalg.matrix_rank(Xr)
        if df < 1:
            continue
        out[t] = (df, np.maximum(rss(Xr, Y) - full, 0))
    return out, (df_res, full), sst


def anova_rows(layer, analysis, feats, n, terms, res, sst):
    df_res, ss_res = res
    rows = []
    for t, (df, ss) in terms.items():
        with np.errstate(divide="ignore", invalid="ignore"):
            F = (ss / df) / (ss_res / df_res)
            r2 = np.where(sst > 0, ss / sst, np.nan)
        P = f_p(F, df, df_res)
        for j, f in enumerate(feats):
            rows.append([layer, analysis, f, t, n, df, ss[j], r2[j], F[j], P[j]])
    with np.errstate(divide="ignore", invalid="ignore"):
        r2r = np.where(sst > 0, ss_res / sst, np.nan)
    for j, f in enumerate(feats):
        rows.append([layer, analysis, f, "(residual)", n, df_res, ss_res[j], r2r[j], np.nan, np.nan])
    return rows


def hat(X):
    return X @ np.linalg.pinv(X.T @ X) @ X.T


def permanova_margin(D, X, owner, perms=999, seed=1):
    """multivariable PERMANOVA with marginal sums of squares per term (Gower-centred distances, like
    vegan::adonis2(by = 'margin')); samples permuted. Returns [(term, df, SS, R2, F, p)], (df_res, SS_res), SST."""
    n = D.shape[0]
    A = -0.5 * D ** 2
    J = np.eye(n) - 1.0 / n
    G = J @ A @ J
    sst = float(np.trace(G))
    H = hat(X)
    rank = np.linalg.matrix_rank(X)
    df_res = n - rank
    tr = lambda M, Gm: float(np.sum(M * Gm))           # trace(M G) for symmetric M, G
    terms = []
    for t in [o for o in dict.fromkeys(owner) if o != "(intercept)"]:
        keep = [j for j, o in enumerate(owner) if o != t]
        Xr = X[:, keep]
        df = rank - np.linalg.matrix_rank(Xr)
        if df >= 1:
            terms.append((t, df, H - hat(Xr)))
    if not terms or df_res < 1 or sst <= 0:
        return [], (df_res, float("nan")), sst

    def stats(Gm):
        res = sst - tr(H, Gm)
        return res, [tr(Hd, Gm) for _, _, Hd in terms]
    res, ss = stats(G)
    F0 = np.array([(s / df) / (res / df_res) for s, (_, df, _) in zip(ss, terms)])
    rng = np.random.default_rng(seed)
    ge = np.zeros(len(terms))
    for _ in range(perms):
        p = rng.permutation(n)
        rp, sp = stats(G[np.ix_(p, p)])
        Fp = np.array([(s / df) / (rp / df_res) for s, (_, df, _) in zip(sp, terms)])
        ge += Fp >= F0 - 1e-12
    out = [(t, df, s, s / sst, f, (g + 1) / (perms + 1)) for (t, df, _), s, f, g in zip(terms, ss, F0, ge)]
    return out, (df_res, res), sst


def response(L, transform, prevalence):
    """(samples x features DataFrame of the response, feature list) of a layer"""
    feats = L.tested(prevalence)
    if L.kind == "counts":
        R = L.rel()[feats]
        if transform == "clr":
            return pd.DataFrame(clr(L.M[feats].values, pseudo=1.0), index=L.M.index, columns=feats)
        return R
    V = L.M[feats].astype(float)
    pos = V.values[V.values > 0]
    pseudo = pos.min() if len(pos) else 1.0
    return pd.DataFrame(clr(np.clip(V.values, 0, None), pseudo=pseudo), index=V.index, columns=feats)


def alpha(L):
    R = L.rel()
    with np.errstate(divide="ignore", invalid="ignore"):
        H = -(R * np.log(R.where(R > 0, 1))).sum(axis=1)
        S = (R > 0).sum(axis=1)
        J = H / np.log(S.where(S > 1, np.nan))
    return pd.DataFrame({"shannon": H, "pielou": J})


def design(C, terms, samples, L):
    """(sample ids, X, column names, owners, info) for the complete cases of the terms in `samples`"""
    kinds, use, info = {}, [], []
    depth = None
    if "depth" in terms and L.kind in ("counts", "auc"):
        depth = np.log10(L.M.sum(axis=1).clip(lower=1))
    for t in terms:
        if t == "depth":
            if depth is not None:
                kinds[t] = "numeric"
                use.append(t)
            continue
        if t == C.group_col:
            kinds[t] = "categorical"
        elif t in C.meta_cols:
            kinds[t] = C.kinds.get(t, "categorical")
            if kinds[t] in ("empty", "text"):
                info.append("%s: no usable values" % t)
                continue
        elif t in C.metric_names:
            kinds[t] = "numeric"
        else:
            info.append("%s: not in the metadata / metrics" % t)
            continue
        use.append(t)
    rows, keys = [], []
    for s in samples:
        r = {}
        for t in use:
            r[t] = repr(float(depth[s])) if t == "depth" else C.meta(s, t)
        if all(v not in ("NA", "") for v in r.values()):
            rows.append(r)
            keys.append(s)
    if len(keys) < MIN_N or not use:
        return None
    ref = {C.group_col: C.ref} if C.group_col else {}
    X, names, owner, inf2 = design_matrix(rows, use, kinds, ref)
    if len(keys) - X.shape[1] < 5:
        return None
    return keys, X, names, owner, info + inf2


def run(C, a):
    """the paper analysis of every count, protein and metabolomics layer; returns result tables"""
    terms = terms_of(a.paper_model)
    drop = qc_drop(C.release, a.paper_qc)
    an, pm, notes = [], [], []
    for key, L in C.layers.items():
        if L.kind not in ("counts", "auc", "metab"):
            continue
        t0 = time.time()
        ql = qc_layer(L)
        bad = drop.get(ql, set()) if ql else set()
        samples = [s for s in L.M.index if s not in bad and L.M.loc[s].sum() > 0]
        d = design(C, terms, samples, L)
        if d is None:
            notes.append("paper model: layer %s skipped (fewer than %d complete samples)" % (key, MIN_N))
            continue
        keys, X, names, owner, info = d
        for i in info:
            notes.append("paper model, %s: %s" % (key, i))
        # alpha diversity
        if L.kind in ("counts", "auc"):
            A = alpha(L).loc[keys]
            for col in ("shannon", "pielou"):
                ok = A[col].notna().values
                if ok.sum() >= MIN_N:
                    tt, res, sst = anova(X[ok], owner, A[col].values[ok][:, None])
                    an += anova_rows(key, col, [col], int(ok.sum()), tt, res, sst)
        # features
        Rf = response(L, a.paper_transform, a.prevalence)
        if Rf.shape[1]:
            Y = Rf.loc[keys].values.astype(float)
            sd = Y.std(axis=0)
            feats = [f for f, s_ in zip(Rf.columns, sd) if s_ > 0]
            Y = Y[:, sd > 0]
            if Y.shape[1]:
                tt, res, sst = anova(X, owner, Y)
                an += anova_rows(key, "feature", feats, len(keys), tt, res, sst)
        # beta diversity
        if L.kind == "counts":
            R = L.rel().loc[keys].values
            D, dist = bray_curtis(R), "Bray-Curtis"
        else:
            V = response(L, "clr", 0.0).loc[keys].values
            sq = (V ** 2).sum(axis=1)
            D = np.sqrt(np.maximum(sq[:, None] + sq[None, :] - 2 * V @ V.T, 0))
            np.fill_diagonal(D, 0)
            dist = "Euclidean (CLR)"
        out, (df_res, ss_res), sst = permanova_margin(D, X, owner, a.permutations)
        for t, df, ss, r2, F, p in out:
            pm.append([key, t, len(keys), dist, df, ss, r2, F, p])
        if out:
            pm.append([key, "(residual)", len(keys), dist, df_res, ss_res, ss_res / sst, np.nan, np.nan])
        C.note("paper model: %s, %d samples, %d features (%.0f s)" % (key, len(keys), Rf.shape[1], time.time() - t0))
    for n_ in notes:
        C.note(n_)
    A = pd.DataFrame(an, columns=["layer", "analysis", "feature", "term", "n", "df", "ss", "r2", "F", "p"])
    A["q"] = np.nan
    for _, idx in A[A["term"] != "(residual)"].groupby(["layer", "analysis", "term"]).groups.items():
        A.loc[idx, "q"] = bh(A.loc[idx, "p"].values)
    P = pd.DataFrame(pm, columns=["layer", "term", "n", "distance", "df", "ss", "r2", "F", "p"])
    P["q"] = np.nan
    i = P["term"] != "(residual)"
    P.loc[i, "q"] = bh(P.loc[i, "p"].values)
    return {"paper_anova": A, "paper_permanova": P, "paper_variance": variance(A, P, a.sig_q, terms)}


def variance(A, P, sig_q, order):
    """per layer and term: r2 of Shannon, of the PERMANOVA, median r2 over features, share of features with q <
    sig_q (the residual: what the model leaves unexplained)"""
    rows = []
    for layer in list(dict.fromkeys(list(A["layer"]) + list(P["layer"]))):
        a = A[A["layer"] == layer]
        p = P[P["layer"] == layer]
        terms = [t for t in order if t in set(a["term"]) | set(p["term"])] + \
            sorted((set(a["term"]) | set(p["term"])) - set(order) - {"(residual)"}) + ["(residual)"]
        for t in terms:
            sh = a[(a["analysis"] == "shannon") & (a["term"] == t)]
            fe = a[(a["analysis"] == "feature") & (a["term"] == t)]
            be = p[p["term"] == t]
            rows.append({"layer": layer, "term": t,
                         "r2_shannon": sh["r2"].iloc[0] if len(sh) else np.nan,
                         "p_shannon": sh["p"].iloc[0] if len(sh) else np.nan,
                         "r2_beta": be["r2"].iloc[0] if len(be) else np.nan,
                         "p_beta": be["p"].iloc[0] if len(be) else np.nan,
                         "r2_features_median": fe["r2"].median() if len(fe) else np.nan,
                         "features": int(len(fe)),
                         "features_q": float((fe["q"] < sig_q).mean()) if len(fe) and t != "(residual)" else np.nan})
    return pd.DataFrame(rows)
