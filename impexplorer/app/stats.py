"""Small statistics helpers (numpy / scipy): Bray-Curtis or Jaccard PCoA, PERMANOVA (marginal, as adonis2),
group tests with BH-FDR."""
import numpy as np
import pandas as pd
from scipy import stats


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    if n == 0:
        return p
    o = np.argsort(p)
    q = p[o] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1)
    return out


def bray_curtis(X):
    X = np.asarray(X, float)
    s = X.sum(axis=1)
    num = np.abs(X[:, None, :] - X[None, :, :]).sum(axis=2)
    den = s[:, None] + s[None, :]
    return np.where(den > 0, num / np.where(den > 0, den, 1), 0.0)


def pcoa(D, k=2):
    n = D.shape[0]
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ (D ** 2) @ J
    w, v = np.linalg.eigh(B)
    idx = np.argsort(w)[::-1]
    w, v = w[idx], v[:, idx]
    pos = w > 1e-12
    expl = w[pos] / w[pos].sum() if pos.any() else np.zeros(k)
    coords = v[:, :k] * np.sqrt(np.clip(w[:k], 0, None))
    return coords, expl[:k]


def jaccard(X, binary=True):
    """Jaccard distance: presence / absence (default) or quantitative (2 BC / (1 + BC)); see cohort.stats.jaccard"""
    from ..cohort.stats import jaccard as _j
    return _j(X, binary)


DISTANCES = {"bray": "Bray-Curtis", "jaccard": "Jaccard (presence / absence)"}


def distance(X, method="bray"):
    """samples x features (relative abundances) -> distance matrix: bray (Bray-Curtis) or jaccard"""
    return bray_curtis(X) if method == "bray" else jaccard(X, binary=True)


def permanova_terms(D, data, perms=499, seed=1, strata=None):
    """multivariable PERMANOVA as vegan::adonis2(by = "margin"): one row per term (marginal SS, R2, F, p from
    permutations under the reduced model), the residual and the total. data: DataFrame, one column per term
    (numeric columns as covariates, others as factors), rows in the order of D."""
    from ..cohort.stats import design_matrix, permanova_margin
    data = pd.DataFrame(data).reset_index(drop=True)
    terms = list(data.columns)
    kinds = {t: "numeric" if pd.api.types.is_numeric_dtype(data[t]) and data[t].nunique() > 2 else "categorical"
             for t in terms}
    rows = [{t: (v if kinds[t] == "numeric" else str(v)) for t, v in r.items()} for r in data.to_dict("records")]
    X, _, owner, _ = design_matrix(rows, terms, kinds)
    out, (df_res, ss_res), sst = permanova_margin(np.asarray(D, float), X, owner, perms, seed, strata)
    cols = ["term", "df", "SS", "R2", "F", "p"]
    if not out:
        return pd.DataFrame(columns=cols)
    res = [list(r) for r in out]
    res.append(["(residual)", df_res, ss_res, ss_res / sst, np.nan, np.nan])
    res.append(["(total)", D.shape[0] - 1, sst, 1.0, np.nan, np.nan])
    return pd.DataFrame(res, columns=cols)


def permanova(D, labels, perms=499, seed=1, covariates=None):
    """PERMANOVA of one grouping (labels), optionally adjusted for covariates (DataFrame, rows as D): the group's
    marginal R2, F and p (reduced-model permutations, as adonis2 by = "margin"). Returns (F, R2, p)."""
    labels = np.asarray(labels)
    n = len(labels)
    if len(np.unique(labels)) < 2 or n <= len(np.unique(labels)):
        return np.nan, np.nan, np.nan
    data = pd.DataFrame({"__group__": labels.astype(str)})
    if covariates is not None:
        data = pd.concat([data, pd.DataFrame(covariates).reset_index(drop=True)], axis=1)
    t = permanova_terms(D, data, perms, seed)
    r = t[t["term"] == "__group__"]
    if r.empty:
        return np.nan, np.nan, np.nan
    return float(r["F"].iloc[0]), float(r["R2"].iloc[0]), float(r["p"].iloc[0])

def group_tests(M, groups, min_prev=0.1):
    """M: samples x features (relative abundances); groups: labels per sample.
    Kruskal-Wallis (2 groups: Mann-Whitney) per feature present in >= min_prev of samples; BH-FDR."""
    groups = pd.Series(np.asarray(groups), index=M.index)
    lv = [g for g in groups.unique() if (groups == g).sum() >= 2]
    rows = []
    if len(lv) < 2:
        return pd.DataFrame(columns=["feature", "p", "q"])
    keep = (M > 0).mean() >= min_prev
    for f in M.columns[keep.values]:
        xs = [M.loc[groups == g, f].values for g in lv]
        try:
            p = stats.mannwhitneyu(*xs).pvalue if len(xs) == 2 else stats.kruskal(*xs).pvalue
        except ValueError:
            p = 1.0
        r = {"feature": f, "p": p}
        for g, x in zip(lv, xs):
            r["median % " + str(g)] = 100 * float(np.median(x))
        rows.append(r)
    df = pd.DataFrame(rows)
    if len(df):
        df["q"] = bh(df["p"].values)
        df = df.sort_values("p")[["feature", "q", "p"] + [c for c in df.columns if c.startswith("median")]]
    return df
