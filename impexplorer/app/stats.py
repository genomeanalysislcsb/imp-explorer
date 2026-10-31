"""Small statistics helpers (numpy / scipy): Bray-Curtis PCoA, group tests with BH-FDR."""
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


def permanova(D, labels, perms=499, seed=1):
    labels = np.asarray(labels)
    groups = np.unique(labels)
    n = len(labels)
    if len(groups) < 2 or n <= len(groups):
        return np.nan, np.nan, np.nan
    D2 = D ** 2

    def f_stat(lab):
        sst = D2[np.triu_indices(n, 1)].sum() / n
        ssw = 0.0
        for g in groups:
            i = np.where(lab == g)[0]
            if len(i) > 1:
                ssw += D2[np.ix_(i, i)][np.triu_indices(len(i), 1)].sum() / len(i)
        ssb = sst - ssw
        return (ssb / (len(groups) - 1)) / (ssw / (n - len(groups))), ssb / sst

    F, r2 = f_stat(labels)
    rng = np.random.default_rng(seed)
    ge = sum(f_stat(rng.permutation(labels))[0] >= F for _ in range(perms))
    return F, r2, (ge + 1) / (perms + 1)


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
