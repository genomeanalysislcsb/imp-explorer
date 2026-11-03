"""Statistics of the cohort analysis, pure numpy (scipy only speeds up t-test p values when installed).

Taken unchanged from the sRNA-IMP cohort analysis (scripts/srna_cohort.py), so that IMP3 and sRNA cohort results
are computed the same way: Kruskal-Wallis, Mann-Whitney (exact for n <= 20), van Elteren (stratified Wilcoxon),
Spearman, BH-FDR, Bray-Curtis, PCoA, NMDS, PERMANOVA (with strata), CLR, design matrices, OLS with HC3 or
cluster-robust standard errors. Only change: PERMANOVA sums of squares as matrix products (same values, faster).
"""
import collections
import gzip
import math
import os
import re
import sys

import numpy as np

MISSING = {"", "na", "n/a", "nan", "none", "null", "-", "?"}
TERM_NAMES = {}       # term -> name, from --term-names


def load_term_names(paths):
    """term -> name from comma-separated files; formats detected per file:
       - KEGG flat files (ko, pathway, module: ENTRY / NAME / DEFINITION blocks; KO: DEFINITION, else NAME)
       - TSV term<TAB>name (e.g. KEGG REST list/ko, list/pathway, list/module; prefixes ko:, path:, md: removed)
       - Pfam-A.clans.tsv(.gz): accession, clan, clan name, short name, description -> description"""
    n0 = len(TERM_NAMES)
    for p in [x.strip() for x in (paths or "").split(",") if x.strip()]:
        if not os.path.exists(p):
            sys.stderr.write("[cohort] term names %s not found\n" % p)
            continue
        op = gzip.open if p.endswith(".gz") else open
        n1 = len(TERM_NAMES)
        with op(p, "rt", errors="replace") as fh:
            first = fh.readline()
            fh.seek(0)
            if first.startswith("ENTRY"):
                entry, name, defin = None, "", ""
                for ln in fh:
                    key = ln[:12].strip()
                    val = ln[12:].strip()
                    if key == "ENTRY":
                        entry, name, defin = val.split()[0], "", ""
                    elif key == "NAME" and entry and not name:
                        name = val
                    elif key == "DEFINITION" and entry and not defin:
                        defin = val
                    elif ln.startswith("///") and entry:
                        label = defin if entry.startswith("K") and defin else name or defin
                        if label:
                            TERM_NAMES.setdefault(entry, re.sub(r"\s+", " ", label))
                        entry = None
            else:
                for ln in fh:
                    if ln.startswith("#"):
                        continue
                    f = ln.rstrip("\n").split("\t")
                    if len(f) < 2:                    # not tab-separated: id, whitespace, name
                        f = ln.strip().split(None, 1)
                    if len(f) < 2 or not f[0]:
                        continue
                    term = re.sub(r"^(ko|path|md|map):", "", f[0].strip())
                    name = f[4] if len(f) >= 5 and term.startswith("PF") else f[1]
                    TERM_NAMES.setdefault(term, name.strip())
        sys.stderr.write("[cohort] %d term names from %s\n" % (len(TERM_NAMES) - n1, p))
    return len(TERM_NAMES) - n0



def rankdata(x):
    x = np.asarray(x, float)
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(len(x))
    xs = x[order]
    i = 0
    while i < len(x):
        j = i
        while j + 1 < len(x) and xs[j + 1] == xs[i]:
            j += 1
        ranks[order[i:j + 1]] = (i + j) / 2.0 + 1
        i = j + 1
    return ranks


def _tie_term(r):
    _, c = np.unique(r, return_counts=True)
    return float(np.sum(c ** 3 - c))


def _gammq(a, x):
    """regularised upper incomplete gamma Q(a, x) (Numerical Recipes)"""
    if x <= 0:
        return 1.0
    gln = math.lgamma(a)
    if x < a + 1:
        ap, s, d = a, 1.0 / a, 1.0 / a
        for _ in range(500):
            ap += 1
            d *= x / ap
            s += d
            if abs(d) < abs(s) * 1e-12:
                break
        return max(0.0, 1.0 - s * math.exp(-x + a * math.log(x) - gln))
    b = x + 1 - a
    c, d = 1e300, 1.0 / b
    h = d
    for i in range(1, 500):
        an = -i * (i - a)
        b += 2
        d = an * d + b
        d = 1e-300 if abs(d) < 1e-300 else d
        c = b + an / c
        c = 1e-300 if abs(c) < 1e-300 else c
        d = 1.0 / d
        de = d * c
        h *= de
        if abs(de - 1) < 1e-12:
            break
    return min(1.0, math.exp(-x + a * math.log(x) - gln) * h)


def kruskal(groups):
    groups = [np.asarray(g, float) for g in groups if len(g) > 0]
    if len(groups) < 2:
        return float("nan"), float("nan")
    allv = np.concatenate(groups)
    n = len(allv)
    if n < 3 or np.all(allv == allv[0]):
        return float("nan"), float("nan")
    r = rankdata(allv)
    h, i = 0.0, 0
    for g in groups:
        rg = r[i:i + len(g)]
        h += len(g) * (rg.mean() - (n + 1) / 2.0) ** 2
        i += len(g)
    h *= 12.0 / (n * (n + 1))
    c = 1 - _tie_term(r) / (n ** 3 - n)
    h = h / c if c > 0 else float("nan")
    return h, _gammq((len(groups) - 1) / 2.0, h / 2.0)


def _mw_exact_p(r, n1, stat):
    """exact two-sided p of the rank sum of n1 of the (mid)ranks r (DP over doubled integer ranks)"""
    rr = [int(round(2 * x)) for x in r]
    tot = sum(rr)
    # dp[k][s] = number of subsets of size k with doubled rank sum s
    dp = [collections.Counter() for _ in range(n1 + 1)]
    dp[0][0] = 1
    for x in rr:
        for k in range(min(n1, len(rr)) - 1, -1, -1):
            for sm, c in dp[k].items():
                dp[k + 1][sm + x] += c
    dist = dp[n1]
    allc = sum(dist.values())
    mean = n1 * tot / len(rr)
    obs = abs(2 * stat - mean)
    return min(1.0, sum(c for sm, c in dist.items() if abs(sm - mean) >= obs - 1e-9) / allc)


def mannwhitney(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    n1, n2 = len(a), len(b)
    if n1 < 2 or n2 < 2:
        return float("nan"), float("nan")
    r = rankdata(np.concatenate([a, b]))
    u = r[:n1].sum() - n1 * (n1 + 1) / 2.0
    n = n1 + n2
    if n <= 20:                                   # exact p for small samples (normal approx. too liberal)
        if np.all(r == r[0]):
            return float("nan"), float("nan")
        return u / (n1 * n2), _mw_exact_p(r, n1, r[:n1].sum())
    var = n1 * n2 / 12.0 * ((n + 1) - _tie_term(r) / (n * (n - 1)))
    if var <= 0:
        return float("nan"), float("nan")
    z = (u - n1 * n2 / 2.0) / math.sqrt(var)
    return u / (n1 * n2), math.erfc(abs(z) / math.sqrt(2))   # effect: P(a > b) (AUC)


def van_elteren(a, ba, b, bb):
    """stratified Wilcoxon (van Elteren) test of a vs b within blocks: values a with block labels ba, values b with
    block labels bb; only blocks with both a and b count. Returns (block-weighted P(a > b), two-sided p, blocks).
    With one informative block it is the Mann-Whitney test (exact for small n)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    ba, bb = np.asarray(ba, object), np.asarray(bb, object)
    blocks = sorted(set(ba) & set(bb))
    if not blocks:
        return float("nan"), float("nan"), 0
    if len(blocks) == 1:
        auc, p = mannwhitney(a[ba == blocks[0]], b[bb == blocks[0]])
        return auc, p, 1
    num = var = 0.0
    aucs, ws = [], []
    for bl in blocks:
        x, y = a[ba == bl], b[bb == bl]
        n1, n2 = len(x), len(y)
        n = n1 + n2
        r = rankdata(np.concatenate([x, y]))
        w = 1.0 / (n + 1)
        v = n1 * n2 / 12.0 * ((n + 1) - _tie_term(r) / (n * (n - 1)))
        num += w * (r[:n1].sum() - n1 * (n + 1) / 2.0)
        var += w * w * v
        aucs.append((r[:n1].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n2))
        ws.append(n1 * n2 / (n + 1.0))
    if var <= 0:
        return float("nan"), float("nan"), len(blocks)
    z = num / math.sqrt(var)
    return float(np.average(aucs, weights=ws)), math.erfc(abs(z) / math.sqrt(2)), len(blocks)


def wilcoxon_paired(a, b):
    d = np.asarray(a, float) - np.asarray(b, float)
    d = d[d != 0]
    n = len(d)
    if n < 5:
        return float("nan"), float("nan"), n
    r = rankdata(np.abs(d))
    w = r[d > 0].sum()
    var = n * (n + 1) * (2 * n + 1) / 24.0 - _tie_term(r) / 48.0
    z = (w - n * (n + 1) / 4.0) / math.sqrt(var)
    return float(np.median(d)), math.erfc(abs(z) / math.sqrt(2)), n


def spearman(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = ~(np.isnan(x) | np.isnan(y))
    x, y = x[ok], y[ok]
    n = len(x)
    if n < 4 or np.all(x == x[0]) or np.all(y == y[0]):
        return float("nan"), float("nan"), n
    rho = float(np.corrcoef(rankdata(x), rankdata(y))[0, 1])
    if abs(rho) >= 1:
        return rho, 0.0, n
    t = rho * math.sqrt((n - 2) / (1 - rho ** 2))         # t approximation, normal tail for df >= 30
    df = n - 2
    p = _t_sf2(abs(t), df)
    return rho, p, n


def _betacf(a, b, x):
    qab, qap, qam = a + b, a + 1, a - 1
    c, d = 1.0, 1 - qab * x / qap
    d = 1.0 / (d if abs(d) > 1e-300 else 1e-300)
    h = d
    for m in range(1, 300):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1 + aa * d
        d = 1.0 / (d if abs(d) > 1e-300 else 1e-300)
        c = 1 + aa / c
        c = c if abs(c) > 1e-300 else 1e-300
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1 + aa * d
        d = 1.0 / (d if abs(d) > 1e-300 else 1e-300)
        c = 1 + aa / c
        c = c if abs(c) > 1e-300 else 1e-300
        de = d * c
        h *= de
        if abs(de - 1) < 1e-12:
            break
    return h


def _betai(a, b, x):
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    bt = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log(1 - x))
    if x < (a + 1) / (a + b + 2):
        return bt * _betacf(a, b, x) / a
    return 1 - bt * _betacf(b, a, 1 - x) / b


def _t_sf2(t, df):
    """two-sided p of Student's t"""
    return _betai(df / 2.0, 0.5, df / (df + t * t))


try:
    from scipy.special import stdtr as _stdtr

    def t_p(T, df):
        """two-sided p of an array of t values"""
        T = np.asarray(T, float)
        return np.where(np.isfinite(T), 2 * _stdtr(df, -np.abs(np.nan_to_num(T))), np.nan)
except ImportError:                       # pragma: no cover
    def t_p(T, df):
        return np.vectorize(lambda t: _t_sf2(abs(t), df) if np.isfinite(t) else np.nan)(T)


def bh(pvals):
    p = np.asarray(pvals, float)
    q = np.full(len(p), np.nan)
    ok = ~np.isnan(p)
    if ok.sum() == 0:
        return q
    pv = p[ok]
    order = np.argsort(pv)
    m = len(pv)
    adj = pv[order] * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(adj, 1)
    q[ok] = out
    return q


def bray_curtis(X):
    """X: samples x features (relative abundances) -> distance matrix"""
    n = X.shape[0]
    D = np.zeros((n, n))
    for i in range(n):
        num = np.abs(X[i] - X).sum(axis=1)
        den = (X[i] + X).sum(axis=1)
        D[i] = np.where(den > 0, num / np.where(den > 0, den, 1), 0)
    return (D + D.T) / 2


def log_values(X):
    """metabolomics values -> log2(x + half the smallest positive value of the feature); negative input (already on
    a log scale) is returned unchanged"""
    X = np.asarray(X, float)
    if (X < 0).any():
        return X
    out = np.empty_like(X)
    for j in range(X.shape[1]):
        pos = X[:, j][X[:, j] > 0]
        out[:, j] = np.log2(X[:, j] + (pos.min() / 2 if len(pos) else 1.0))
    return out


def euclidean_z(X):
    """Euclidean distances of autoscaled (z-scored) log values: the usual beta diversity of metabolite profiles"""
    L = log_values(X)
    sd = L.std(axis=0)
    L = (L[:, sd > 0] - L[:, sd > 0].mean(axis=0)) / sd[sd > 0]
    G = (L ** 2).sum(axis=1)
    D = np.sqrt(np.maximum(G[:, None] + G[None, :] - 2 * L @ L.T, 0))
    np.fill_diagonal(D, 0)
    return D


def pcoa(D, k=2):
    n = D.shape[0]
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ (D ** 2) @ J
    w, v = np.linalg.eigh(B)
    idx = np.argsort(w)[::-1]
    w, v = w[idx], v[:, idx]
    pos = w[w > 0].sum() or 1
    coords = v[:, :k] * np.sqrt(np.maximum(w[:k], 0))
    return coords, w[:k] / pos


# ----------------------------------------------------------------------------- adjusted linear models
def clr(C, pseudo=0.5):
    """centred log-ratio of a count matrix (samples x features) with a pseudocount"""
    L = np.log(C + pseudo)
    return L - L.mean(axis=1, keepdims=True)


def design_matrix(rows, terms, kinds, levels_first=None):
    """rows: list of dicts; terms: variable names; kinds: name -> numeric / categorical.
    Returns (X, column names, column -> term, info) with an intercept; categorical terms as treatment
    dummies (reference = levels_first[term] or the most frequent level); numeric terms scaled to SD units.
    Linearly dependent columns (e.g. a covariate fully confounded with the group) are dropped and reported."""
    cols, names, owner, info = [np.ones(len(rows))], ["(intercept)"], ["(intercept)"], []
    for t in terms:
        vals = [r.get(t, "NA") for r in rows]
        if kinds.get(t) == "numeric":
            x = np.array([float(v) for v in vals])
            sd = x.std(ddof=1) or 1.0
            cols.append((x - x.mean()) / sd)
            names.append("%s (per SD %s)" % (t, fmt(sd, 3)))
            owner.append(t)
        else:
            lv = sorted(set(vals), key=lambda v: (-vals.count(v), v))
            ref = (levels_first or {}).get(t)
            if ref not in lv:                       # natural reference for yes/no or normal/abnormal variables
                ref = next((v for v in lv if str(v).strip().lower() in ("no", "none", "normal", "no problems", "0",
                                                                         "false", "never")), lv[0])
            for l in [x for x in lv if x != ref]:
                cols.append(np.array([1.0 if v == l else 0.0 for v in vals]))
                names.append("%s: %s vs %s" % (t, l, ref))
                owner.append(t)
    X = np.column_stack(cols)
    keep = [0]
    for j in range(1, X.shape[1]):
        if np.linalg.matrix_rank(X[:, keep + [j]]) > len(keep):
            keep.append(j)
        else:
            info.append("dropped (collinear): %s" % names[j])
    return X[:, keep], [names[j] for j in keep], [owner[j] for j in keep], info


MODEL_MIN_N = 20          # minimum libraries per model fit (and >= 3 x terms)
MODEL_MIN_N_ASIDX = 40    # antisense index: only libraries with >= 10 reads on the term count -> require more
MODEL_MIN_LEVEL = 3       # minimum libraries in and outside every 0/1 design column (e.g. a group or site level)
MODEL_MAX_LEVERAGE = 0.95 # fits with a library of leverage above this are skipped (HC3 SE undefined / collapses)


def design_ok(X, off=0):
    """enough libraries for a stable fit: n, libraries per categorical level, no library that alone carries a term"""
    n, p = X.shape
    if n < max(MODEL_MIN_N, 3 * (p - off)):
        return False
    for c in range(1 if off == 0 else off, p):
        col = X[:, c]
        u = np.unique(col[~np.isnan(col)])
        if len(u) <= 2:                                  # dummy (0/1) or a binary numeric term
            k = int(np.sum(col == u.max())) if len(u) else 0
            if min(k, n - k) < MODEL_MIN_LEVEL:
                return False
    try:
        Xc = X[:, off:]
        h = np.einsum("ij,ji->i", Xc, np.linalg.pinv(Xc.T @ Xc) @ Xc.T)
        if h.max() > MODEL_MAX_LEVERAGE:
            return False
    except np.linalg.LinAlgError:
        return False
    return True


def ols_hc3(X, Y):
    """OLS of every column of Y (n x m) on X (n x p): estimates, HC3 robust SE, t, two-sided p (t, df n-p)"""
    n, p = X.shape
    XtXi = np.linalg.pinv(X.T @ X)
    A = XtXi @ X.T                                   # p x n
    B = A @ Y                                        # p x m
    E = Y - X @ B
    h = np.einsum("ij,ji->i", X, A)                  # leverages
    w = 1.0 / np.maximum(1 - h, 1e-8) ** 2
    V = (A ** 2) @ (E ** 2 * w[:, None])             # p x m: HC3 variances
    SE = np.sqrt(np.maximum(V, 0))
    with np.errstate(divide="ignore", invalid="ignore"):
        T = np.where(SE > 0, B / SE, np.nan)
    df = n - p
    P = t_p(T, df) if df > 0 else T * np.nan
    return B, SE, T, P, df


def ols_cluster(X, Y, clusters):
    """OLS with cluster-robust (CR1) standard errors: residuals of the same cluster (e.g. family) may be
    correlated; t-test with G - 1 degrees of freedom"""
    n, p = X.shape
    cl = np.asarray(clusters, object)
    levels = list(set(cl))
    G = len(levels)
    XtXi = np.linalg.pinv(X.T @ X)
    B = XtXi @ X.T @ Y
    E = Y - X @ B
    meat = np.zeros((p, Y.shape[1]))
    A = XtXi @ X.T                                      # p x n
    for g in levels:
        idx = np.where(cl == g)[0]
        S = A[:, idx] @ E[idx]                          # p x m: (X'X)^-1 X_g' e_g
        meat += S ** 2
    c = G / max(G - 1, 1) * (n - 1) / max(n - p, 1)
    SE = np.sqrt(np.maximum(meat * c, 0))
    with np.errstate(divide="ignore", invalid="ignore"):
        T = np.where(SE > 0, B / SE, np.nan)
    df = G - 1
    P = t_p(T, df) if df > 0 else T * np.nan
    return B, SE, T, P, df


def demean_within(M, blocks):
    """subtract the block mean from every column (fixed effects by the within transformation)"""
    M = np.array(M, float)
    bl = np.asarray(blocks, object)
    for g in set(bl):
        idx = np.where(bl == g)[0]
        M[idx] -= M[idx].mean(axis=0)
    return M


def _pava(y):
    """isotonic (non-decreasing) regression, pool-adjacent-violators"""
    y = list(map(float, y))
    vals, wts, lens = [], [], []
    for v in y:
        vals.append(v); wts.append(1.0); lens.append(1)
        while len(vals) > 1 and vals[-2] > vals[-1]:
            w = wts[-2] + wts[-1]
            v2 = (vals[-2] * wts[-2] + vals[-1] * wts[-1]) / w
            l2 = lens[-2] + lens[-1]
            vals[-2:], wts[-2:], lens[-2:] = [v2], [w], [l2]
    return np.repeat(vals, lens)


def _stress1(D, X):
    """Kruskal stress-1 of configuration X for dissimilarities D (monotone fit, secondary ties)"""
    iu = np.triu_indices(D.shape[0], 1)
    d = D[iu]
    dv = np.sqrt(((X[:, None, :] - X[None, :, :]) ** 2).sum(axis=2))[iu]
    order = np.lexsort((dv, d))
    dh = np.empty_like(dv)
    dh[order] = _pava(dv[order])
    return float(np.sqrt(((dv - dh) ** 2).sum() / max((dv ** 2).sum(), 1e-12)))


def nmds(D, k=2, starts=10, max_iter=300, eps=1e-7, seed=1):
    """non-metric MDS (Kruskal) by SMACOF with monotone regression on the dissimilarities D (n x n).
    Ties in D: secondary approach (tied dissimilarities keep the order of the fitted distances).
    Returns (coordinates n x k, Kruskal stress-1); best of a PCoA start and `starts` random starts."""
    n = D.shape[0]
    if n < k + 2:
        return None, float("nan")
    iu = np.triu_indices(n, 1)
    d = D[iu]
    rng = np.random.default_rng(seed)
    inits = [pcoa(D, k)[0]] + [rng.normal(size=(n, k)) for _ in range(starts)]
    best = (None, float("inf"))
    for X in inits:
        X = X - X.mean(axis=0)
        sx = np.sqrt((X ** 2).sum() / n) or 1.0
        X = X / sx
        prev = float("inf")
        for _ in range(max_iter):
            diff = X[:, None, :] - X[None, :, :]
            dist = np.sqrt((diff ** 2).sum(axis=2))
            dv = dist[iu]
            order = np.lexsort((dv, d))                 # by dissimilarity, ties by current distance
            dh = np.empty_like(dv)
            dh[order] = _pava(dv[order])
            dh *= np.sqrt(len(dh) / max((dh ** 2).sum(), 1e-12))   # normalise disparities
            stress = np.sqrt(((dv - dh) ** 2).sum() / max((dv ** 2).sum(), 1e-12))
            if prev - stress < eps:
                break
            prev = stress
            Dh = np.zeros((n, n))
            Dh[iu] = dh
            Dh = Dh + Dh.T
            with np.errstate(divide="ignore", invalid="ignore"):
                Bm = -np.where(dist > 0, Dh / dist, 0.0)
            Bm[np.diag_indices(n)] = 0
            Bm[np.diag_indices(n)] = -Bm.sum(axis=1)
            X = Bm @ X / n                               # Guttman transform
        stress = _stress1(D, X)                          # stress of the configuration actually returned
        if stress < best[1]:
            best = (X.copy(), stress)
    X, st = best
    # orient: first axis = direction of largest variance
    w, v = np.linalg.eigh(np.cov(X.T))
    X = X @ v[:, np.argsort(w)[::-1]]
    return X, float(st)


def isolated_samples(D, max_frac=0.05):
    """indices of samples that are further from every other sample than any two remaining samples are from each
    other (min distance to the rest > max distance within the rest). Such a sample makes NMDS degenerate: the rest
    collapse to one point and stress goes to 0. Removed one at a time, at most max_frac of the samples."""
    keep = list(range(D.shape[0]))
    out = []
    while len(keep) > 6 and len(out) < max(1, int(max_frac * D.shape[0])):
        sub = D[np.ix_(keep, keep)]
        off = sub + np.diag([np.inf] * len(keep))
        nearest = off.min(axis=1)
        i = int(np.argmax(nearest))
        rest = [j for j in range(len(keep)) if j != i]
        if nearest[i] <= sub[np.ix_(rest, rest)].max():
            break
        out.append(keep.pop(i))
    return out


def _ssw_fn(D2):
    """within-group sum of squares of a labelling: sum over groups of (sum of squared distances within the group) /
    group size, as one matrix product (same value as the pairwise loop of srna_cohort.py, much faster)"""
    def ssw(lab):
        u, inv = np.unique(np.asarray(lab).astype(str), return_inverse=True)
        G = np.zeros((len(inv), len(u)))
        G[np.arange(len(inv)), inv] = 1.0
        n = G.sum(axis=0)
        within = np.einsum("ij,ij->j", G, D2 @ G) / 2.0
        return float(np.sum(np.where(n > 1, within / np.maximum(n, 1), 0.0)))
    return ssw


def permanova(D, labels, perms=999, seed=1, strata=None):
    labels = np.asarray(labels)
    lv = sorted(set(labels))
    n, a = len(labels), len(lv)
    if a < 2 or n <= a or min(collections.Counter(labels).values()) < 2:
        return float("nan"), float("nan"), float("nan")
    D2 = D ** 2
    sst = D2[np.triu_indices(n, 1)].sum() / n
    ssw = _ssw_fn(D2)

    w = ssw(labels)
    f = ((sst - w) / (a - 1)) / (w / (n - a)) if w > 0 else float("nan")
    rng = np.random.default_rng(seed)
    groups_ = None
    if strata is not None:                       # permutations only within strata (blocks)
        strata = np.asarray(strata, object)
        groups_ = [np.where(strata == s_)[0] for s_ in set(strata)]

    def shuffle():
        if groups_ is None:
            return rng.permutation(labels)
        out = labels.copy()
        for idx in groups_:
            out[idx] = labels[rng.permutation(idx)]
        return out
    ge = 0
    for _ in range(perms):
        wp = ssw(shuffle())
        fp = ((sst - wp) / (a - 1)) / (wp / (n - a)) if wp > 0 else 0
        ge += fp >= f
    return f, (sst - w) / sst if sst else float("nan"), (ge + 1) / (perms + 1)


def permanova_partial(D, labels, strata, perms=999, seed=1):
    """PERMANOVA of `labels` after `strata` (e.g. block after group): partial R² = SS(labels | strata) / SS total;
    labels permuted within strata. Returns (F, partial R², p)."""
    labels, strata = np.asarray(labels, object), np.asarray(strata, object)
    n = len(labels)
    D2 = D ** 2
    sst = D2[np.triu_indices(n, 1)].sum() / n
    ssw = _ssw_fn(D2)
    cell = lambda lab: np.array(["%s|%s" % (a_, b_) for a_, b_ in zip(strata, lab)], object)
    w_s = ssw(strata)
    a_s, a_c = len(set(strata)), len(set(cell(labels)))
    df1, df2 = a_c - a_s, n - a_c
    if df1 < 1 or df2 < 1 or not sst:
        return float("nan"), float("nan"), float("nan")
    w_c = ssw(cell(labels))
    stat = lambda wc: ((w_s - wc) / df1) / (wc / df2) if wc > 0 else float("inf")
    f = stat(w_c)
    rng = np.random.default_rng(seed)
    groups_ = [np.where(strata == s_)[0] for s_ in set(strata)]
    ge = 0
    for _ in range(perms):
        out = labels.copy()
        for idx in groups_:
            out[idx] = labels[rng.permutation(idx)]
        ge += stat(ssw(cell(out))) >= f
    return f, (w_s - w_c) / sst, (ge + 1) / (perms + 1)


# ----------------------------------------------------------------------------- metadata
def is_missing(v):
    return v is None or str(v).strip().lower() in MISSING



def column_kind(values):
    vals = [v for v in values if not is_missing(v)]
    if not vals:
        return "empty"
    try:
        [float(v) for v in vals]
        return "numeric" if len(set(vals)) > 2 else "categorical"
    except ValueError:
        return "categorical"



def fmt(v, digits=4):
    if isinstance(v, float):
        if math.isnan(v):
            return "NA"
        return ("%.*g" % (digits, v))
    return str(v)



def chi2_test(tab):
    tab = tab[tab.sum(axis=1) > 0][:, tab.sum(axis=0) > 0]
    if tab.shape[0] < 2 or tab.shape[1] < 2:
        return float("nan")
    e = tab.sum(axis=1, keepdims=True) * tab.sum(axis=0, keepdims=True) / tab.sum()
    x2 = float(((tab - e) ** 2 / np.where(e > 0, e, 1)).sum())
    return _gammq((tab.shape[0] - 1) * (tab.shape[1] - 1) / 2.0, x2 / 2.0)


