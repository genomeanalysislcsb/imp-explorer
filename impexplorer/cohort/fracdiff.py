"""Intracellular vs extracellular per feature (--fraction-pairs): which taxa and functions are enriched in one
fraction of the same samples, as in a paper's figure comparing the cellular and the extracellular fraction.

  --fraction-pairs auto
  --fraction-pairs "kraken.genus.mg,kraken.genus.mg_ext; prot.genus.intra,prot.genus.extra; srna.genus.srna_int,srna.genus.srna_ext"

Per pair of layers with the same features (a = first, b = second), on the samples with both (samples flagged FAIL in a
layer's sample_qc left out, as --paper-qc):

  features   detected in >= --prevalence of the shared samples in either layer; both profiles closed to these features
             and CLR-transformed (counts: pseudo-count 1; protein AUC: half the smallest positive value), so a value is
             the feature's log-ratio to the geometric mean of the shared features in that fraction
  test       paired Wilcoxon signed-rank test of CLR(b) - CLR(a) over the samples (normal approximation with tie and
             continuity correction; zero differences dropped, fewer than 6 non-zero ones: p = 1), BH-FDR per pair
  effect     matched-pairs rank-biserial correlation r = (W+ - W-) / (W+ + W-), -1..1: > 0 higher in b (e.g.
             extracellular), the paired counterpart of Cliff's delta; median log2 ratio b / a of the closed profiles
  enriched   q < --sig-q and |r| >= --fraction-min-r (default 0.5: with ~180 pairs almost any shift is significant)
"""
import numpy as np
import pandas as pd
from scipy.stats import norm

from .concord import matrices
from .stats import bh, clr

AUTO = [("kraken.genus.mg", "kraken.genus.mg_ext"), ("contig.genus.mg", "contig.genus.mg_ext"),
        ("prot.genus.intra", "prot.genus.extra"),
        ("srna.genus.srna_int", "srna.genus.srna_ext"),
        ("kegg_ko.mg", "kegg_ko.mg_ext"), ("kegg_module.mg", "kegg_module.mg_ext"), ("cazy.mg", "cazy.mg_ext"),
        ("pfam.mg", "pfam.mg_ext"), ("prot.kegg_ko.intra", "prot.kegg_ko.extra"),
        ("prot.pfam.intra", "prot.pfam.extra")]
MIN_PAIRS = 10


def parse(spec):
    spec = (spec or "").strip()
    if not spec:
        return []
    if spec == "auto":
        return list(AUTO)
    out = []
    for part in spec.split(";"):
        x = [k.strip() for k in part.split(",") if k.strip()]
        if len(x) == 2:
            out.append(tuple(x))
    return out


def needed_layers(pairs):
    return [k for p in pairs for k in p if not k.startswith("srna.")]


def wilcoxon_paired(D):
    """D: samples x features of paired differences -> W+, W-, n (non-zero pairs), z, p (two-sided) per feature"""
    k = D.shape[1]
    wp, wm, nn, z = (np.full(k, np.nan) for _ in range(4))
    for j in range(k):
        d = D[:, j]
        d = d[np.isfinite(d) & (np.abs(d) > 1e-12)]
        n = len(d)
        nn[j] = n
        if n < 6:                      # (almost) no paired differences: not testable, p = 1, r = 0
            wp[j] = wm[j] = 0.0
            z[j] = 0.0
            continue
        a = np.abs(d)
        order = np.argsort(a, kind="mergesort")
        r = np.empty(n)
        r[order] = np.arange(1, n + 1)
        u, inv, cnt = np.unique(a, return_inverse=True, return_counts=True)
        if (cnt > 1).any():                                  # average ranks of ties
            s = np.bincount(inv, weights=r)
            r = (s / cnt)[inv]
        wp[j] = r[d > 0].sum()
        wm[j] = r[d < 0].sum()
        mu = n * (n + 1) / 4
        var = n * (n + 1) * (2 * n + 1) / 24 - ((cnt ** 3 - cnt).sum()) / 48
        if var > 0:
            diff = wp[j] - mu
            z[j] = (diff - 0.5 * np.sign(diff)) / np.sqrt(var)
    p = 2 * norm.sf(np.abs(z))
    return wp, wm, nn, z, p


def compare(A, B, prevalence):
    """paired comparison of two layers (samples x features, raw counts / AUC): one row per feature"""
    ss = sorted(set(A.index) & set(B.index))
    f = sorted(set(A.columns) & set(B.columns))
    cols = ["feature", "n", "prev_a", "prev_b", "median_rel_a", "median_rel_b", "median_log2_ratio", "rbc", "z", "p"]
    if len(ss) < MIN_PAIRS or len(f) < 2:
        return pd.DataFrame(columns=cols)
    X, Y = A.loc[ss, f].astype(float), B.loc[ss, f].astype(float)
    keep = [c for c in f if (X[c] > 0).mean() >= prevalence or (Y[c] > 0).mean() >= prevalence]
    if len(keep) < 2:
        return pd.DataFrame(columns=cols)
    X, Y = X[keep], Y[keep]
    ok = (X.sum(axis=1) > 0) & (Y.sum(axis=1) > 0)
    X, Y = X[ok], Y[ok]
    if len(X) < MIN_PAIRS:
        return pd.DataFrame(columns=cols)

    def closed_clr(M):
        v = M.values
        pos = v[v > 0]
        pseudo = 1.0 if v.max() > 1 and np.allclose(v, np.round(v)) else (pos.min() / 2 if len(pos) else 1e-9)
        return clr(v, pseudo=pseudo)

    D = closed_clr(Y) - closed_clr(X)
    wp, wm, nn, z, p = wilcoxon_paired(D)
    RX, RY = X.div(X.sum(axis=1), axis=0), Y.div(Y.sum(axis=1), axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        rbc = np.where(wp + wm > 0, (wp - wm) / (wp + wm), 0.0)
    return pd.DataFrame({"feature": keep, "n": len(X), "prev_a": (X > 0).mean().values, "prev_b": (Y > 0).mean().values,
                         "median_rel_a": RX.median().values, "median_rel_b": RY.median().values,
                         "median_log2_ratio": np.median(D, axis=0) / np.log(2), "rbc": rbc, "z": z, "p": p})


def run(C, a):
    pairs = parse(a.fraction_pairs)
    Ms = matrices(C, a, sorted({k for pr in pairs for k in pr}))
    rows, summ = [], []
    q_sig = getattr(a, "sig_q", 0.1)
    min_r = float(getattr(a, "fraction_min_r", 0.5) or 0)
    for la, lb in pairs:
        if la not in Ms or lb not in Ms:
            continue
        d = compare(Ms[la], Ms[lb], a.prevalence)
        if d.empty:
            C.note("fraction pair %s / %s: fewer than %d shared samples or no shared features" % (la, lb, MIN_PAIRS))
            continue
        d["q"] = bh(d["p"].values)
        d.insert(0, "layer_b", lb)
        d.insert(0, "layer_a", la)
        rows.append(d)
        sig = d["q"] < q_sig
        strong = sig & (d["rbc"].abs() >= min_r)
        summ.append({"layer_a": la, "layer_b": lb, "samples": int(d["n"].iloc[0]), "features": len(d),
                     "higher_in_a": int((strong & (d["rbc"] < 0)).sum()),
                     "higher_in_b": int((strong & (d["rbc"] > 0)).sum()),
                     "higher_in_a_any": int((sig & (d["rbc"] < 0)).sum()),
                     "higher_in_b_any": int((sig & (d["rbc"] > 0)).sum()),
                     "median_rbc": float(np.nanmedian(d["rbc"])), "sig_q": q_sig, "min_r": min_r})
        C.note("fraction pair %s / %s: %d samples, %d features, %d higher in %s, %d higher in %s (q < %g, |r| >= %g; "
               "q only: %d / %d)" % (la, lb, summ[-1]["samples"], len(d), summ[-1]["higher_in_a"], la,
                                    summ[-1]["higher_in_b"], lb, q_sig, min_r, summ[-1]["higher_in_a_any"],
                                    summ[-1]["higher_in_b_any"]))
    if not rows:
        return {}
    return {"fraction_diff": pd.concat(rows, ignore_index=True), "fraction_summary": pd.DataFrame(summ)}
