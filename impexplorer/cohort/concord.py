"""Agreement between layers (--concordance): how the omics of the same samples relate, as in a paper's figure on
consistency and individuality of the layers.

Layers are compared in sets of the same feature space, e.g. genera in the metagenome, metatranscriptome, metaproteome
and sRNA, or Pfam domains in mg, mt and the proteins:

  --concordance "genus=kraken.genus.mg,kraken.genus.mt,prot.genus.intra,prot.genus.extra,srna.genus.srna_int;
                 pfam=pfam.mg,pfam.mt,prot.pfam.intra,prot.pfam.extra"
  ('auto' = these two sets with the layers present; srna.<set>.<fraction> = an imported sRNA-IMP profile)

Per set, on the samples two layers share (samples flagged FAIL in a layer's sample_qc left out, as --paper-qc):

  presence   per layer and feature: samples with the feature, prevalence, median relative abundance and median CLR
             (overlap of the layers, singletons, shared by all samples, top features across layers)
  alpha      Spearman correlation of Shannon diversity between every two layers
  mantel     Mantel test (Pearson correlation of the Bray-Curtis distances, samples permuted) between every two layers
  features   Spearman correlation of the same feature between every two layers (features detected in >= --prevalence of
             the shared samples in both)
BH-FDR over the pairs of a set (alpha, mantel) and per pair (features).
"""
import numpy as np
import pandas as pd

from ..release import read_table, table_path
from .paper import qc_drop, qc_layer
from .stats import bh, bray_curtis, clr, rankdata, spearman, t_p

AUTO = {"genus": ["kraken.genus.mg", "kraken.genus.mt", "prot.genus.intra", "prot.genus.extra",
                  "srna.genus.srna_int", "srna.genus.srna_ext"],
        "pfam": ["pfam.mg", "pfam.mt", "prot.pfam.intra", "prot.pfam.extra"]}
MIN_SHARED = 10


def parse(spec):
    """'name=a,b,c; name2=d,e' -> {name: [keys]}; 'auto' -> AUTO"""
    spec = (spec or "").strip()
    if not spec:
        return {}
    if spec == "auto":
        return dict(AUTO)
    out = {}
    for part in spec.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k.strip()] = [x.strip() for x in v.split(",") if x.strip()]
    return out


def needed_layers(sets):
    """cohort layers (not srna.*) the sets need, so the cohort command loads them"""
    return [k for keys in sets.values() for k in keys if not k.startswith("srna.")]


def srna_matrix(release, st, omic):
    if not table_path(release, "srna_features"):
        return None
    d = read_table(release, "srna_features")
    d = d[(d["set"] == st) & (d["omic"] == omic)]
    if d.empty:
        return None
    M = d.pivot_table(index="sample_id", columns="feature", values="reads", aggfunc="sum", fill_value=0)
    M.index = M.index.astype(str)
    return M.astype(float)


def matrices(C, a, keys):
    """{key: samples x features counts / AUC (rows with a total > 0, QC-flagged samples of the layer left out)}"""
    drop = qc_drop(C.release, getattr(a, "paper_qc", "FAIL"))
    out = {}
    for k in keys:
        if k.startswith("srna."):
            _, st, omic = k.split(".", 2)
            M = srna_matrix(C.release, st, omic)
            bad = drop.get(omic, set())
            if M is not None:
                M = M[[s in C.samples.index for s in M.index]]
        else:
            L = C.layers.get(k)
            if L is None or L.kind not in ("counts", "auc"):
                continue
            M = L.M
            ql = qc_layer(L)
            bad = drop.get(ql, set()) if ql else set()
        if M is None or M.empty:
            continue
        M = M.loc[[s for s in M.index if s not in bad]]
        M = M.loc[M.sum(axis=1) > 0, M.sum(axis=0) > 0]
        if M.shape[0] >= MIN_SHARED and M.shape[1] >= 2:
            out[k] = M
    return out


def rel(M):
    return M.div(M.sum(axis=1), axis=0)


def shannon(M):
    R = rel(M)
    with np.errstate(divide="ignore", invalid="ignore"):
        return -(R * np.log(R.where(R > 0, 1))).sum(axis=1)


def mantel(D1, D2, perms=999, seed=1):
    """Pearson correlation of the upper triangles; one-sided p (r >= observed) from permuting the samples of D2"""
    n = D1.shape[0]
    iu = np.triu_indices(n, 1)
    x = D1[iu]
    x = (x - x.mean()) / (x.std() or 1)
    r0 = float(np.mean(x * ((D2[iu] - D2[iu].mean()) / (D2[iu].std() or 1))))
    rng = np.random.default_rng(seed)
    ge = 0
    for _ in range(perms):
        p = rng.permutation(n)
        y = D2[np.ix_(p, p)][iu]
        r = float(np.mean(x * ((y - y.mean()) / (y.std() or 1))))
        ge += r >= r0 - 1e-12
    return r0, (ge + 1) / (perms + 1)


def feature_cor(A, B, prevalence):
    """Spearman of every shared feature between two layers (same samples): rho, p, n"""
    ss = sorted(set(A.index) & set(B.index))
    f = sorted(set(A.columns) & set(B.columns))
    if len(ss) < MIN_SHARED or not f:
        return pd.DataFrame(columns=["feature", "n", "rho", "p"])
    X, Y = rel(A.loc[ss])[f], rel(B.loc[ss])[f]
    keep = [c for c in f if (X[c] > 0).mean() >= prevalence and (Y[c] > 0).mean() >= prevalence]
    if not keep:
        return pd.DataFrame(columns=["feature", "n", "rho", "p"])
    RX = np.apply_along_axis(rankdata, 0, X[keep].values)
    RY = np.apply_along_axis(rankdata, 0, Y[keep].values)
    RX, RY = RX - RX.mean(axis=0), RY - RY.mean(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        rho = (RX * RY).sum(axis=0) / np.sqrt((RX ** 2).sum(axis=0) * (RY ** 2).sum(axis=0))
        n = len(ss)
        T = rho * np.sqrt((n - 2) / np.maximum(1 - rho ** 2, 1e-12))
    return pd.DataFrame({"feature": keep, "n": n, "rho": rho, "p": t_p(T, n - 2)})


def run(C, a):
    sets = parse(a.concordance)
    pres, alpha, man, feat = [], [], [], []
    for name, keys in sets.items():
        Ms = matrices(C, a, keys)
        if len(Ms) < 2:
            C.note("concordance %s: fewer than two layers available (%s)" % (name, ", ".join(Ms) or "none"))
            continue
        ks = [k for k in keys if k in Ms]
        C.note("concordance %s: %s" % (name, ", ".join("%s (%d x %d)" % (k, Ms[k].shape[0], Ms[k].shape[1])
                                                        for k in ks)))
        for k in ks:
            M = Ms[k]
            pos = M.values[M.values > 0]
            Z = pd.DataFrame(clr(M.values, pseudo=1.0 if M.values.max() > 1 else (pos.min() / 2 if len(pos) else 1e-9)),
                             index=M.index, columns=M.columns)
            R = rel(M)
            pres.append(pd.DataFrame({"set": name, "layer": k, "feature": M.columns,
                                      "samples": (M > 0).sum(axis=0).values, "of": M.shape[0],
                                      "prevalence": (M > 0).mean(axis=0).values,
                                      "median_rel": R.median(axis=0).values, "median_clr": Z.median(axis=0).values}))
        H = {k: shannon(Ms[k]) for k in ks}
        rows_a, rows_m = [], []
        for i, k1 in enumerate(ks):
            for k2 in ks[i + 1:]:
                ss = sorted(set(Ms[k1].index) & set(Ms[k2].index))
                if len(ss) < MIN_SHARED:
                    continue
                rho, p, n = spearman(H[k1].loc[ss].values, H[k2].loc[ss].values)
                rows_a.append([name, k1, k2, n, rho, p])
                D1, D2 = bray_curtis(rel(Ms[k1].loc[ss]).values), bray_curtis(rel(Ms[k2].loc[ss]).values)
                r, p = mantel(D1, D2, a.permutations)
                rows_m.append([name, k1, k2, len(ss), r, p])
                fc = feature_cor(Ms[k1], Ms[k2], a.prevalence)
                if len(fc):
                    fc.insert(0, "layer_b", k2)
                    fc.insert(0, "layer_a", k1)
                    fc.insert(0, "set", name)
                    fc["q"] = bh(fc["p"].values)
                    feat.append(fc)
        for rows, out in ((rows_a, alpha), (rows_m, man)):
            if rows:
                q = bh([r[-1] for r in rows])
                out += [r + [qq] for r, qq in zip(rows, q)]
    cols = ["set", "layer_a", "layer_b", "n"]
    return {"concord_presence": pd.concat(pres, ignore_index=True) if pres else pd.DataFrame(),
            "concord_alpha": pd.DataFrame(alpha, columns=cols + ["rho", "p", "q"]),
            "concord_mantel": pd.DataFrame(man, columns=cols + ["r", "p", "q"]),
            "concord_features": pd.concat(feat, ignore_index=True) if feat else pd.DataFrame()}
