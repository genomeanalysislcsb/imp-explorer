"""Helpers shared by the pages: selected samples, groups, matrices, standard plots."""
import dash_ag_grid as dag
import numpy as np
import pandas as pd
import plotly.graph_objects as go

from ..ui import colors, empty_fig, style


def selection(data, group_col, fcol, fval):
    """(samples DataFrame, group Series indexed by sample_id, colour map)"""
    s = data.subset(fcol, fval)
    g = data.group_of(group_col, s)
    return s, g, colors(g.unique())


def ids_param(samples):
    return [list(samples["sample_id"].astype(str))]


def matrix(df, value="rel_abundance", feature="taxon"):
    """long -> samples x features (0 filled)"""
    if df.empty:
        return pd.DataFrame()
    return df.pivot_table(index="sample_id", columns=feature, values=value, aggfunc="sum", fill_value=0.0)


def order_samples(M, groups, by=None):
    """samples sorted by group, then by the first feature column (or `by`)"""
    g = groups.reindex(M.index).fillna("NA")
    key = M[by] if by is not None and by in M.columns else (M.iloc[:, 0] if M.shape[1] else 0)
    return pd.DataFrame({"g": g, "k": key}).sort_values(["g", "k"], ascending=[True, False]).index


def box_by_group(values, groups, cmap, ytitle, title=None, log=False, height=360, pseudo=None):
    """box + points per group; values: Series indexed by sample_id"""
    v = values.dropna()
    if log and pseudo is not None:          # zeros shown at the pseudo-count line
        v = v.clip(lower=pseudo)
    elif log:
        v = v[v > 0]
    if v.empty:
        return empty_fig()
    g = groups.reindex(v.index).fillna("NA")
    fig = go.Figure()
    for lv in sorted(g.unique()):
        x = v[g == lv]
        fig.add_trace(go.Box(y=x.values, name="%s (n=%d)" % (lv, len(x)), boxpoints="all", jitter=0.4,
                             pointpos=0, marker=dict(size=5, color=cmap.get(lv, "#555"), opacity=0.75),
                             line=dict(color=cmap.get(lv, "#555"), width=1.4), fillcolor="rgba(0,0,0,0)",
                             text=list(x.index), hovertemplate="%{text}<br>%{y:.4g}<extra></extra>"))
    fig.update_yaxes(title=ytitle, type="log" if log else "linear")
    fig.update_layout(title=title)
    return style(fig, height, legend=False)


def stacked(M, groups, top=15, ytitle="relative abundance", title=None, height=430):
    """stacked bars per sample of the `top` features (rest = other), samples grouped"""
    if M.empty:
        return empty_fig()
    R = M.div(M.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
    keep = R.mean().sort_values(ascending=False).index[:top]
    R2 = R[keep].copy()
    R2["other"] = 1 - R2.sum(axis=1)
    order = order_samples(R2, groups, keep[0] if len(keep) else None)
    R2 = R2.loc[order]
    pal = ["#2a6f97", "#d1495b", "#edae49", "#00798c", "#8f5aa8", "#66a182", "#e07a5f", "#3d405b", "#81b29a",
           "#f2cc8f", "#6d597a", "#b56576", "#457b9d", "#a8dadc", "#ffb4a2", "#9c6644", "#52796f", "#cad2c5",
           "#e5989b", "#118ab2"]
    fig = go.Figure()
    for i, f in enumerate(R2.columns):
        fig.add_trace(go.Bar(x=list(R2.index), y=R2[f].values, name=str(f),
                             marker_color="#d9dde2" if f == "other" else pal[i % len(pal)],
                             hovertemplate="%{x}<br>" + str(f).replace("%", "%%") + ": %{y:.1%}<extra></extra>"))
    g = groups.reindex(R2.index).fillna("NA")
    # group separators + labels
    pos = 0
    for lv, n in g.groupby(g, sort=False).size().items():
        fig.add_annotation(x=pos + n / 2 - 0.5, y=1.04, xref="x", yref="paper", text="<b>%s</b>" % lv,
                           showarrow=False, font=dict(size=11))
        if pos:
            fig.add_vline(x=pos - 0.5, line=dict(color="#1f2933", width=1))
        pos += n
    fig.update_layout(barmode="stack", bargap=0.08, title=title)
    fig.update_yaxes(title=ytitle, tickformat=".0%", range=[0, 1])
    fig.update_xaxes(showticklabels=len(R2) <= 60, tickangle=-90, tickfont=dict(size=9))
    style(fig, height)
    fig.update_layout(legend=dict(orientation="v", y=1, x=1.01, font=dict(size=10)), margin=dict(r=180, t=50))
    return fig


def ordination(M, groups, cmap, title="PCoA (Bray-Curtis)", height=430):
    from ..stats import bray_curtis, pcoa, permanova
    if M.shape[0] < 3 or M.shape[1] < 2:
        return empty_fig("need at least 3 samples")
    R = M.div(M.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
    D = bray_curtis(R.values)
    X, expl = pcoa(D)
    g = groups.reindex(R.index).fillna("NA")
    F, r2, p = permanova(D, g.values) if g.nunique() > 1 else (np.nan, np.nan, np.nan)
    fig = go.Figure()
    for lv in sorted(g.unique()):
        i = np.where(g.values == lv)[0]
        fig.add_trace(go.Scatter(x=X[i, 0], y=X[i, 1], mode="markers", name=str(lv), text=list(R.index[i]),
                                 marker=dict(size=9, color=cmap.get(lv, "#555"), line=dict(width=0.5, color="#fff")),
                                 hovertemplate="%{text}<extra>" + str(lv) + "</extra>"))
    sub = "" if np.isnan(r2) else "  ·  PERMANOVA R² = %.3f, p = %.3g" % (r2, p)
    fig.update_layout(title=title + sub)
    fig.update_xaxes(title="PCo1 (%.0f%%)" % (100 * expl[0]) if len(expl) else "PCo1", zeroline=False)
    fig.update_yaxes(title="PCo2 (%.0f%%)" % (100 * expl[1]) if len(expl) > 1 else "PCo2", zeroline=False)
    return style(fig, height)


def grid(id_, df, height=380, number_cols=(), page_size=None, fit=False):
    cols = []
    for c in df.columns:
        d = {"field": c, "headerName": c, "filter": True, "sortable": True, "resizable": True}
        if c in number_cols or pd.api.types.is_float_dtype(df[c]):
            d["valueFormatter"] = {"function": "params.value == null ? '' : d3.format('.3~g')(params.value)"}
            d["filter"] = "agNumberColumnFilter"
        if c == "link":
            d.update(cellRenderer="markdown", linkTarget="_blank", filter=False, sortable=False, width=90,
                     minWidth=90)
        if fit:
            d.setdefault("minWidth", 95)
        if fit and c in ("name", "taxon", "term", "kind", "feature", "reasons", "comparison"):
            d["flex"] = 2
            d["minWidth"] = 160
        cols.append(d)
    return dag.AgGrid(id=id_, rowData=df.to_dict("records"), columnDefs=cols, className="ag-theme-alpine",
                      columnSize="responsiveSizeToFit" if fit else "autoSize", style={"height": height},
                      dashGridOptions={"rowSelection": "single", "pagination": bool(page_size),
                                       "paginationPageSize": page_size or 100, "animateRows": False})


# ----------------------------------------------------------------------------- clustered heatmap
HEAT_ROWS = [{"label": "by variance", "value": "variance"}, {"label": "clustered", "value": "cluster"}]
HEAT_COLS = [{"label": "by group, then top feature", "value": "top"},
             {"label": "clustered within group", "value": "within"},
             {"label": "clustered (all samples)", "value": "all"}]
HEAT_VALUES = [{"label": "z of log10 abundance", "value": "log"}, {"label": "z of CLR", "value": "clr"}]
ANNOT_PAL = ["#8f5aa8", "#66a182", "#e07a5f", "#3d405b", "#81b29a", "#f2cc8f", "#6d597a", "#b56576", "#457b9d",
             "#9c6644", "#52796f", "#cad2c5"]


def _linkage(X, metric):
    """average-linkage tree of the rows of X with optimal leaf ordering; None if not possible"""
    from scipy.cluster.hierarchy import linkage
    if X.shape[0] < 3:
        return None
    X = np.nan_to_num(np.asarray(X, float))
    if metric == "correlation":                     # constant rows have no correlation: tiny jitter keeps them
        X = X + np.random.default_rng(0).normal(0, 1e-9, X.shape)
    try:
        return linkage(X, method="average", metric=metric, optimal_ordering=True)
    except Exception:
        return None


def _leaves(Zl, n):
    from scipy.cluster.hierarchy import leaves_list
    return list(leaves_list(Zl)) if Zl is not None else list(range(n))


def _dendro_lines(Zl, orient):
    """dendrogram segments in heatmap coordinates (leaf i at position i); orient 'left' (rows) or 'top' (cols)"""
    from scipy.cluster.hierarchy import dendrogram
    d = dendrogram(Zl, no_plot=True)
    xs, ys = [], []
    for ic, dc in zip(d["icoord"], d["dcoord"]):
        pos = [(v - 5) / 10.0 for v in ic]
        if orient == "left":
            xs += [-h for h in dc] + [None]
            ys += pos + [None]
        else:
            xs += pos + [None]
            ys += list(dc) + [None]
    return xs, ys


def heatmap_matrix(M, values="log", top=40):
    """M: samples x features (relative abundances) -> z-scores of the `top` most variable features"""
    if values == "clr":
        pseudo = M[M > 0].min().min() / 2 if (M > 0).any().any() else 1e-6
        L = np.log(M + pseudo)
        L = L.sub(L.mean(axis=1), axis=0)
    else:
        L = np.log10(M + 1e-6)
    feats = L.var().sort_values(ascending=False).index[:top]
    Z = (L[feats] - L[feats].mean()) / L[feats].std().replace(0, 1)
    return Z, M[feats]


def clustered_heatmap(M, groups, cmap, rows="variance", cols="top", values="log", top=40, annot=None,
                      ylabels=None, height=None, title=None):
    """heatmap of the most variable features (rows) x samples (columns) with optional clustering.
    rows: 'variance' | 'cluster' (correlation distance, average linkage, dendrogram on the left);
    cols: 'top' (group, then the first feature) | 'within' (clustered within each group) | 'all' (clustered over
    all samples, dendrogram on top); annot: {name: Series sample -> value} colour bars above the heatmap (the group
    bar is always shown); ylabels: tick text per feature (e.g. links)."""
    from plotly.subplots import make_subplots
    if M.empty or M.shape[1] < 2:
        return empty_fig()
    Z, Mt = heatmap_matrix(M, values, top)
    g = groups.reindex(Z.index).fillna("NA")
    # feature order
    rZ = _linkage(Z.T.values, "correlation") if rows == "cluster" else None
    forder = [Z.columns[i] for i in _leaves(rZ, Z.shape[1])]
    # sample order
    cZ = None
    if cols == "all":
        cZ = _linkage(Z.values, "euclidean")
        sorder = [Z.index[i] for i in _leaves(cZ, Z.shape[0])]
    elif cols == "within":
        sorder = []
        for lv in sorted(g.unique()):
            ids = list(g.index[g == lv])
            sub = Z.loc[ids]
            sorder += [ids[i] for i in _leaves(_linkage(sub.values, "euclidean"), len(ids))]
    else:
        sorder = list(order_samples(Mt, g, Mt.columns[0]))
    Z = Z.loc[sorder, forder]
    g = g.loc[sorder]
    bars = [("group", g)] + [(k, v.reindex(sorder).fillna("NA").astype(str)) for k, v in (annot or {}).items()]
    nb = len(bars)
    show_cd = cZ is not None
    heights = ([0.12] if show_cd else []) + [0.035] * nb + [1.0]
    fig = make_subplots(rows=len(heights), cols=2, shared_xaxes=True, shared_yaxes=False,
                        row_heights=heights, column_widths=[0.09 if rZ is not None else 0.001, 0.91],
                        horizontal_spacing=0.004, vertical_spacing=0.006)
    hr = len(heights)
    xs = list(range(len(sorder)))          # plain lists: no int64 typed arrays (plotly >= 6 with dash 2)
    if show_cd:
        dx, dy = _dendro_lines(cZ, "top")
        fig.add_trace(go.Scatter(x=dx, y=dy, mode="lines", line=dict(color="#7b8794", width=0.8),
                                 hoverinfo="skip", showlegend=False), row=1, col=2)
        fig.update_yaxes(visible=False, row=1, col=2)
    # annotation bars
    for i, (name, ser) in enumerate(bars):
        r = (2 if show_cd else 1) + i
        levels = sorted(ser.unique())
        pal = {lv: (cmap.get(lv) if name == "group" else ANNOT_PAL[j % len(ANNOT_PAL)]) or "#999"
               for j, lv in enumerate(levels)}
        code = [int(x) for x in ser.map({lv: j for j, lv in enumerate(levels)})]
        k = max(len(levels), 1)
        scale = []
        for j, lv in enumerate(levels):
            scale += [[j / k, pal[lv]], [(j + 1) / k, pal[lv]]]
        fig.add_trace(go.Heatmap(z=[code], x=xs, y=[name], zmin=-0.5, zmax=k - 0.5, colorscale=scale,
                                 showscale=False, customdata=[[f"{s} · {name}: {v}" for s, v in zip(sorder, ser)]],
                                 hovertemplate="%{customdata}<extra></extra>"), row=r, col=2)
        fig.update_yaxes(tickfont=dict(size=10), row=r, col=2)
        fig.add_annotation(text="  ".join('<span style="color:%s">■</span> %s' % (pal[lv], lv) for lv in levels[:12]),
                           xref="paper", yref="y%d" % (2 * r), x=1.005, xanchor="left", y=name, showarrow=False,
                           font=dict(size=10))
    # heatmap
    text = [[f"{f} · {s}" for s in sorder] for f in forder]
    fig.add_trace(go.Heatmap(z=Z.T.values.astype(float).tolist(), x=xs, y=list(range(len(forder))), colorscale="RdBu_r", zmid=0, zmin=-3,
                             zmax=3, colorbar=dict(title=dict(text="z", side="top"), orientation="h", thickness=8, len=160, lenmode="pixels",
                                           x=1.0, xanchor="right", y=0, yanchor="top", ypad=6,
                                           tickfont=dict(size=9)), text=text,
                             hovertemplate="%{text}<br>z %{z:.2f}<extra></extra>"), row=hr, col=2)
    lab = ylabels or {}
    fig.update_yaxes(tickmode="array", tickvals=list(range(len(forder))), ticktext=[lab.get(f, f) for f in forder],
                     tickfont=dict(size=9), autorange="reversed", side="right" if rZ is not None else "left",
                     row=hr, col=2)
    if rZ is not None:
        dx, dy = _dendro_lines(rZ, "left")
        fig.add_trace(go.Scatter(x=dx, y=dy, mode="lines", line=dict(color="#7b8794", width=0.8),
                                 hoverinfo="skip", showlegend=False), row=hr, col=1)
        fig.update_yaxes(range=[len(forder) - 0.5, -0.5], visible=False, row=hr, col=1)
        fig.update_xaxes(visible=False, row=hr, col=1)
    for r in range(1, hr + 1):
        fig.update_xaxes(visible=False, range=[-0.5, len(sorder) - 0.5], row=r, col=2)
        if r < hr:
            fig.update_xaxes(visible=False, row=r, col=1)
            fig.update_yaxes(visible=False, row=r, col=1)
    if cols != "all":                                # group separators and labels
        pos = 0
        for lv, n in g.groupby(g, sort=False).size().items():
            if pos:
                fig.add_vline(x=pos - 0.5, line=dict(color="#1f2933", width=1.5), row=hr, col=2)
            pos += n
    style(fig, height or (220 + 14 * len(forder) + 22 * nb + (70 if show_cd else 0)), legend=False)
    fig.update_layout(title=title, margin=dict(l=10 if rZ is not None else 80, r=230 if rZ is not None else 150,
                                               t=40 if title else 20, b=50))   # z colour bar under the heatmap
    return fig
