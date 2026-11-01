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
