"""Concordance: agreement between the layers of the same samples ('imp-explorer cohort --concordance'): overlap of
the features, top features across layers, sparsity, Shannon correlation, Mantel tests and correlations of the same
feature between layers. Shown when the release has these results."""
import itertools

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from dash import Input, Output, html
from plotly.subplots import make_subplots

from ..ui import Page, card, control, dd, empty_fig, graph, style
from ._common import grid
from .cohort import layer_label

TOP = 30


def label(data, key):
    if key.startswith("srna."):
        _, st, omic = key.split(".", 2)
        fr = data.manifest.get("srna", {}).get("fractions", {})
        return "sRNA %s, %s" % (fr.get(omic, omic), st)
    return layer_label(data, key)


def short(data, key):
    return label(data, key).replace("Kraken2/Bracken ", "").replace("Metaproteome ", "proteins ")


def layers_of(data, s):
    d = data.q("SELECT DISTINCT layer FROM cohort_concord_presence WHERE \"set\" = ?", [s])["layer"].tolist()
    order = [k for keys in [data.manifest.get("cohort", {}).get("concordance", "")] for k in
             keys.replace(";", ",").replace("=", ",").split(",")]
    return [k for k in order if k in d] + [k for k in d if k not in order]


def pair_matrix(data, table, col, s, keys):
    d = data.q("SELECT * FROM %s WHERE \"set\" = ?" % table, [s])
    n = len(keys)
    Z = np.full((n, n), np.nan)
    T = [["" for _ in keys] for _ in keys]
    idx = {k: i for i, k in enumerate(keys)}
    for r in d.itertuples(index=False):
        i, j = idx.get(r.layer_a), idx.get(r.layer_b)
        if i is None or j is None:
            continue
        v, q = getattr(r, col), r.q
        for a, b in ((i, j), (j, i)):
            Z[a, b] = v
            T[a][b] = "%.2f%s" % (v, "*" if q < 0.05 else "") if v == v else ""
    for i in range(n):
        Z[i, i] = 1.0
        T[i][i] = ""
    labs = [short(data, k) for k in keys]
    fig = go.Figure(go.Heatmap(z=Z, x=labs, y=labs, text=T, texttemplate="%{text}", zmin=-1, zmax=1,
                               colorscale="RdBu", reversescale=True, colorbar=dict(thickness=12),
                               hovertemplate="%{y} · %{x}: %{z:.3f}<extra></extra>"))
    fig.update_yaxes(autorange="reversed", dtick=1, automargin=True)
    fig.update_xaxes(dtick=1, tickangle=-45, automargin=True)
    style(fig, 200 + 34 * n, legend=False)
    fig.update_layout(margin=dict(l=200, r=20, t=20, b=140))
    return fig


def upset(data, s, keys, min_samples=2, top=15):
    """overlap of the layers: features detected in >= min_samples samples of a layer; the most frequent
    combinations of layers (bars) and which layers they contain (dots)"""
    p = data.q("SELECT layer, feature FROM cohort_concord_presence WHERE \"set\" = ? AND samples >= ?",
               [s, min_samples])
    if p.empty:
        return empty_fig()
    combo = p.groupby("feature")["layer"].apply(lambda x: tuple(k for k in keys if k in set(x)))
    cnt = combo.value_counts().head(top)
    fig = make_subplots(rows=2, cols=1, shared_xaxes=True, row_heights=[0.62, 0.38], vertical_spacing=0.02)
    xs = list(range(len(cnt)))
    fig.add_trace(go.Bar(x=xs, y=cnt.values, marker_color="#2a6f97", text=cnt.values, textposition="outside",
                         hovertemplate="%{y} features<extra></extra>"), row=1, col=1)
    for j, k in enumerate(keys):
        on = [k in c for c in cnt.index]
        fig.add_trace(go.Scatter(x=xs, y=[j] * len(xs), mode="markers",
                                 marker=dict(size=10, color=["#1f2933" if o else "#e4e7eb" for o in on]),
                                 hoverinfo="skip"), row=2, col=1)
    for x, c in zip(xs, cnt.index):
        js = [j for j, k in enumerate(keys) if k in c]
        if len(js) > 1:
            fig.add_trace(go.Scatter(x=[x, x], y=[min(js), max(js)], mode="lines", line=dict(color="#1f2933"),
                                     hoverinfo="skip"), row=2, col=1)
    fig.update_yaxes(tickvals=list(range(len(keys))), ticktext=[short(data, k) for k in keys], autorange="reversed",
                     row=2, col=1)
    fig.update_yaxes(title="features", row=1, col=1)
    fig.update_xaxes(visible=False)
    style(fig, 300 + 26 * len(keys), legend=False)
    fig.update_layout(margin=dict(l=200, r=20, t=20, b=20))
    return fig


def top_features(data, s, keys, top=TOP):
    """the features ranked highest across the layers: per layer the percentile of the feature's median relative
    abundance among that layer's features (100 = the most abundant), so layers of different depth, sparsity and
    scale (counts, protein AUC) are comparable; a feature missing from a layer is blank"""
    p = data.q("SELECT layer, feature, median_rel, median_clr FROM cohort_concord_presence WHERE \"set\" = ?", [s])
    if p.empty:
        return empty_fig()
    p["pct"] = 100 * p.groupby("layer")["median_clr"].rank(pct=True)
    Z = p.pivot_table(index="feature", columns="layer", values="pct").reindex(columns=keys)
    order = Z.fillna(0).median(axis=1).sort_values(ascending=False).index[:top]
    Z = Z.loc[order]
    R = p.pivot_table(index="feature", columns="layer", values="median_rel").reindex(index=order, columns=keys)
    fig = go.Figure(go.Heatmap(z=Z.values, x=[short(data, k) for k in keys], y=list(Z.index), colorscale="Viridis",
                               zmin=0, zmax=100, customdata=100 * R.values,
                               colorbar=dict(title="percentile<br>in the layer", thickness=12),
                               hovertemplate="%{y} · %{x}: percentile %{z:.0f}, median %{customdata:.3g} %"
                                             "<extra></extra>"))
    fig.update_yaxes(autorange="reversed", tickfont=dict(size=10), dtick=1, automargin=True)
    fig.update_xaxes(dtick=1, tickangle=-45, automargin=True)
    style(fig, 120 + 18 * len(order), legend=False)
    fig.update_layout(margin=dict(l=200, r=20, t=20, b=140))
    return fig


def sparsity(data, s, keys):
    p = data.q("SELECT layer, samples, \"of\", prevalence FROM cohort_concord_presence WHERE \"set\" = ?", [s])
    if p.empty:
        return empty_fig()
    g = p.groupby("layer").agg(all_=("prevalence", lambda x: float((x >= 1).mean())),
                               single=("samples", lambda x: float((x == 1).mean())),
                               n=("samples", "size")).reindex(keys)
    labs = [short(data, k) for k in keys]
    fig = go.Figure([go.Bar(x=labs, y=100 * g["all_"], name="in every sample", marker_color="#2a6f97"),
                     go.Bar(x=labs, y=100 * g["single"], name="in one sample only (singleton)",
                            marker_color="#edae49")])
    fig.update_yaxes(title="% of the layer's features")
    fig.update_xaxes(tickangle=-45, automargin=True)
    fig.update_layout(barmode="group")
    style(fig, 420)
    fig.update_layout(margin=dict(b=40, t=50), legend=dict(orientation="h", x=0, y=1.02, yanchor="bottom"))
    return fig


def feature_summary(data, s, keys, q=0.05):
    d = data.q("SELECT layer_a, layer_b, count(*) AS n, avg(CASE WHEN q < ? AND rho > 0 THEN 1.0 ELSE 0 END) AS pos, "
               "avg(CASE WHEN q < ? AND rho < 0 THEN 1.0 ELSE 0 END) AS neg, median(rho) AS med "
               "FROM cohort_concord_features WHERE \"set\" = ? GROUP BY 1, 2", [q, q, s])
    if d.empty:
        return empty_fig("no shared features")
    d["pair"] = [short(data, a) + " · " + short(data, b) for a, b in zip(d["layer_a"], d["layer_b"])]
    d = d.sort_values("pos")
    fig = go.Figure([go.Bar(y=d["pair"], x=100 * d["pos"], orientation="h", name="positive (q < %g)" % q,
                            marker_color="#2a6f97", customdata=d["n"],
                            hovertemplate="%{y}<br>%{x:.0f} % of %{customdata} features<extra></extra>"),
                     go.Bar(y=d["pair"], x=100 * d["neg"], orientation="h", name="negative (q < %g)" % q,
                            marker_color="#d1495b", customdata=d["n"],
                            hovertemplate="%{y}<br>%{x:.0f} % of %{customdata} features<extra></extra>")])
    fig.update_layout(barmode="stack")
    fig.update_xaxes(title="% of the shared features (Spearman)", range=[0, 100], showticklabels=True,
                     side="bottom")
    fig.update_yaxes(automargin=True, dtick=1)
    style(fig, 160 + 24 * len(d))
    fig.update_layout(margin=dict(l=20, r=20, t=50, b=50), legend=dict(orientation="h", x=0, y=1.0,
                                                                      yanchor="bottom"))
    return fig


def within_fig(data, s):
    """per pair of layers: Spearman of the two profiles of the same sample, against the median with the other
    samples (individuality)"""
    d = data.q("SELECT layer_a, layer_b, sample_id, n_features, rho, rho_others FROM cohort_concord_within "
               "WHERE \"set\" = ?", [s])
    if d.empty:
        return empty_fig("no shared samples / features")
    d["pair"] = [short(data, a) + " · " + short(data, b) for a, b in zip(d["layer_a"], d["layer_b"])]
    order = list(d.groupby("pair")["rho"].median().sort_values(ascending=False).index)
    fig = go.Figure()
    for col, name, color in (("rho", "same sample", "#2a6f97"), ("rho_others", "other samples (median)", "#9aa5b1")):
        fig.add_trace(go.Box(x=d[col], y=d["pair"], name=name, orientation="h", boxpoints="all", jitter=0.4,
                             pointpos=0, marker=dict(size=3, color=color, opacity=0.6), line=dict(color=color),
                             fillcolor="rgba(0,0,0,0)", text=d["sample_id"],
                             hovertemplate="%{text}: ρ = %{x:.2f}<extra>" + name + "</extra>"))
    fig.update_layout(boxmode="group")
    fig.update_yaxes(categoryorder="array", categoryarray=order[::-1])
    fig.update_xaxes(title="Spearman ρ of the profiles", range=[-0.2, 1], showticklabels=True, mirror="allticks",
                     side="top")
    fig.update_yaxes(automargin=True)
    style(fig, 200 + 46 * len(order))
    fig.update_layout(margin=dict(l=20, r=20, t=60, b=40), legend=dict(orientation="h", x=0, y=-0.02,
                                                                      yanchor="top"))
    return fig


def layout(data):
    sets = data.distinct("cohort_concord_presence", '"set"')
    s0 = sets[0] if sets else None
    keys = layers_of(data, s0) if s0 else []
    pairs = [{"label": "%s · %s" % (short(data, a), short(data, b)), "value": "%s|%s" % (a, b)}
             for a, b in itertools.combinations(keys, 2)]
    return html.Div([
        html.H2("Concordance of the layers"),
        html.P("How the omics of the same samples agree ('imp-explorer cohort --concordance'): overlap of the "
               "features, top features across layers, sparsity, correlation of Shannon diversity, Mantel tests of the "
               "Bray-Curtis distances, and correlation of the same feature between layers. * q < 0.05. Samples "
               "flagged FAIL in a layer's QC are left out.", className="lead"),
        html.Div([control("feature set", dd("cc-set", sets, s0))], className="controls"),
        html.Div([card("Overlap of the layers", graph("cc-upset"),
                       note="features detected in at least 2 samples of a layer; the most frequent combinations"),
                  card("Top features across layers", graph("cc-top"),
                       note="the %d features ranked highest across the layers: percentile of the feature's median "
                            "abundance within each layer (layers of different depth and scale made comparable)"
                            % TOP)],
                 className="grid g2"),
        html.Div([card("Shannon diversity", graph("cc-alpha"), note="Spearman ρ between layers"),
                  card("Mantel tests", graph("cc-mantel"), note="Pearson r of the Bray-Curtis distances")],
                 className="grid g2"),
        card("Sparsity", graph("cc-sparse"), note="features present in every sample and singletons"),
        card("Same feature in two layers", graph("cc-feat"),
             note="share of the shared features with a significant Spearman correlation across samples (blue "
                  "positive, red negative)"),
        card("Within each sample", graph("cc-within"),
             note="Spearman ρ of the two profiles of the same sample over their shared features (blue), and of its "
                  "profile with the other samples' profiles in the other layer (grey, median): blue above grey = "
                  "the layers keep the individual's signature"),
        card("Features of one pair of layers", html.Div([
            html.Div([control("layers", dd("cc-pair", pairs, pairs[0]["value"] if pairs else None))],
                     className="controls"),
            html.Div(id="cc-feat-table")]), note="Spearman ρ of each shared feature, sorted by q"),
    ])


def register(app, data):
    @app.callback(Output("cc-upset", "figure"), Output("cc-top", "figure"), Output("cc-alpha", "figure"),
                  Output("cc-mantel", "figure"), Output("cc-sparse", "figure"), Output("cc-feat", "figure"),
                  Output("cc-pair", "options"), Output("cc-within", "figure"), Input("cc-set", "value"))
    def _set(s):
        if not s:
            return (empty_fig(),) * 6 + ([], empty_fig())
        keys = layers_of(data, s)
        pairs = [{"label": "%s · %s" % (short(data, a), short(data, b)), "value": "%s|%s" % (a, b)}
                 for a, b in itertools.combinations(keys, 2)]
        return (upset(data, s, keys), top_features(data, s, keys),
                pair_matrix(data, "cohort_concord_alpha", "rho", s, keys) if data.has("cohort_concord_alpha")
                else empty_fig(),
                pair_matrix(data, "cohort_concord_mantel", "r", s, keys) if data.has("cohort_concord_mantel")
                else empty_fig(),
                sparsity(data, s, keys),
                feature_summary(data, s, keys) if data.has("cohort_concord_features") else empty_fig(), pairs,
                within_fig(data, s) if data.has("cohort_concord_within") else empty_fig())

    @app.callback(Output("cc-feat-table", "children"), Input("cc-pair", "value"), Input("cc-set", "value"))
    def _pair(pair, s):
        if not pair or not s or not data.has("cohort_concord_features"):
            return html.P("no pair", className="hint")
        a, b = pair.split("|")
        d = data.q("SELECT feature, n, rho, p, q FROM cohort_concord_features WHERE \"set\" = ? AND layer_a = ? AND "
                   "layer_b = ? ORDER BY q, p LIMIT 1000", [s, a, b])
        if d.empty:
            d = data.q("SELECT feature, n, rho, p, q FROM cohort_concord_features WHERE \"set\" = ? AND layer_a = ? "
                       "AND layer_b = ? ORDER BY q, p LIMIT 1000", [s, b, a])
        return grid("cc-grid", d.round(5), height=400, fit=True) if len(d) else html.P("no shared features",
                                                                                         className="hint")


PAGE = Page(path="/concordance", name="Concordance", order=72, layout=layout, register=register,
            requires=("cohort_concord_presence",))
