"""Variance explained: results of the paper analysis ('imp-explorer cohort --paper-model'): how much of each layer's
variance each term of the model explains (Shannon, beta diversity, features), which terms matter most across layers,
how many features each term affects, and how much stays unexplained. Shown when the release has these results."""
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from dash import Input, Output, html

from ..ui import PALETTE, Page, card, control, dd, empty_fig, graph, style
from ._common import grid
from .cohort import layer_label

MEASURES = [("r2_beta", "beta diversity (PERMANOVA R²)", "p_beta"),
            ("r2_shannon", "Shannon diversity (ANOVA r²)", "p_shannon"),
            ("r2_features_median", "features (median r² per feature)", None),
            ("features_q", "features with q < threshold (share)", None)]
KIND_COLOR = {"taxa": "#8f5aa8", "mags": "#d1495b", "functions": "#2a6f97", "prot": "#66a182",
              "metab": "#edae49", "ratio": "#7b8794"}
KINDS = [("taxa", "taxa"), ("mags", "MAGs"), ("functions", "functions"), ("prot", "proteins"),
         ("metab", "metabolites"), ("ratio", "mt/mg ratios")]


def kind_of(layer):
    p = str(layer).split(".")[0]
    if p in ("kraken", "motus"):
        return "taxa"
    if p == "mag":
        return "mags"
    if p in ("prot", "metab", "ratio"):
        return p
    return "functions"


def terms(data):
    model = data.manifest.get("cohort", {}).get("paper", {}).get("model", "")
    order = [t.strip() for t in model.split("+") if t.strip()]
    have = data.distinct("cohort_paper_variance", "term")
    return [t for t in order if t in have] + [t for t in have if t not in order and t != "(residual)"]


def table(data, kinds):
    v = data.q("SELECT * FROM cohort_paper_variance")
    v["kind"] = v["layer"].map(kind_of)
    if kinds:
        v = v[v["kind"].isin(kinds)]
    order = list(data.manifest.get("cohort", {}).get("layers", {}))
    v["o"] = v["layer"].map({k: i for i, k in enumerate(order)}).fillna(len(order))
    return v.sort_values(["o", "layer"])


def heatmap(data, v, measure, terms_):
    mcol, mlab, pcol = next(m for m in MEASURES if m[0] == measure)
    layers = list(dict.fromkeys(v["layer"]))
    Z = v.pivot_table(index="layer", columns="term", values=mcol, aggfunc="first").reindex(index=layers,
                                                                                          columns=terms_)
    if Z.dropna(how="all").empty:
        return empty_fig("no %s for these layers" % mlab)
    P = v.pivot_table(index="layer", columns="term", values=pcol, aggfunc="first").reindex(
        index=layers, columns=terms_) if pcol else None
    txt = []
    for ly in layers:
        row = []
        for t in terms_:
            x = Z.at[ly, t]
            s = "–" if pd.isna(x) else "%.1f" % (100 * x)       # – : term not in the layer's model
            if P is not None and not pd.isna(P.at[ly, t]):
                s += "*" if P.at[ly, t] < 0.05 else ""
            row.append(s)
        txt.append(row)
    labels = ["%s  [%s]" % (layer_label(data, ly), ly) for ly in layers]
    fig = go.Figure(go.Heatmap(z=(100 * Z).values, x=terms_, y=labels, text=txt, texttemplate="%{text}",
                               textfont=dict(size=10), colorscale="Blues", zmin=0,
                               colorbar=dict(title="%", thickness=12),
                               hovertemplate="%{y}<br>%{x}: %{z:.2f} %<extra></extra>"))
    fig.update_yaxes(autorange="reversed", tickfont=dict(size=10))
    fig.update_xaxes(side="top")
    style(fig, 120 + 22 * len(layers), legend=False)
    fig.update_layout(margin=dict(l=320, r=20, t=60, b=20))
    return fig


def ranking(v, measure, terms_):
    mcol, mlab, _ = next(m for m in MEASURES if m[0] == measure)
    d = v[v["term"].isin(terms_)].dropna(subset=[mcol])
    if d.empty:
        return empty_fig()
    med = d.groupby("term")[mcol].median().sort_values(ascending=False)
    fig = go.Figure()
    kinds = sorted(d["kind"].unique())
    col = {k: KIND_COLOR.get(k, PALETTE[i % len(PALETTE)]) for i, k in enumerate(kinds)}
    for k in kinds:
        x = d[d["kind"] == k]
        fig.add_trace(go.Box(x=x["term"], y=100 * x[mcol], name=dict(KINDS).get(k, k), boxpoints="all", jitter=0.5,
                             pointpos=0, marker=dict(size=5, color=col[k], opacity=0.75), line=dict(width=0),
                             fillcolor="rgba(0,0,0,0)", text=x["layer"], showlegend=False, legendgroup=k,
                             hovertemplate="%{text}<br>%{x}: %{y:.2f} %<extra></extra>"))
        fig.add_trace(go.Scatter(x=[None], y=[None], mode="markers", name=dict(KINDS).get(k, k), legendgroup=k,
                                 marker=dict(size=8, color=col[k])))          # legend entry with the colour
    fig.add_trace(go.Scatter(x=list(med.index), y=100 * med.values, mode="markers", name="median",
                             marker=dict(symbol="line-ew-open", size=28, color="#1f2933", line=dict(width=2)),
                             hovertemplate="median %{y:.2f} %<extra></extra>"))
    fig.update_xaxes(categoryorder="array", categoryarray=list(med.index))
    fig.update_yaxes(title="% explained, " + mlab.split(" (")[0])
    fig.update_layout(title="Terms ranked by the median over layers")
    return style(fig, 400)


def unexplained(data, v, measure):
    mcol, mlab, _ = next(m for m in MEASURES if m[0] == measure)
    if mcol == "features_q":
        mcol = "r2_features_median"
    d = v[v["term"] == "(residual)"].dropna(subset=[mcol])
    if d.empty:
        return empty_fig()
    d = d.sort_values(mcol)
    fig = go.Figure(go.Bar(x=100 * d[mcol], y=[layer_label(data, x) for x in d["layer"]], orientation="h",
                           marker_color="#9aa5b1", text=["%.0f %%" % (100 * x) for x in d[mcol]],
                           textposition="outside", hovertemplate="%{y}: %{x:.1f} %<extra></extra>"))
    fig.update_xaxes(title="% unexplained (residual)", range=[0, 118], tickvals=[0, 25, 50, 75, 100])
    fig.update_yaxes(tickfont=dict(size=10))
    style(fig, 120 + 20 * len(d), legend=False)
    fig.update_layout(margin=dict(l=300, r=40, t=30, b=40))
    return fig


def sig_features(data, layer, term, sig_q):
    """share of features with p < 0.05 and with q < sig_q per term (Fig. 4C), and the hits of one term"""
    s = data.q("SELECT term, avg(CASE WHEN p < 0.05 THEN 1.0 ELSE 0 END) AS p05, "
               "avg(CASE WHEN q < ? THEN 1.0 ELSE 0 END) AS q_sig, count(*) AS n FROM cohort_paper_anova "
               "WHERE layer = ? AND analysis = 'feature' AND term <> '(residual)' GROUP BY term", [sig_q, layer])
    hits = data.q("SELECT feature, r2, F, p, q, n FROM cohort_paper_anova WHERE layer = ? AND "
                  "analysis = 'feature' AND term = ? ORDER BY p LIMIT 500", [layer, term]) if term else None
    return s, hits


def layout(data):
    st = data.manifest.get("cohort", {})
    pa = st.get("paper", {})
    sig = st.get("sig_q", 0.1)
    layers = list(dict.fromkeys(data.distinct("cohort_paper_variance", "layer")))
    have_k = sorted({kind_of(x) for x in layers})
    t = terms(data)
    return html.Div([
        html.H2("Variance explained"),
        html.P(["Paper analysis of 'imp-explorer cohort --paper-model': one linear model per layer, ",
                html.Code(pa.get("model", "")), ". Type II ANOVA per term for Shannon diversity and every feature, "
                "PERMANOVA with all terms (marginal R²); the share of variance each term explains given the others. "
                "Samples flagged %s in the layer's QC are left out." % pa.get("qc", "FAIL")], className="lead"),
        html.Div([control("measure", dd("va-measure", [{"label": lab, "value": m} for m, lab, _ in MEASURES],
                                        "r2_beta")),
                  control("layers", dd("va-kinds", [{"label": dict(KINDS).get(k, k), "value": k} for k in have_k],
                                       have_k, multi=True))], className="controls"),
        card("Variance explained per layer and term", graph("va-heat"),
             note="% of the variance explained by each term given all others; * p < 0.05 (beta diversity, "
                  "Shannon); – term not in that layer's model (depth: count and protein layers only). Correlated "
                  "terms share variance that neither gets."),
        html.Div([card("Terms across layers", graph("va-rank"), note="each point a layer; terms sorted by the "
                       "median"),
                  card("Unexplained variance", graph("va-resid"), note="residual share per layer (for the feature "
                       "measures: median over features)")], className="grid g2"),
        card("Features per term", html.Div([
            html.Div([control("layer", dd("va-layer", [{"label": "%s [%s]" % (layer_label(data, x), x), "value": x}
                                                       for x in layers], layers[0] if layers else None)),
                      control("term", dd("va-term", t, t[0] if t else None))], className="controls"),
            html.Div([graph("va-sig"), html.Div(id="va-hits")], className="grid g2")]),
            note="share of the layer's features with p < 0.05 and with q < %g per term; the features of one term "
                 "by p" % sig),
    ])


def register(app, data):
    sig = data.manifest.get("cohort", {}).get("sig_q", 0.1)

    @app.callback(Output("va-heat", "figure"), Output("va-rank", "figure"), Output("va-resid", "figure"),
                  Input("va-measure", "value"), Input("va-kinds", "value"))
    def _main(measure, kinds):
        v = table(data, kinds or [])
        if v.empty:
            return empty_fig(), empty_fig(), empty_fig()
        t = terms(data)
        return heatmap(data, v, measure, t), ranking(v, measure, t), unexplained(data, v, measure)

    @app.callback(Output("va-sig", "figure"), Output("va-hits", "children"), Input("va-layer", "value"),
                  Input("va-term", "value"))
    def _features(layer, term):
        if not layer:
            return empty_fig(), None
        s, hits = sig_features(data, layer, term, sig)
        if s.empty:
            return empty_fig("no feature tests for this layer"), None
        order = [t for t in terms(data) if t in set(s["term"])]
        s = s.set_index("term").reindex(order)
        fig = go.Figure([go.Bar(x=order, y=100 * s["p05"], name="p < 0.05", marker_color="#8fb8d0"),
                         go.Bar(x=order, y=100 * s["q_sig"], name="q < %g" % sig, marker_color="#2a6f97")])
        fig.update_yaxes(title="%% of %d features" % int(s["n"].max()))
        fig.update_layout(barmode="group", title=layer_label(data, layer))
        tab = grid("va-grid", hits.round(5), height=380, fit=True) if hits is not None and len(hits) else \
            html.P("no features", className="hint")
        return style(fig, 380), tab


PAGE = Page(path="/variance", name="Variance", order=75, layout=layout, register=register,
            requires=("cohort_paper_variance",))
