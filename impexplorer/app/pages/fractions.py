"""Fractions: intracellular vs extracellular per feature ('imp-explorer cohort --fraction-pairs'): for each pair of
layers of the same samples (e.g. metagenome / extracellular metagenome, intracellular / extracellular proteins,
sRNA-Int / sRNA-Ext), the features enriched in one fraction (paired Wilcoxon on CLR, matched-pairs rank-biserial
correlation, BH-FDR). Shown when the release has these results."""
import numpy as np
import plotly.graph_objects as go
from dash import Input, Output, html

from ..ui import Page, card, control, dd, empty_fig, graph, style
from ._common import grid
from .concordance import short

A_COL, B_COL, NS_COL = "#2a6f97", "#d1495b", "#bcccdc"
TOP = 25


def pairs(data):
    d = data.q("SELECT layer_a, layer_b FROM cohort_fraction_summary")
    return [{"label": "%s → %s" % (short(data, a), short(data, b)), "value": "%s|%s" % (a, b)}
            for a, b in zip(d["layer_a"], d["layer_b"])]


def summary_fig(data):
    d = data.q("SELECT * FROM cohort_fraction_summary")
    if d.empty:
        return empty_fig()
    d["pair"] = ["%s → %s" % (short(data, a), short(data, b)) for a, b in zip(d["layer_a"], d["layer_b"])]
    d = d.iloc[::-1]
    fa, fb = 100 * d["higher_in_a"] / d["features"], 100 * d["higher_in_b"] / d["features"]
    fig = go.Figure([go.Bar(y=d["pair"], x=-fa, orientation="h", name="higher in the first (e.g. intracellular)",
                            marker_color=A_COL, customdata=np.c_[d["higher_in_a"], d["features"]],
                            hovertemplate="%{y}<br>%{customdata[0]} of %{customdata[1]} features<extra></extra>"),
                     go.Bar(y=d["pair"], x=fb, orientation="h", name="higher in the second (e.g. extracellular)",
                            marker_color=B_COL, customdata=np.c_[d["higher_in_b"], d["features"]],
                            hovertemplate="%{y}<br>%{customdata[0]} of %{customdata[1]} features<extra></extra>")])
    m = float(max(fa.max(), fb.max(), 5))
    fig.update_layout(barmode="relative")
    fig.update_xaxes(title="% of the tested features (q < %g)" % d["sig_q"].iloc[0], range=[-m * 1.1, m * 1.1],
                     zeroline=True, zerolinecolor="#52606d")
    fig.update_yaxes(automargin=True, dtick=1)
    style(fig, 180 + 30 * len(d))
    fig.update_layout(margin=dict(l=20, r=20, t=50, b=50), legend=dict(orientation="h", x=0, y=1.0, yanchor="bottom"))
    return fig


def load(data, a, b):
    return data.q("SELECT * FROM cohort_fraction_diff WHERE layer_a = ? AND layer_b = ?", [a, b])


def volcano(d, q, la, lb):
    if d.empty:
        return empty_fig()
    y = -np.log10(d["q"].clip(lower=1e-300))
    col = np.where(d["q"] >= q, NS_COL, np.where(d["rbc"] > 0, B_COL, A_COL))
    fig = go.Figure(go.Scattergl(x=d["rbc"], y=y, mode="markers", text=d["feature"],
                                 marker=dict(size=5, color=col, opacity=0.75),
                                 customdata=np.c_[d["median_log2_ratio"], d["n"]],
                                 hovertemplate="%{text}<br>r = %{x:.2f}, q = 10^-%{y:.1f}<br>median log2 ratio "
                                               "%{customdata[0]:.2f}, n = %{customdata[1]}<extra></extra>"))
    fig.add_hline(y=-np.log10(q), line=dict(color="#9aa5b1", dash="dot"))
    fig.update_xaxes(title="rank-biserial r   (← higher in %s · higher in %s →)" % (la, lb), range=[-1.05, 1.05])
    fig.update_yaxes(title="−log10 q")
    style(fig, 460, legend=False)
    return fig


def top_fig(d, q, la, lb, top=TOP):
    s = d[d["q"] < q]
    if s.empty:
        return empty_fig("no feature with q < %g" % q)
    s = s.reindex(s["rbc"].abs().sort_values(ascending=False).index)
    s = s.groupby(s["rbc"] > 0, group_keys=False).head(top).sort_values("rbc")
    fig = go.Figure(go.Bar(y=s["feature"], x=s["rbc"], orientation="h",
                           marker_color=np.where(s["rbc"] > 0, B_COL, A_COL),
                           customdata=np.c_[s["median_log2_ratio"], s["q"]],
                           hovertemplate="%{y}: r = %{x:.2f}, median log2 ratio %{customdata[0]:.2f}, "
                                         "q = %{customdata[1]:.2g}<extra></extra>"))
    fig.update_xaxes(title="rank-biserial r  (← %s · %s →)" % (la, lb), range=[-1.05, 1.05])
    fig.update_yaxes(automargin=True, dtick=1, tickfont=dict(size=10))
    style(fig, 140 + 16 * len(s), legend=False)
    fig.update_layout(margin=dict(l=20, r=20, t=20, b=50))
    return fig


def layout(data):
    ps = pairs(data)
    q = float(data.q("SELECT max(sig_q) AS q FROM cohort_fraction_summary")["q"].iloc[0] or 0.1)
    return html.Div([
        html.H2("Intracellular vs extracellular"),
        html.P("Per feature, which fraction of the same samples has more of it ('imp-explorer cohort "
               "--fraction-pairs'): paired Wilcoxon signed-rank test of the CLR values (each fraction closed to the "
               "features of the pair), matched-pairs rank-biserial correlation r as effect size (−1 all samples higher "
               "in the first fraction, +1 all higher in the second), BH-FDR per pair (q < %g). Samples flagged FAIL "
               "in a layer's QC are left out." % q, className="lead"),
        card("Features enriched per fraction", graph("fr-summary"),
             note="share of each pair's tested features significantly higher in the first (left) or the second "
                  "(right) fraction"),
        html.Div([control("pair of layers", dd("fr-pair", ps, ps[0]["value"] if ps else None))],
                 className="controls"),
        html.Div([card("Volcano", graph("fr-volcano"), note="each point a feature; dotted line q = %g" % q),
                  card("Strongest differences", graph("fr-top"),
                       note="the %d features with the largest |r| on each side (q < %g)" % (TOP, q))],
                 className="grid g2"),
        card("Features of the pair", html.Div(id="fr-table"),
             note="prevalence and median relative abundance per fraction, median log2 ratio (second / first), r, "
                  "p and q; sorted by q"),
    ])


def register(app, data):
    q = float(data.q("SELECT max(sig_q) AS q FROM cohort_fraction_summary")["q"].iloc[0] or 0.1) \
        if data.has("cohort_fraction_summary") else 0.1

    @app.callback(Output("fr-summary", "figure"), Input("fr-pair", "id"))
    def _summary(_):
        return summary_fig(data)

    @app.callback(Output("fr-volcano", "figure"), Output("fr-top", "figure"), Output("fr-table", "children"),
                  Input("fr-pair", "value"))
    def _pair(pair):
        if not pair:
            return empty_fig(), empty_fig(), html.P("no pair", className="hint")
        a, b = pair.split("|")
        d = load(data, a, b)
        la, lb = short(data, a), short(data, b)
        t = d[["feature", "n", "prev_a", "prev_b", "median_rel_a", "median_rel_b", "median_log2_ratio", "rbc", "p",
               "q"]].sort_values(["q", "p"]).head(2000)
        t = t.rename(columns={"rbc": "r"})
        return (volcano(d, q, la, lb), top_fig(d, q, la, lb),
                grid("fr-grid", t.round(5), height=420, fit=True) if len(t) else html.P("no features", className="hint"))


PAGE = Page(path="/fractions", name="Intra vs extra", order=73, layout=layout, register=register,
            requires=("cohort_fraction_summary",))
