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


MIN_R = [0, 0.3, 0.5, 0.7, 0.9]


def default_min_r(data):
    d = data.q("SELECT * FROM cohort_fraction_summary")
    return float(d["min_r"].iloc[0]) if "min_r" in d.columns and len(d) else 0.5


def summary_fig(data, q, min_r):
    """per pair the share of tested features enriched in each fraction: q < q and |r| >= min_r"""
    d = data.q("SELECT layer_a, layer_b, count(*) AS features, "
               "sum(CASE WHEN q < ? AND rbc <= -? THEN 1 ELSE 0 END) AS higher_in_a, "
               "sum(CASE WHEN q < ? AND rbc >= ? THEN 1 ELSE 0 END) AS higher_in_b "
               "FROM cohort_fraction_diff GROUP BY 1, 2", [q, min_r, q, min_r])
    if d.empty:
        return empty_fig()
    order = data.q("SELECT layer_a, layer_b FROM cohort_fraction_summary")
    d = order.merge(d, on=["layer_a", "layer_b"], how="left").fillna(0)
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
    fig.update_xaxes(title="%% of the tested features (q < %g, |r| ≥ %g)" % (q, min_r), range=[-m * 1.1, m * 1.1],
                     zeroline=True, zerolinecolor="#52606d")
    fig.update_yaxes(automargin=True, dtick=1)
    style(fig, 180 + 30 * len(d))
    fig.update_layout(margin=dict(l=20, r=20, t=50, b=50), legend=dict(orientation="h", x=0, y=1.0, yanchor="bottom"))
    return fig


def load(data, a, b):
    return data.q("SELECT * FROM cohort_fraction_diff WHERE layer_a = ? AND layer_b = ?", [a, b])


def volcano(d, q, la, lb, min_r=0.0):
    if d.empty:
        return empty_fig()
    y = -np.log10(d["q"].clip(lower=1e-300))
    weak = (d["q"] >= q) | (d["rbc"].abs() < min_r)
    col = np.where(weak, NS_COL, np.where(d["rbc"] > 0, B_COL, A_COL))
    # point size = abundance: median relative abundance in the fraction where the feature is higher (1e-6 .. 1e-1)
    ab = np.maximum(d["median_rel_a"].fillna(0), d["median_rel_b"].fillna(0)).astype(float)
    size = 3 + 3 * np.clip(np.log10(np.maximum(ab, 1e-6)) + 6, 0, 5)
    fig = go.Figure(go.Scattergl(x=d["rbc"], y=y, mode="markers", text=d["feature"],
                                 marker=dict(size=size, color=col, opacity=0.7, line=dict(width=0)),
                                 customdata=np.c_[d["median_log2_ratio"], d["n"], 100 * ab],
                                 hovertemplate="%{text}<br>r = %{x:.2f}, q = 10^-%{y:.1f}<br>median log2 ratio "
                                               "%{customdata[0]:.2f}, n = %{customdata[1]}<br>median abundance "
                                               "%{customdata[2]:.3g} %<extra></extra>"))
    fig.add_hline(y=-np.log10(q), line=dict(color="#9aa5b1", dash="dot"))
    if min_r > 0:
        for x in (-min_r, min_r):
            fig.add_vline(x=x, line=dict(color="#9aa5b1", dash="dot"))
    fig.update_xaxes(title="rank-biserial r   (← higher in %s · higher in %s →)" % (la, lb), range=[-1.05, 1.05])
    fig.update_yaxes(title="−log10 q")
    style(fig, 460, legend=False)
    return fig


def top_fig(d, q, la, lb, min_r=0.0, top=TOP):
    s = d[(d["q"] < q) & (d["rbc"].abs() >= min_r)]
    if s.empty:
        return empty_fig("no feature with q < %g and |r| ≥ %g" % (q, min_r))
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
    r0 = default_min_r(data)
    q = float(data.q("SELECT max(sig_q) AS q FROM cohort_fraction_summary")["q"].iloc[0] or 0.1)
    return html.Div([
        html.H2("Intracellular vs extracellular"),
        html.P("Per feature, which fraction of the same samples has more of it ('imp-explorer cohort "
               "--fraction-pairs'): paired Wilcoxon signed-rank test of the CLR values (each fraction closed to the "
               "features of the pair), matched-pairs rank-biserial correlation r as effect size (−1 all samples higher "
               "in the first fraction, +1 all higher in the second), BH-FDR per pair (q < %g). A feature counts as "
               "enriched with q < %g and |r| at least the chosen minimum: with ~180 paired samples almost any shift "
               "is significant. Samples flagged FAIL in a layer's QC are left out." % (q, q), className="lead"),
        html.Div([control("minimum |r|", dd("fr-minr", [{"label": "%g%s" % (x, " (q only)" if x == 0 else ""),
                                                          "value": x} for x in sorted(set(MIN_R + [r0]))], r0))],
                 className="controls"),
        card("Features enriched per fraction", graph("fr-summary"),
             note="share of each pair's tested features enriched in the first (left) or the second (right) fraction"),
        html.Div([control("pair of layers", dd("fr-pair", ps, ps[0]["value"] if ps else None))],
                 className="controls"),
        html.Div([card("Volcano", graph("fr-volcano"),
                       note="each point a feature, its size the median abundance in the fraction where it is higher "
                            "(1e-6 to 10 %%); dotted lines q = %g and the minimum |r|" % q),
                  card("Strongest differences", graph("fr-top"),
                       note="the %d features with the largest |r| on each side (enriched)" % TOP)],
                 className="grid g2"),
        card("Features of the pair", html.Div([
            html.Div([control("show", dd("fr-show", [{"label": "enriched features (q and minimum |r|)", "value": "enriched"},
                                                     {"label": "all tested features", "value": "all"}], "enriched"))],
                     className="controls"),
            html.Div(id="fr-table")]),
             note="prevalence and median relative abundance per fraction, median log2 ratio (second / first), r, "
                  "p, q and where the feature is enriched; sorted by q"),
    ])


def register(app, data):
    q = float(data.q("SELECT max(sig_q) AS q FROM cohort_fraction_summary")["q"].iloc[0] or 0.1) \
        if data.has("cohort_fraction_summary") else 0.1

    @app.callback(Output("fr-summary", "figure"), Output("fr-volcano", "figure"), Output("fr-top", "figure"),
                  Output("fr-table", "children"), Input("fr-pair", "value"), Input("fr-minr", "value"),
                  Input("fr-show", "value"))
    def _pair(pair, min_r, show):
        min_r = float(min_r or 0)
        summ = summary_fig(data, q, min_r)
        if not pair:
            return summ, empty_fig(), empty_fig(), html.P("no pair", className="hint")
        a, b = pair.split("|")
        d = load(data, a, b)
        la, lb = short(data, a), short(data, b)
        strong = (d["q"] < q) & (d["rbc"].abs() >= min_r)
        d = d.assign(enriched=np.where(~strong, "", np.where(d["rbc"] > 0, "second", "first")))
        t = d if show == "all" else d[strong]
        t = t[["feature", "enriched", "n", "prev_a", "prev_b", "median_rel_a", "median_rel_b", "median_log2_ratio",
               "rbc", "p", "q"]].sort_values(["q", "p"]).head(2000)
        t = t.rename(columns={"rbc": "r"})
        pq = t[["p", "q"]].copy()
        t = t.round(5)
        t[["p", "q"]] = pq.apply(lambda c: c.map(lambda x: float("%.3g" % x) if x == x else x))   # 3.2e-27, not 0
        return (summ, volcano(d, q, la, lb, min_r), top_fig(d, q, la, lb, min_r),
                grid("fr-grid", t, height=420, fit=True) if len(t) else
                html.P("no feature with q < %g and |r| ≥ %g" % (q, min_r), className="hint"))


PAGE = Page(path="/fractions", name="Intra vs extra", order=73, layout=layout, register=register,
            requires=("cohort_fraction_summary",))
