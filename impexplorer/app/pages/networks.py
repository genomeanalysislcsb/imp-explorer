"""Metabolic networks (IMP-metabolic): topology of the community networks by group, MAG coverage, reactions that
differ between groups, pan-reactomes of species / genera / families, hub metabolites."""
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from dash import Input, Output, html
from scipy import stats as st

from ..stats import bh
from ..ui import GLOBAL, Page, card, control, dd, empty_fig, graph, style
from ._common import box_by_group, grid, selection

METRICS = {"reactions": "reactions", "metabolites": "metabolites", "components": "connected components",
           "giant_share": "share in the largest component", "mean_degree": "mean degree",
           "clustering": "clustering coefficient", "mean_path_length": "mean path length",
           "scope_metabolites": "scope: metabolites reachable", "scope_reaction_share": "scope: share of reactions"}
SETS = [{"label": "specific reactions", "value": "specific"}, {"label": "all reactions (incl. broad ECs)",
                                                               "value": "all"}]


def topo(data, rset, level="community"):
    d = data.q("SELECT * FROM met_topology WHERE \"set\"=? AND level=?", [rset, level])
    for c in METRICS:
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors="coerce")
    return d


def tests(M, g):
    """per column: Kruskal-Wallis (2 groups: Mann-Whitney), BH-FDR, median per group"""
    g = g.reindex(M.index).fillna("NA")
    lv = [x for x in sorted(g.unique()) if (g == x).sum() >= 2]
    rows = []
    if len(lv) < 2:
        return pd.DataFrame()
    for c in M.columns:
        xs = [M.loc[g == x, c].dropna().values for x in lv]
        try:
            p = st.mannwhitneyu(*xs).pvalue if len(xs) == 2 else st.kruskal(*xs).pvalue
        except ValueError:
            p = 1.0
        r = {"feature": c, "p": p}
        r.update({"median " + str(x): float(np.median(v)) if len(v) else np.nan for x, v in zip(lv, xs)})
        rows.append(r)
    df = pd.DataFrame(rows)
    df["q"] = bh(df["p"].values)
    return df.sort_values("p")[["feature", "q", "p"] + [c for c in df.columns if c.startswith("median")]]


def layout(data):
    ranks = [x for x in data.distinct("met_pan_summary", "rank")]
    return html.Div([
        html.H2("Metabolic networks"),
        html.P("Metabolic networks reconstructed by IMP-metabolic from the Mantis annotation: genes -> MetaNetX "
               "reactions (through KEGG reactions and ECs), one network per community (all genes of the sample), "
               "per MAG and per taxon (pan-reactome of the MAGs of a species, genus or family across the cohort). "
               "Specific reactions have KEGG-reaction evidence or come from ECs with few reactions; 'all' adds the "
               "substrate variants of broad ECs. Topology is computed without currency metabolites (ATP, NAD, ...); "
               "the scope is what a network can make from a minimal glucose medium.", className="lead"),
        html.Div([control("reactions", dd("nw-set", SETS, "specific")),
                  control("metric", dd("nw-metric", [{"label": v, "value": k} for k, v in METRICS.items()],
                                       "reactions")),
                  control("min. prevalence (reactions)", dd("nw-prev", [0.1, 0.25, 0.5], 0.25))],
                 className="controls"),
        html.Div([card("Community networks by group", graph("nw-box")),
                  card("Share of the community's reactions in its MAGs", graph("nw-cov"),
                       note="reactions of the high- and medium-quality MAGs / reactions of the community")],
                 className="grid g2"),
        html.Div(style={"height": "16px"}),
        card("Network properties: differences between groups", html.Div(id="nw-tests"),
             note="community networks; Kruskal-Wallis (2 groups: Mann-Whitney), BH-FDR"),
        html.Div(style={"height": "16px"}),
        card("Reactions present more often in one group", html.Div(id="nw-rx"),
             note="presence of each reaction group in the community networks; reactions present in at least the "
                  "chosen share and at most 95 % of the samples; Kruskal-Wallis / Mann-Whitney on presence, BH-FDR; "
                  "share of samples with the reaction per group"),
        html.Div(style={"height": "16px"}),
        html.Div([control("taxonomic rank", dd("nw-rank", ranks, "species" if "species" in ranks else
                                               (ranks[0] if ranks else None)))], className="controls"),
        html.Div([card("Pan-reactomes: core share vs number of MAGs", graph("nw-pan"),
                       note="core = reaction in >= 90 % of the taxon's MAGs; with few MAGs the core share is "
                            "inflated (2 MAGs: shared by both), so compare taxa with similar MAG numbers"),
                  card("Pan-reactomes", html.Div(id="nw-pantab"))], className="grid g2"),
        html.Div(style={"height": "16px"}),
        card("Hub metabolites of the community networks", html.Div(id="nw-hubs"),
             note="metabolites with the highest betweenness (top 20 per network): in how many of the selected "
                  "communities, median betweenness and degree"),
    ])


def register(app, data):
    rinfo = data.q("SELECT mnx_id, equation, ec_list FROM met_reactions") if data.has("met_reactions") else None

    @app.callback(Output("nw-box", "figure"), Output("nw-cov", "figure"), Output("nw-tests", "children"),
                  Input("nw-set", "value"), Input("nw-metric", "value"),
                  Input(GLOBAL[0], "value"), Input(GLOBAL[1], "value"), Input(GLOBAL[2], "value"),
                  Input(GLOBAL[3], "value"))
    def _topo(rset, metric, group, fcol, fval, qcm=None):
        s, g, cmap = selection(data, group, fcol, fval, qcm, "mg")
        d = topo(data, rset or "specific")
        d = d[d["id"].astype(str).isin(s["sample_id"].astype(str))]
        if d.empty:
            return empty_fig(), empty_fig(), None
        M = d.set_index("id")[[c for c in METRICS if c in d.columns]]
        box = box_by_group(M[metric], g, cmap, METRICS.get(metric, metric), height=400)
        c = data.q("SELECT sample_id, share_in_mags FROM met_mag_coverage WHERE \"set\"=?", [rset or "specific"])
        c = pd.Series(pd.to_numeric(c["share_in_mags"], errors="coerce").values, index=c["sample_id"].astype(str))
        cov = box_by_group(c[c.index.isin(s["sample_id"].astype(str))], g, cmap, "share in MAGs", height=400)
        if g.nunique() > 1:
            t = tests(M, g)
            t["feature"] = t["feature"].map(lambda x: METRICS.get(x, x))
            tab = grid("nw-tgrid", t.rename(columns={"feature": "property"}), height=330, fit=True)
        else:
            tab = html.P("choose a 'group by' variable with at least two groups", className="hint")
        return box, cov, tab

    @app.callback(Output("nw-rx", "children"), Input("nw-set", "value"), Input("nw-prev", "value"),
                  Input(GLOBAL[0], "value"), Input(GLOBAL[1], "value"), Input(GLOBAL[2], "value"),
                  Input(GLOBAL[3], "value"))
    def _rx(rset, prev, group, fcol, fval, qcm=None):
        s, g, _ = selection(data, group, fcol, fval, qcm, "mg")
        if g.nunique() < 2:
            return html.P("choose a 'group by' variable with at least two groups", className="hint")
        w = "sample_id IN (SELECT UNNEST(?::VARCHAR[]))" + (" AND confidence='specific'" if rset == "specific" else "")
        d = data.q("SELECT sample_id, mnx_id FROM met_community_reactions WHERE " + w,
                   [list(s["sample_id"].astype(str))])
        if d.empty:
            return None
        P = pd.crosstab(d["sample_id"].astype(str), d["mnx_id"]).clip(upper=1).astype(float)
        share = P.mean()
        P = P.loc[:, (share >= float(prev or 0.25)) & (share <= 0.95)]
        if P.empty:
            return html.P("no reaction in the chosen prevalence range", className="hint")
        gg = g.reindex(P.index).fillna("NA")
        lv = [x for x in sorted(gg.unique()) if (gg == x).sum() >= 2]
        X = P.values
        ps = []
        for j in range(X.shape[1]):
            xs = [X[(gg == x).values, j] for x in lv]
            try:
                ps.append(st.mannwhitneyu(*xs).pvalue if len(xs) == 2 else st.kruskal(*xs).pvalue)
            except ValueError:
                ps.append(1.0)
        t = pd.DataFrame({"reaction": P.columns, "p": ps})
        t["q"] = bh(t["p"].values)
        for x in lv:
            t["% " + str(x)] = (100 * P[(gg == x).values].mean()).round(1).values
        t = t.sort_values("p").head(500)
        if rinfo is not None:
            t = t.merge(rinfo.rename(columns={"mnx_id": "reaction", "ec_list": "EC"}), on="reaction", how="left")
        first = ["reaction", "q", "p"] + ["% " + str(x) for x in lv]
        t = t[first + [c for c in t.columns if c not in first]]
        return grid("nw-rxgrid", t, height=420, fit=True)

    @app.callback(Output("nw-pan", "figure"), Output("nw-pantab", "children"), Input("nw-set", "value"),
                  Input("nw-rank", "value"))
    def _pan(rset, rank):
        if not rank:
            return empty_fig(), None
        p = data.q("SELECT * FROM met_pan_summary WHERE \"set\"=? AND rank=?", [rset or "specific", rank])
        if p.empty:
            return empty_fig(), None
        for c in ("n_mags", "pan", "core", "accessory", "rare"):
            p[c] = pd.to_numeric(p[c], errors="coerce")
        p["core_share"] = p["core"] / p["pan"]
        fig = go.Figure(go.Scatter(x=p["n_mags"], y=p["core_share"], mode="markers", text=p["taxon"],
                                   marker=dict(size=7, color="#2a6f97", opacity=0.6),
                                   hovertemplate="%{text}<br>%{x} MAGs<br>core share %{y:.2f}<extra></extra>"))
        fig.update_xaxes(title="MAGs of the taxon", type="log")
        fig.update_yaxes(title="core share of the pan-reactome", range=[0, 1.02])
        tab = p.sort_values("n_mags", ascending=False)[["taxon", "n_mags", "pan", "core", "accessory", "rare",
                                                        "core_share"]]
        return style(fig, 420, legend=False), grid("nw-pgrid", tab, height=420, fit=True)

    @app.callback(Output("nw-hubs", "children"), Input("nw-set", "value"),
                  Input(GLOBAL[1], "value"), Input(GLOBAL[2], "value"), Input(GLOBAL[3], "value"))
    def _hubs(rset, fcol, fval, qcm=None):
        s, _, _ = selection(data, None, fcol, fval, qcm, "mg")
        h = data.q("SELECT * FROM met_hubs WHERE \"set\"=? AND level='community'", [rset or "specific"])
        h = h[h["id"].astype(str).isin(s["sample_id"].astype(str))]
        if h.empty:
            return None
        for c in ("betweenness", "degree"):
            h[c] = pd.to_numeric(h[c], errors="coerce")
        keys = ["metabolite"] + (["name"] if "name" in h.columns else [])
        t = h.groupby(keys, as_index=False).agg(communities=("id", "nunique"), betweenness=("betweenness", "median"),
                                                 degree=("degree", "median"))
        return grid("nw-hgrid", t.sort_values(["communities", "betweenness"], ascending=False).head(200),
                    height=380, fit=True)


PAGE = Page(path="/networks", name="Networks", order=47, layout=layout, register=register,
            requires=("met_topology",))
