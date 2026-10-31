"""Taxonomy: read-based profiles (Kraken2/Bracken, mOTUs) per rank and omic; composition, ordination, group tests,
one taxon in detail (metagenome vs metatranscriptome)."""
import numpy as np
import plotly.graph_objects as go
from dash.exceptions import PreventUpdate
from dash import Input, Output, State, html

from ..stats import group_tests
from ..links import anchor
from ..ui import GLOBAL, Page, card, control, dd, empty_fig, graph, style
from ._common import box_by_group, grid, matrix, ordination, selection, stacked

SOURCE = {"kraken2_bracken": "Kraken2/Bracken", "motus": "mOTUs"}
RANK_ORDER = ["domain", "phylum", "class", "order", "family", "genus", "species", "motu"]


def profiles(data):
    d = data.q("SELECT DISTINCT source, omic, rank FROM taxonomy")
    d["o"] = d["rank"].map(lambda r: RANK_ORDER.index(r) if r in RANK_ORDER else 99)
    d = d.sort_values(["source", "omic", "o"])
    return [{"label": "%s · %s · %s" % (SOURCE.get(s, s), o, r), "value": "%s|%s|%s" % (s, o, r)}
            for s, o, r in zip(d["source"], d["omic"], d["rank"])]


def load(data, prof, samples):
    s, o, r = prof.split("|")
    d = data.q("SELECT sample_id, taxon, rel_abundance FROM taxonomy WHERE source=? AND omic=? AND rank=? "
               "AND sample_id IN (SELECT UNNEST(?::VARCHAR[]))", [s, o, r, list(samples["sample_id"])])
    return matrix(d)


def layout(data):
    opts = profiles(data)
    v0 = next((o["value"] for o in opts if o["value"] == "kraken2_bracken|mg|genus"), opts[0]["value"] if opts else None)
    return html.Div([
        html.H2("Taxonomy"),
        html.P("Read-based profiles of every sample: Kraken2 + Bracken per rank (database as configured in IMP3) and mOTUs, for the metagenome "
               "(who is there) and the metatranscriptome (who is active). Abundances are relative within each profile.",
               className="lead"),
        html.Div([control("profile", dd("tx-prof", opts, v0)),
                  control("taxa shown", dd("tx-top", [10, 15, 20, 30], 15))], className="controls"),
        card("Composition", graph("tx-bar"), note="top taxa by mean relative abundance; samples sorted by group, "
             "then by the most abundant taxon"),
        html.Div(style={"height": "16px"}),
        html.Div([
            card("Beta diversity", graph("tx-pcoa"), note="Bray-Curtis on relative abundances; PERMANOVA by the "
                 "'group by' variable (499 permutations)"),
            card("Differences between groups", html.Div(id="tx-tests"),
                 note="Kruskal-Wallis (2 groups: Mann-Whitney) per taxon in ≥ 10 % of samples, BH-FDR; click a row "
                      "to plot it"),
        ], className="grid g2"),
        html.Div(style={"height": "16px"}),
        html.Div([control("taxon", dd("tx-taxon", [], None, placeholder="type to search")),
                  control("scale", dd("tx-scale", ["log", "linear"], "log")),
                  html.Div(id="tx-link", className="linkbox")], className="controls"),
        html.Div([card("Taxon by group", graph("tx-box")),
                  card("Present vs active", graph("tx-mgmt"), note="relative abundance in the metagenome vs the "
                       "metatranscriptome of the same profile and rank; above the diagonal = more active than "
                       "abundant")], className="grid g2"),
    ])


def register(app, data):
    @app.callback(Output("tx-bar", "figure"), Output("tx-pcoa", "figure"), Output("tx-tests", "children"),
                  Output("tx-taxon", "options"), Output("tx-taxon", "value"),
                  Input("tx-prof", "value"), Input("tx-top", "value"),
                  Input(GLOBAL[0], "value"), Input(GLOBAL[1], "value"), Input(GLOBAL[2], "value"))
    def _main(prof, top, group, fcol, fval):
        if not prof:
            return empty_fig(), empty_fig(), None, [], None
        s, g, cmap = selection(data, group, fcol, fval)
        M = load(data, prof, s)
        if M.empty:
            return empty_fig(), empty_fig(), None, [], None
        bar = stacked(M, g, top=int(top or 15))
        pc = ordination(M, g, cmap)
        t = group_tests(M, g.reindex(M.index).fillna("NA").values) if g.nunique() > 1 else None
        if t is not None and len(t):
            t = t.rename(columns={"feature": "taxon"})
            tests = grid("tx-grid", t.head(500), height=430)
            first = t["taxon"].iloc[0]
        else:
            tests = html.P("choose a 'group by' variable with at least two groups", className="hint")
            first = M.mean().idxmax()
        order = M.mean().sort_values(ascending=False).index
        return bar, pc, tests, [{"label": x, "value": x} for x in order], first

    @app.callback(Output("tx-taxon", "value", allow_duplicate=True), Input("tx-grid", "selectedRows"),
                  prevent_initial_call=True)
    def _pick(rows):
        if not rows:
            raise PreventUpdate
        return rows[0]["taxon"]

    @app.callback(Output("tx-link", "children"), Input("tx-taxon", "value"), Input("tx-prof", "value"))
    def _link(taxon, prof):
        # GTDB search for Kraken2 names (GTDB or NCBI database, names mostly shared); mOTUs names: no link
        return anchor("gtdb", taxon) if taxon and prof and prof.startswith("kraken2") else None

    @app.callback(Output("tx-box", "figure"), Output("tx-mgmt", "figure"), Input("tx-taxon", "value"),
                  State("tx-prof", "value"), Input("tx-scale", "value"),
                  Input(GLOBAL[0], "value"), Input(GLOBAL[1], "value"), Input(GLOBAL[2], "value"))
    def _taxon(taxon, prof, scale, group, fcol, fval):
        if not taxon or not prof:
            return empty_fig("choose a taxon"), empty_fig("choose a taxon")
        s, g, cmap = selection(data, group, fcol, fval)
        src, omic, rank = prof.split("|")
        d = data.q("SELECT sample_id, omic, rel_abundance FROM taxonomy WHERE source=? AND rank=? AND taxon=? "
                   "AND sample_id IN (SELECT UNNEST(?::VARCHAR[]))", [src, rank, taxon, list(s["sample_id"])])
        cur = d[d["omic"] == omic].set_index("sample_id")["rel_abundance"].reindex(s["sample_id"]).fillna(0)
        lg = scale == "log"
        pc = (cur * 100)
        pseudo = max(pc[pc > 0].min() / 2, 1e-5) if (pc > 0).any() else 1e-5
        box = box_by_group(pc, g, cmap, "%% of %s reads (%s)%s" % (omic, rank, ", 0 at the bottom line" if lg else ""),
                           title="%s · detected in %d of %d samples" % (taxon, int((pc > 0).sum()), len(pc)),
                           log=lg, pseudo=pseudo if lg else None, height=400)
        p = d.pivot_table(index="sample_id", columns="omic", values="rel_abundance", fill_value=0)
        if not {"mg", "mt"} <= set(p.columns):
            return box, empty_fig("needs both mg and mt profiles")
        fig = go.Figure()
        gg = g.reindex(p.index).fillna("NA")
        lo = max(min(p[p > 0].min().min() * 100, 1e-3), 1e-5)
        hi = max(p.max().max() * 100, lo * 10)
        for lv in sorted(gg.unique()):
            q = p[gg.values == lv]
            fig.add_trace(go.Scatter(x=q["mg"] * 100, y=q["mt"] * 100, mode="markers", name=lv, text=list(q.index),
                                     marker=dict(size=8, color=cmap.get(lv), line=dict(width=0.5, color="#fff")),
                                     hovertemplate="%{text}<br>mg %{x:.3g} %<br>mt %{y:.3g} %<extra></extra>"))
        fig.add_trace(go.Scatter(x=[lo, hi], y=[lo, hi], mode="lines", line=dict(color="#9aa5b1", dash="dot"),
                                 showlegend=False, hoverinfo="skip"))
        fig.update_xaxes(type="log", title="% of metagenome reads")
        fig.update_yaxes(type="log", title="% of metatranscriptome reads")
        r = np.corrcoef(p["mg"].rank(), p["mt"].rank())[0, 1] if len(p) > 3 else np.nan
        fig.update_layout(title="%s · Spearman ρ = %.2f" % (taxon, r))
        return box, style(fig, 400)


PAGE = Page(path="/taxonomy", name="Taxonomy", order=20, layout=layout, register=register, requires=("taxonomy",))
