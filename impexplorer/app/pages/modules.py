"""KEGG modules: completeness of each module in the community (all genes of a sample, Mantis) and in the MAGs
(which species carry a module, in which samples); species profiles across the cohort as a first pan-genome view."""
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from dash import Input, Output, html
from dash.exceptions import PreventUpdate

from ...model import GTDB_RANKS
from ..links import anchor, md
from ..stats import group_tests
from ..ui import GLOBAL, Page, card, control, dd, empty_fig, graph, style
from ._common import box_by_group, grid, order_samples, selection

SEQ = [[0, "#f4f7fa"], [0.5, "#7fa9c9"], [1, "#1f4f6e"]]          # one hue, light -> dark (completeness 0..1)
THRESH = [{"label": "≥ 50 %", "value": 0.5}, {"label": "≥ 75 %", "value": 0.75}, {"label": "≥ 90 %", "value": 0.9},
          {"label": "100 %", "value": 1.0}]


def names(data):
    if not data.has("module_info"):
        return {}
    d = data.q("SELECT module, name FROM module_info")
    return dict(zip(d["module"], d["name"].fillna("")))


def community(data, samples):
    d = data.q("SELECT sample_id, module, completeness FROM modules WHERE scope='community' AND sample_id IN "
               "(SELECT UNNEST(?::VARCHAR[]))", [list(samples["sample_id"])])
    if d.empty:
        return pd.DataFrame()
    return d.pivot_table(index="sample_id", columns="module", values="completeness", fill_value=0.0)


def mag_rows(data, samples, tiers, module=None):
    w = "m.scope='bin' AND m.sample_id IN (SELECT UNNEST(?::VARCHAR[]))"
    p = [list(samples["sample_id"])]
    if module:
        w += " AND m.module=?"
        p.append(module)
    d = data.q("SELECT m.sample_id, m.bin_id, m.module, m.completeness, m.missing, b.quality, %s FROM modules m "
               "JOIN bins b ON b.bin_id = m.bin_id WHERE %s" % (", ".join("b." + r for r in GTDB_RANKS), w), p)
    return d[d["quality"].isin(tiers)] if tiers else d


def heat(M, xlab, ylab, hover, height, zmax=1.0, ztitle="completeness"):
    fig = go.Figure(go.Heatmap(z=M.values, x=xlab, y=ylab, colorscale=SEQ, zmin=0, zmax=zmax, xgap=1, ygap=1,
                               text=hover, hovertemplate="%{text}<extra></extra>",
                               colorbar=dict(title=ztitle, thickness=10, len=0.6)))
    fig.update_yaxes(autorange="reversed", tickfont=dict(size=10))
    fig.update_xaxes(tickfont=dict(size=10), tickangle=-60)
    return style(fig, height, legend=False)


def layout(data):
    has_mag = bool(data.q("SELECT count(*) AS n FROM modules WHERE scope='bin'")["n"].iloc[0]) \
        if data.has("modules") else False
    return html.Div([
        html.H2("KEGG modules"),
        html.P("Completeness of the KEGG modules: of the community (all genes of the sample, Mantis) and of each MAG "
               "(the same score on the KOs of the MAG's genes). Completeness = share of the KOs present on the best "
               "path through the module's definition. A MAG's module can be incomplete because the MAG is "
               "incomplete: compare with its CheckM2 completeness.", className="lead"),
        html.Div([control("MAG quality", dd("mo-tier", ["high", "medium", "low"], ["high", "medium"], multi=True)),
                  control("module counted as present at", dd("mo-thr", THRESH, 1.0)),
                  control("min. prevalence", dd("mo-prev", [0.1, 0.25, 0.5, 0.9], 0.25))], className="controls"),
        card("Community: most variable modules", graph("mo-heat"),
             note="completeness of the 40 modules varying most between the selected samples; samples ordered by "
                  "group"),
        html.Div(style={"height": "16px"}),
        card("Differences between groups", html.Div(id="mo-tests"),
             note="Kruskal-Wallis (2 groups: Mann-Whitney) on the community completeness of each module present "
                  "(> 0) in at least the chosen share of samples, BH-FDR; medians in %; click a row to plot it"),
        html.Div(style={"height": "16px"}),
        html.Div([control("module", dd("mo-mod", [], None, placeholder="type to search")),
                  control("taxonomic rank", dd("mo-rank", list(GTDB_RANKS[1:]), "species")),
                  html.Div(id="mo-link", className="linkbox")], className="controls"),
        html.Div([card("Community completeness by group", graph("mo-box")),
                  card("Which MAGs carry it", graph("mo-carriers"),
                       note="highest completeness of the module among the MAGs of each taxon in each sample "
                            "(blank: no MAG of the taxon); taxa with the module present in most samples first"
                       if has_mag else "needs MAG modules: imp-explorer modules --definitions FILE")],
                 className="grid g2"),
        html.Div(style={"height": "16px"}),
        card("Module profiles of the taxa (pan-genome view)", graph("mo-species"),
             note="share of the MAGs of each taxon (rows: the 40 taxa with most MAGs, at least 2) in which the "
                  "module is present at the chosen completeness; columns: the 40 modules varying most between "
                  "taxa. Core modules of a taxon are near 1, accessory ones in between."),
    ])


def register(app, data):
    nm = names(data)

    @app.callback(Output("mo-heat", "figure"), Output("mo-tests", "children"), Output("mo-mod", "options"),
                  Output("mo-mod", "value"), Input("mo-prev", "value"),
                  Input(GLOBAL[0], "value"), Input(GLOBAL[1], "value"), Input(GLOBAL[2], "value"),
                  Input(GLOBAL[3], "value"))
    def _main(prev, group, fcol, fval, qcm=None):
        s, g, cmap = selection(data, group, fcol, fval, qcm, "mg")
        M = community(data, s)
        if M.empty:
            return empty_fig(), None, [], None
        top = M.var().sort_values(ascending=False).index[:40]
        top = M[top].mean().sort_values(ascending=False).index
        order = order_samples(M[top], g)
        H = M.loc[order, top].T
        hover = [["%s · %s<br>%s<br>completeness %.2f" % (m, nm.get(m, ""), sid, H.loc[m, sid]) for sid in H.columns]
                 for m in H.index]
        ylab = ["%s %s" % (m, (nm.get(m, "")[:38] + "…") if len(nm.get(m, "")) > 39 else nm.get(m, ""))
                for m in H.index]
        xlab = ["%s · %s" % (g.get(x, ""), x) if g.nunique() > 1 else x for x in H.columns]
        fig = heat(H, xlab, ylab, hover, max(420, 18 * len(H) + 160))
        if g.nunique() > 1:
            t = group_tests(M, g.reindex(M.index).fillna("NA").values, min_prev=float(prev or 0.25))
            t = t.rename(columns={"feature": "module"})
            t.insert(1, "name", [nm.get(x, "") for x in t["module"]])
            t.insert(2, "link", [md("kegg_module", x) for x in t["module"]])
            tests = grid("mo-grid", t.head(1000), height=420, fit=True)
            first = t["module"].iloc[0] if len(t) else top[0]
        else:
            tests, first = html.P("choose a 'group by' variable with at least two groups", className="hint"), top[0]
        opts = [{"label": "%s %s" % (m, nm.get(m, "")), "value": m}
                for m in M.mean().sort_values(ascending=False).index]
        return fig, tests, opts, first

    @app.callback(Output("mo-mod", "value", allow_duplicate=True), Input("mo-grid", "selectedRows"),
                  prevent_initial_call=True)
    def _pick(rows):
        if not rows:
            raise PreventUpdate
        return rows[0]["module"]

    @app.callback(Output("mo-link", "children"), Input("mo-mod", "value"))
    def _link(m):
        return anchor("kegg_module", m) if m else None

    @app.callback(Output("mo-box", "figure"), Output("mo-carriers", "figure"), Input("mo-mod", "value"),
                  Input("mo-rank", "value"), Input("mo-tier", "value"), Input("mo-thr", "value"),
                  Input(GLOBAL[0], "value"), Input(GLOBAL[1], "value"), Input(GLOBAL[2], "value"),
                  Input(GLOBAL[3], "value"))
    def _module(m, rank, tiers, thr, group, fcol, fval, qcm=None):
        if not m:
            return empty_fig("choose a module"), empty_fig("choose a module")
        s, g, cmap = selection(data, group, fcol, fval, qcm, "mg")
        M = community(data, s)
        v = M[m] if m in M.columns else pd.Series(0.0, index=M.index)
        box = box_by_group(v.reindex(s["sample_id"]).fillna(0), g, cmap, "completeness (community)",
                           title="%s %s" % (m, nm.get(m, "")), height=420)
        d = mag_rows(data, s, tiers or [], m)
        if d.empty:
            return box, empty_fig("no MAG with this module")
        rank = rank or "species"
        d["taxon"] = d[rank].replace("", np.nan).fillna(d["genus"].replace("", np.nan).map(
            lambda x: "%s (unnamed %s)" % (x, rank) if isinstance(x, str) else x)).fillna("unclassified")
        best = d.sort_values("completeness").groupby(["taxon", "sample_id"]).tail(1)
        P = best.pivot_table(index="taxon", columns="sample_id", values="completeness")
        thr = float(thr or 1.0)
        P = P.loc[(P >= thr - 1e-9).sum(axis=1).sort_values(ascending=False).index[:40]]
        order = [x for x in order_samples(pd.DataFrame(index=s["sample_id"]).assign(k=0), g) if x in P.columns]
        P = P[order]
        bins = best.set_index(["taxon", "sample_id"])["bin_id"]
        hover = [["%s<br>%s<br>%s" % (t, sid, "%s: %.2f" % (bins.get((t, sid), ""), P.loc[t, sid])
                                           if P.loc[t, sid] == P.loc[t, sid] else "no MAG / module absent")
                  for sid in P.columns] for t in P.index]
        xlab = ["%s · %s" % (g.get(x, ""), x) if g.nunique() > 1 else x for x in P.columns]
        fig = heat(P, xlab, list(P.index), hover, max(420, 18 * len(P) + 160))
        fig.update_layout(title="%s in the MAGs (%s)" % (m, rank))
        return box, fig

    @app.callback(Output("mo-species", "figure"), Input("mo-rank", "value"), Input("mo-tier", "value"),
                  Input("mo-thr", "value"), Input(GLOBAL[1], "value"), Input(GLOBAL[2], "value"),
                  Input(GLOBAL[3], "value"))
    def _species(rank, tiers, thr, fcol, fval, qcm=None):
        s, _, _ = selection(data, None, fcol, fval, qcm, "mg")
        d = mag_rows(data, s, tiers or [])
        if d.empty:
            return empty_fig("needs MAG modules: imp-explorer modules --definitions FILE")
        rank = rank or "species"
        d = d[d[rank].fillna("") != ""]
        nmag = d.groupby(rank)["bin_id"].nunique()
        taxa = nmag[nmag >= 2].sort_values(ascending=False).index[:40]
        if not len(taxa):
            return empty_fig("no %s with at least 2 MAGs" % rank)
        d = d[d[rank].isin(taxa)]
        thr = float(thr or 1.0)
        pres = d[d["completeness"] >= thr - 1e-9].groupby([rank, "module"])["bin_id"].nunique()
        P = (pres.unstack(fill_value=0).reindex(taxa).fillna(0).T / nmag.reindex(taxa)).T
        if P.empty:
            return empty_fig("no module present at this completeness")
        cols = P.var().sort_values(ascending=False).index[:40]
        cols = P[cols].mean().sort_values(ascending=False).index
        P = P[cols]
        hover = [["%s (%d MAGs)<br>%s %s<br>present in %.0f %%" % (t, nmag[t], m, nm.get(m, ""), 100 * P.loc[t, m])
                  for m in P.columns] for t in P.index]
        fig = heat(P, list(P.columns), ["%s (%d)" % (t, nmag[t]) for t in P.index], hover,
                   max(460, 18 * len(P) + 180), ztitle="share of MAGs")
        return fig


PAGE = Page(path="/modules", name="Modules", order=45, layout=layout, register=register, requires=("modules",))
