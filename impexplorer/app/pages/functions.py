"""Functions: KEGG / CAZy / EC terms in the metagenome (potential) and metatranscriptome (expression), group tests,
one term in detail with its expression ratio."""
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from dash.exceptions import PreventUpdate
from dash import Input, Output, html

from ..links import anchor, md, url
from ..stats import group_tests
from ..ui import GLOBAL, Page, card, control, dd, empty_fig, graph, style
from ._common import box_by_group, grid, matrix, order_samples, selection

DB_LABEL = {"kegg_ko": "KEGG orthologs", "kegg_module": "KEGG modules", "kegg_pathway": "KEGG pathways",
            "cazy": "CAZymes", "enzyme_ec": "EC numbers"}


def load(data, db, omic, samples):
    d = data.q("SELECT sample_id, feature, reads FROM functions WHERE db=? AND omic=? AND sample_id IN "
               "(SELECT UNNEST(?::VARCHAR[]))", [db, omic, list(samples["sample_id"])])
    M = matrix(d, value="reads", feature="feature")
    return M.div(M.sum(axis=1).replace(0, np.nan), axis=0).fillna(0) if not M.empty else M


def layout(data):
    dbs = data.distinct("functions", "db")
    omics = data.distinct("functions", "omic")
    return html.Div([
        html.H2("Functions"),
        html.P("Reads per functional term (featureCounts on the Mantis annotation of the sample's genes): the "
               "metagenome gives the functional potential, the metatranscriptome the expression. Values are "
               "relative within each database and omic.", className="lead"),
        html.Div([control("database", dd("fn-db", [{"label": DB_LABEL.get(d, d), "value": d} for d in dbs],
                                         "kegg_module" if "kegg_module" in dbs else (dbs[0] if dbs else None))),
                  control("omic", dd("fn-omic", omics, "mt" if "mt" in omics else (omics[0] if omics else None))),
                  control("min. prevalence", dd("fn-prev", [0.1, 0.25, 0.5, 0.9], 0.5))], className="controls"),
        html.Div([card("Most variable terms", graph("fn-heat"), note="z-score of log10 relative abundance, top 40 "
                       "terms by variance; samples sorted by group"),
                  card("Differences between groups", html.Div(id="fn-tests"),
                       note="Kruskal-Wallis (2 groups: Mann-Whitney) per term at the chosen prevalence, BH-FDR; "
                            "click a row to plot it")], className="grid g2"),
        html.Div(style={"height": "16px"}),
        html.Div([control("term", dd("fn-term", [], None, placeholder="type to search")),
                  html.Div(id="fn-link", className="linkbox")], className="controls"),
        html.Div([card("Term by group", graph("fn-box")),
                  card("Expression ratio", graph("fn-ratio"), note="log2(share in mt / share in mg) per sample: "
                       "> 0 = the term is expressed more than its gene copy number")], className="grid g2"),
    ])


def register(app, data):
    @app.callback(Output("fn-heat", "figure"), Output("fn-tests", "children"), Output("fn-term", "options"),
                  Output("fn-term", "value"), Input("fn-db", "value"), Input("fn-omic", "value"),
                  Input("fn-prev", "value"),
                  Input(GLOBAL[0], "value"), Input(GLOBAL[1], "value"), Input(GLOBAL[2], "value"))
    def _main(db, omic, prev, group, fcol, fval):
        if not db or not omic:
            return empty_fig(), None, [], None
        s, g, cmap = selection(data, group, fcol, fval)
        M = load(data, db, omic, s)
        if M.empty:
            return empty_fig(), None, [], None
        M = M.loc[:, (M > 0).mean() >= float(prev or 0.5)]
        L = np.log10(M + 1e-6)
        top = L.var().sort_values(ascending=False).index[:40]
        Z = (L[top] - L[top].mean()) / L[top].std().replace(0, 1)
        order = order_samples(M[top], g, None)
        Z = Z.loc[order]
        gg = g.reindex(Z.index).fillna("NA")
        heat = go.Figure(go.Heatmap(z=Z.T.values, x=list(Z.index), y=list(Z.columns), colorscale="RdBu_r", zmid=0,
                                    zmin=-3, zmax=3, colorbar=dict(title="z", thickness=10),
                                    hovertemplate="%{y} · %{x}<br>z %{z:.2f}<extra></extra>"))
        pos = 0
        for lv, n in gg.groupby(gg, sort=False).size().items():
            heat.add_annotation(x=pos + n / 2 - 0.5, y=1.02, yref="paper", text="<b>%s</b>" % lv, showarrow=False,
                                font=dict(size=11, color=cmap.get(lv)))
            if pos:
                heat.add_vline(x=pos - 0.5, line=dict(color="#1f2933", width=1.5))
            pos += n
        heat.update_xaxes(showticklabels=False)
        heat.update_yaxes(tickfont=dict(size=9), autorange="reversed", tickmode="array", tickvals=list(Z.columns),
                         ticktext=['<a href="%s" target="_blank">%s</a>' % (url(db, c), c) if url(db, c) else c
                                   for c in Z.columns])
        style(heat, 560, legend=False)
        if g.nunique() > 1:
            t = group_tests(M, g.reindex(M.index).fillna("NA").values,
                            min_prev=float(prev or 0.5))
            t = t.rename(columns={"feature": "term"})
            t.insert(1, "link", [md(db, x) for x in t["term"]])
            tests = grid("fn-grid", t.head(1000), height=560)
            first = t["term"].iloc[0] if len(t) else M.mean().idxmax()
        else:
            tests, first = html.P("choose a 'group by' variable with at least two groups", className="hint"), \
                M.mean().idxmax()
        order = M.mean().sort_values(ascending=False).index[:3000]
        return heat, tests, [{"label": x, "value": x} for x in order], first

    @app.callback(Output("fn-term", "value", allow_duplicate=True), Input("fn-grid", "selectedRows"),
                  prevent_initial_call=True)
    def _pick(rows):
        if not rows:
            raise PreventUpdate
        return rows[0]["term"]

    @app.callback(Output("fn-link", "children"), Input("fn-term", "value"), Input("fn-db", "value"))
    def _link(term, db):
        return anchor(db, term) if term and db else None

    @app.callback(Output("fn-box", "figure"), Output("fn-ratio", "figure"), Input("fn-term", "value"),
                  Input("fn-db", "value"), Input(GLOBAL[0], "value"), Input(GLOBAL[1], "value"),
                  Input(GLOBAL[2], "value"))
    def _term(term, db, group, fcol, fval):
        if not term or not db:
            return empty_fig("choose a term"), empty_fig("choose a term")
        s, g, cmap = selection(data, group, fcol, fval)
        tot = data.q("SELECT sample_id, omic, sum(reads) AS t FROM functions WHERE db=? GROUP BY 1, 2", [db])
        d = data.q("SELECT sample_id, omic, reads FROM functions WHERE db=? AND feature=?", [db, term])
        d = d.merge(tot, on=["sample_id", "omic"])
        d["share"] = d["reads"] / d["t"]
        p = d.pivot_table(index="sample_id", columns="omic", values="share").reindex(s["sample_id"]).fillna(0)
        fig = go.Figure()
        gl = sorted(g.unique())
        for omic in [o for o in ("mg", "mt") if o in p.columns]:
            for lv in gl:
                v = p.loc[g.reindex(p.index).values == lv, omic] * 100
                fig.add_trace(go.Box(y=v.values, x=[lv] * len(v), name=omic, legendgroup=omic, offsetgroup=omic,
                                     showlegend=(lv == gl[0]), boxpoints="all", jitter=0.4, pointpos=0,
                                     marker=dict(size=4, color=cmap.get(lv), opacity=0.7),
                                     line=dict(color="#1f2933" if omic == "mg" else "#d1495b", width=1.3),
                                     fillcolor="rgba(0,0,0,0)", text=list(v.index),
                                     hovertemplate="%{text}<br>%{y:.3g} %<extra>" + omic + "</extra>"))
        fig.update_layout(boxmode="group", title=term)
        fig.update_yaxes(title="%% of %s reads (mg black, mt red)" % DB_LABEL.get(db, db))
        style(fig, 420)
        if not {"mg", "mt"} <= set(p.columns):
            return fig, empty_fig("needs mg and mt")
        q = p[p["mg"] > 0]
        r = np.log2((q["mt"] + 1e-9) / q["mg"])
        return fig, box_by_group(r, g, cmap, "log2(mt / mg)", title=term, height=420)


PAGE = Page(path="/functions", name="Functions", order=40, layout=layout, register=register, requires=("functions",))
