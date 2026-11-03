"""Functions: KEGG / CAZy / EC / Pfam / COG terms in the metagenome (potential) and metatranscriptome (expression), group tests,
one term in detail with its expression ratio."""
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from dash.exceptions import PreventUpdate
from dash import Input, Output, html

from ..links import anchor, md, url
from ..stats import group_tests
from ..ui import GLOBAL, Page, card, control, dd, empty_fig, graph, style
from ._common import (HEAT_COLS, HEAT_ROWS, HEAT_VALUES, OMIC_COLOR, box_by_group, clustered_heatmap, grid,
                      matrix, selection)

DB_LABEL = {"kegg_ko": "KEGG orthologs", "kegg_module": "KEGG modules", "kegg_pathway": "KEGG pathways",
            "cazy": "CAZymes", "enzyme_ec": "EC numbers", "pfam": "Pfam domains", "cog": "COGs"}


def load(data, db, omic, samples):
    d = data.q("SELECT sample_id, feature, reads FROM functions WHERE db=? AND omic=? AND sample_id IN "
               "(SELECT UNNEST(?::VARCHAR[]))", [db, omic, list(samples["sample_id"])])
    M = matrix(d, value="reads", feature="feature")
    return M.div(M.sum(axis=1).replace(0, np.nan), axis=0).fillna(0) if not M.empty else M


def annot_options(data):
    return [{"label": "(none)", "value": ""}] + [{"label": c, "value": c} for c in data.group_cols]


def default_annot(data):
    return next((c for c in ("site", "Site", "centre", "center", "batch") if c in data.group_cols), "")


def annot_series(data, col, samples):
    if not col or col not in samples.columns:
        return None
    return {col: pd.Series(samples[col].astype(str).values, index=samples["sample_id"])}


def layout(data):
    dbs = data.distinct("functions", "db")
    omics = data.distinct("functions", "omic")
    return html.Div([
        html.H2("Functions"),
        html.P("Reads per functional term (featureCounts on the Mantis annotation of the sample's genes): the "
               "metagenome gives the functional potential, the metatranscriptome the expression. Values are "
               "relative within each database and omic. With metaproteomics imported, the omics mp_<fraction> are the "
               "protein abundances (AUC) per KEGG term of the proteins of each fraction.", className="lead"),
        html.Div([control("database", dd("fn-db", [{"label": DB_LABEL.get(d, d), "value": d} for d in dbs],
                                         "kegg_module" if "kegg_module" in dbs else (dbs[0] if dbs else None))),
                  control("omic", dd("fn-omic", omics, "mt" if "mt" in omics else (omics[0] if omics else None))),
                  control("min. prevalence", dd("fn-prev", [0.1, 0.25, 0.5, 0.9], 0.5))], className="controls"),
        card("Most variable terms",
             html.Div([control("terms (rows)", dd("fn-hrows", HEAT_ROWS, "cluster")),
                       control("samples (columns)", dd("fn-hcols", HEAT_COLS, "within")),
                       control("values", dd("fn-hval", HEAT_VALUES, "log")),
                       control("colour bar", dd("fn-hann", annot_options(data), default_annot(data)))],
                      className="controls"),
             graph("fn-heat"),
             note="z-scores of the 40 most variable terms. Clustered: average linkage with optimal leaf order; "
                  "correlation distance between terms, Euclidean distance between samples. CLR: log-ratio to the "
                  "sample's mean term (compositional), then z per term. A tree always forms, also from noise: "
                  "check blocks against the group tests and the ordination."),
        html.Div(style={"height": "16px"}),
        card("Differences between groups", html.Div(id="fn-tests"),
             note="Kruskal-Wallis (2 groups: Mann-Whitney) per term at the chosen prevalence, BH-FDR; "
                  "click a row to plot it"),
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
                  Input("fn-prev", "value"), Input("fn-hrows", "value"), Input("fn-hcols", "value"),
                  Input("fn-hval", "value"), Input("fn-hann", "value"),
                  Input(GLOBAL[0], "value"), Input(GLOBAL[1], "value"), Input(GLOBAL[2], "value"))
    def _main(db, omic, prev, hrows, hcols, hval, hann, group, fcol, fval):
        if not db or not omic:
            return empty_fig(), None, [], None
        s, g, cmap = selection(data, group, fcol, fval)
        M = load(data, db, omic, s)
        if M.empty:
            return empty_fig(), None, [], None
        M = M.loc[:, (M > 0).mean() >= float(prev or 0.5)]
        heat = clustered_heatmap(M, g, cmap, rows=hrows, cols=hcols, values=hval, top=40,
                                 annot=annot_series(data, hann, s) if hann != group else None,
                                 ylabels={c: '<a href="%s" target="_blank">%s</a>' % (url(db, c), c) if url(db, c)
                                          else c for c in M.columns})
        if g.nunique() > 1:
            t = group_tests(M, g.reindex(M.index).fillna("NA").values,
                            min_prev=float(prev or 0.5))
            t = t.rename(columns={"feature": "term"})
            t.insert(1, "link", [md(db, x) for x in t["term"]])
            tests = grid("fn-grid", t.head(1000), height=420, fit=True)
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
        prot = set(data.distinct("prot_samples", "omic")) if data.has("prot_samples") else set()
        for omic in [o for o in ("mg", "mt") if o in p.columns] + sorted(c for c in p.columns if c.startswith("mp_")):
            # proteins: only the samples with that fraction (others are missing, not 0)
            have = set(data.q("SELECT sample_id FROM prot_samples WHERE omic=?", [omic])["sample_id"]) \
                if omic in prot else None
            for lv in gl:
                v = p.loc[g.reindex(p.index).values == lv, omic] * 100
                if have is not None:
                    v = v[v.index.isin(have)]
                fig.add_trace(go.Box(y=v.values, x=[lv] * len(v), name=omic, legendgroup=omic, offsetgroup=omic,
                                     showlegend=(lv == gl[0]), boxpoints="all", jitter=0.4, pointpos=0,
                                     marker=dict(size=4, color=cmap.get(lv), opacity=0.7),
                                     line=dict(color=OMIC_COLOR.get(omic, "#8f5aa8"), width=1.3),
                                     fillcolor="rgba(0,0,0,0)", text=list(v.index),
                                     hovertemplate="%{text}<br>%{y:.3g} %<extra>" + omic + "</extra>"))
        fig.update_layout(boxmode="group", title=term)
        fig.update_yaxes(title="%% of %s (reads; proteins: AUC)" % DB_LABEL.get(db, db))
        style(fig, 420)
        if not {"mg", "mt"} <= set(p.columns):
            return fig, empty_fig("needs mg and mt")
        q = p[p["mg"] > 0]
        r = np.log2((q["mt"] + 1e-9) / q["mg"])
        return fig, box_by_group(r, g, cmap, "log2(mt / mg)", title=term, height=420)


PAGE = Page(path="/functions", name="Functions", order=40, layout=layout, register=register, requires=("functions",))
