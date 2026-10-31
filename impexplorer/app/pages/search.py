"""Search: one name (taxon or functional term) across every layer of the release, by group."""
import pandas as pd
from ..links import md
import plotly.graph_objects as go
from dash import Input, Output, html, dcc

from ..ui import GLOBAL, Page, card, empty_fig, graph, style
from ._common import grid, selection

SOURCE = {"kraken2_bracken": "Kraken2", "motus": "mOTUs"}
COLS = ["layer", "kind", "name", "link", "samples", "mean_pct"]


def find(data, text, limit=200):
    like = "%" + text.lower() + "%"
    parts = []
    if data.has("taxonomy"):
        parts.append(data.q(
            "SELECT 'taxonomy' AS layer, source || ' · ' || omic || ' · ' || rank AS kind, taxon AS name, "
            "count(DISTINCT sample_id) AS samples, round(100 * avg(rel_abundance), 4) AS mean_pct FROM taxonomy "
            "WHERE lower(taxon) LIKE ? GROUP BY 1, 2, 3", [like]))
    if data.has("bins"):
        for r in ("phylum", "family", "genus", "species"):
            parts.append(data.q(
                "SELECT 'MAGs' AS layer, 'GTDB · %s' AS kind, %s AS name, count(DISTINCT sample_id) AS samples, "
                "NULL AS mean_pct FROM bins WHERE lower(%s) LIKE ? GROUP BY 1, 2, 3" % (r, r, r), [like]))
    if data.has("functions"):
        parts.append(data.q(
            "SELECT 'functions' AS layer, db || ' · ' || omic AS kind, feature AS name, count(DISTINCT sample_id) AS samples, "
            "NULL AS mean_pct FROM functions WHERE lower(feature) LIKE ? GROUP BY 1, 2, 3", [like]))
    d = pd.concat([p for p in parts if len(p)], ignore_index=True) if any(len(p) for p in parts) else pd.DataFrame()
    if not len(d):
        return d
    d = d.sort_values("samples", ascending=False).head(limit)
    d["link"] = [md(k.split(" · ")[0], n) if L == "functions" else (md("gtdb", n) if L == "MAGs" else "")
                 for L, k, n in zip(d["layer"], d["kind"], d["name"])]
    return d[COLS]


def series_for(data, layer, kind, name, samples):
    """{label: Series(sample_id -> % value)} for every layer where the name occurs (same rank)"""
    ids = list(samples["sample_id"])
    out = {}
    if layer == "taxonomy" or layer == "MAGs":
        rank = kind.split(" · ")[-1]
        d = data.q("SELECT sample_id, source, omic, rel_abundance FROM taxonomy WHERE taxon = ? AND rank = ? AND "
                   "sample_id IN (SELECT UNNEST(?::VARCHAR[]))", [name, rank, ids]) if data.has("taxonomy") else None
        if d is not None:
            for (src, omic), x in d.groupby(["source", "omic"]):
                out["%s %s" % (SOURCE.get(src, src), omic)] = x.set_index("sample_id")["rel_abundance"] * 100
        if data.has("bins") and data.has("bin_abundance") and rank in ("phylum", "family", "genus", "species"):
            a = data.q("SELECT a.sample_id, a.omic, sum(a.share_mapped) AS s FROM bin_abundance a JOIN bins b "
                       "USING (bin_id) WHERE b.%s = ? AND a.sample_id IN (SELECT UNNEST(?::VARCHAR[])) "
                       "GROUP BY 1, 2" % rank, [name, ids])
            for omic, x in a.groupby("omic"):
                out["MAGs %s" % omic] = x.set_index("sample_id")["s"] * 100
    else:
        db = kind.split(" · ")[0]
        tot = data.q("SELECT sample_id, omic, sum(reads) AS t FROM functions WHERE db = ? GROUP BY 1, 2", [db])
        d = data.q("SELECT sample_id, omic, reads FROM functions WHERE db = ? AND feature = ?", [db, name])
        d = d.merge(tot, on=["sample_id", "omic"])
        for omic, x in d.groupby("omic"):
            out["%s %s" % (db, omic)] = (x.set_index("sample_id")["reads"] / x.set_index("sample_id")["t"]) * 100
    return {k: v.reindex(ids).fillna(0) for k, v in out.items()}


def layout(data):
    return html.Div([
        html.H2("Search across layers"),
        html.P("Type part of a taxon name (any rank, Kraken2 / mOTUs / GTDB of the MAGs) or a functional term "
               "(KO, module, pathway, CAZyme, EC). Pick a hit to see it in every layer where it occurs.",
               className="lead"),
        dcc.Input(id="se-text", type="text", placeholder="e.g. Akkermansia, K00873, GH13 ...", debounce=True,
                  style={"width": "420px", "padding": "8px 10px", "border": "1px solid #cbd2d9",
                         "borderRadius": "6px", "marginBottom": "12px"}),
        html.Div([card("Hits", html.P(id="se-msg", children="no search yet", className="hint"),
                       grid("se-grid", pd.DataFrame(columns=COLS), height=480, fit=True), note="click a row"),
                  card("Across layers", graph("se-fig"), note="% of the reads of each layer, by group; "
                       "samples without the feature count as 0")], className="grid g2"),
    ])


def register(app, data):
    @app.callback(Output("se-grid", "rowData"), Output("se-msg", "children"), Input("se-text", "value"))
    def _find(text):
        if not text or len(text.strip()) < 2:
            return [], "type at least 2 characters"
        d = find(data, text.strip())
        if d.empty:
            return [], "nothing found"
        return d.to_dict("records"), "%d hits%s" % (len(d), " (first 200)" if len(d) >= 200 else "")

    @app.callback(Output("se-fig", "figure"), Input("se-grid", "selectedRows"),
                  Input(GLOBAL[0], "value"), Input(GLOBAL[1], "value"), Input(GLOBAL[2], "value"))
    def _show(rows, group, fcol, fval):
        if not rows:
            return empty_fig("pick a hit")
        r = rows[0]
        s, g, cmap = selection(data, group, fcol, fval)
        ser = series_for(data, r["layer"], r["kind"], r["name"], s)
        if not ser:
            return empty_fig()
        fig = go.Figure()
        gl = sorted(g.unique())
        pal = ["#1f2933", "#d1495b", "#2a6f97", "#edae49", "#66a182", "#8f5aa8"]
        for i, (lab, v) in enumerate(ser.items()):
            for lv in gl:
                x = v[g.reindex(v.index).values == lv]
                fig.add_trace(go.Box(y=x.values, x=[lv] * len(x), name=lab, legendgroup=lab, offsetgroup=lab,
                                     showlegend=(lv == gl[0]), boxpoints="all", jitter=0.4, pointpos=0,
                                     marker=dict(size=3, color=pal[i % len(pal)], opacity=0.6),
                                     line=dict(color=pal[i % len(pal)], width=1.3), fillcolor="rgba(0,0,0,0)",
                                     text=list(x.index), hovertemplate="%{text}<br>%{y:.3g} %<extra>" + lab +
                                     "</extra>"))
        fig.update_layout(boxmode="group", title="%s (%s)" % (r["name"], r["kind"]))
        fig.update_yaxes(title="% of reads of the layer")
        return style(fig, 480)


PAGE = Page(path="/search", name="Search", order=90, layout=layout, register=register, requires=())
