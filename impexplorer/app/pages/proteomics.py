"""Proteomics: metaproteomics imported with 'imp-explorer proteomics' (optional layer). The protein profiles per
taxon and KEGG term are also on the Taxonomy and Functions pages (omics mp_<fraction>); this page shows the library QC
and what only proteins give: which taxa make the proteins of a function, which functions a taxon's proteins have, and
the proteins themselves."""
import pandas as pd
import plotly.graph_objects as go
from dash import Input, Output, html

from ..links import anchor
from ..ui import GLOBAL, Page, card, control, dd, empty_fig, graph, style
from ._common import box_by_group, column_defs, grid, omic_label, selection

DBS = [{"label": "KEGG module", "value": "kegg_module"}, {"label": "KEGG ortholog", "value": "kegg_ko"},
       {"label": "KEGG pathway", "value": "kegg_pathway"}, {"label": "Pfam domain", "value": "pfam"},
       {"label": "COG", "value": "cog"}]


def dbs(data):
    have = set(data.distinct("prot_functions", "db"))
    return [d for d in DBS if d["value"] in have]
METRICS = [{"label": "proteins identified", "value": "n_proteins"}, {"label": "peptides", "value": "n_peptides"},
           {"label": "total protein abundance (AUC)", "value": "total_auc"}]
GENUS = "regexp_extract(p.lineage, '(?:^|;)g_([^;]+)', 1)"
PAL = ["#2a6f97", "#d1495b", "#edae49", "#66a182", "#8f5aa8", "#00798c", "#e07a5f", "#3d405b", "#81b29a",
       "#f2cc8f", "#9aa5b1"]


def fractions(data):
    os_ = data.distinct("prot_samples", "omic")
    return [{"label": omic_label(data, o), "value": o} for o in os_]


def groups(data, group, fcol, fval, omic):
    s, g, cmap = selection(data, group, fcol, fval)
    have = set(data.q("SELECT sample_id FROM prot_samples WHERE omic=?", [omic])["sample_id"].astype(str))
    keep = [x for x in s["sample_id"].astype(str) if x in have]
    return keep, g.reindex(keep).fillna("NA"), cmap


def term_genera(data, omic, db, term, ids):
    """samples x genera: share of the sample's AUC of the term's proteins made by each genus"""
    d = data.q("SELECT a.sample_id, coalesce(nullif(%s, ''), 'unclassified') AS genus, sum(a.auc) AS auc "
               "FROM prot_abundance a JOIN (SELECT * FROM prot_proteins WHERE list_contains(string_split(%s, ','), ?)) p "
               "USING (protein_id) WHERE a.omic = ? AND a.sample_id IN (SELECT UNNEST(?::VARCHAR[])) GROUP BY 1, 2"
               % (GENUS, db), [term, omic, list(ids)])
    if d.empty:
        return pd.DataFrame()
    M = d.pivot_table(index="sample_id", columns="genus", values="auc", aggfunc="sum", fill_value=0)
    return M.div(M.sum(axis=1), axis=0)


def genus_terms(data, omic, db, genus, ids):
    """samples x terms: share of the genus's protein AUC per term"""
    d = data.q("SELECT a.sample_id, t.term, sum(a.auc) AS auc FROM prot_abundance a JOIN (SELECT protein_id, "
               "UNNEST(string_split(%s, ',')) AS term FROM prot_proteins p WHERE %s = ? AND %s <> '') t "
               "USING (protein_id) WHERE a.omic = ? AND a.sample_id IN (SELECT UNNEST(?::VARCHAR[])) GROUP BY 1, 2"
               % (db, GENUS, db), [genus, omic, list(ids)])
    # denominator: the genus's proteins with a term of this database (each protein once)
    tot = data.q("SELECT a.sample_id, sum(a.auc) AS t FROM prot_abundance a JOIN prot_proteins p USING (protein_id) "
                 "WHERE %s = ? AND p.%s <> '' AND a.omic = ? GROUP BY 1" % (GENUS, db), [genus, omic])
    if d.empty:
        return pd.DataFrame()
    M = d.pivot_table(index="sample_id", columns="term", values="auc", aggfunc="sum", fill_value=0)
    return M.div(tot.set_index("sample_id")["t"].reindex(M.index), axis=0)


def group_bars(M, g, top, title, xtitle):
    """mean share per group of the top columns of M (+ other), horizontal stacked bars"""
    if M.empty:
        return empty_fig("nothing for this selection")
    gg = g.reindex(M.index).fillna("NA")
    mean = M.groupby(gg.values).mean()
    cols = list(M.mean().sort_values(ascending=False).index[:top])
    other = 1 - mean[cols].sum(axis=1)
    fig = go.Figure()
    for i, c in enumerate(cols + ["other"]):
        y = mean[c] if c != "other" else other.clip(lower=0)
        fig.add_trace(go.Bar(y=["%s (n=%d)" % (lv, (gg == lv).sum()) for lv in mean.index], x=y * 100, name=str(c),
                             orientation="h", marker_color=PAL[i % len(PAL)] if c != "other" else "#d9dde2",
                             hovertemplate="%{y}<br>" + str(c).replace("%", "%%") + ": %{x:.1f} %<extra></extra>"))
    fig.update_xaxes(title=xtitle, range=[0, 100])
    style(fig, 230 + 46 * len(mean))
    fig.update_layout(barmode="stack", title=title, showlegend=True,
                      legend=dict(font=dict(size=10), orientation="h", x=0, xanchor="left", y=-0.32, yanchor="top"),
                      margin=dict(b=120))
    return fig


def term_bars(M, g, cmap, top, title, xtitle):
    """mean share per group of the top terms of M: grouped horizontal bars, one colour per group"""
    if M.empty:
        return empty_fig("nothing for this selection")
    gg = g.reindex(M.index).fillna("NA")
    mean = M.groupby(gg.values).mean()
    cols = list(M.mean().sort_values(ascending=False).index[:top])[::-1]
    fig = go.Figure()
    for lv in mean.index:
        fig.add_trace(go.Bar(y=[str(c) for c in cols], x=mean.loc[lv, cols] * 100, name="%s (n=%d)" % (
            lv, (gg == lv).sum()), orientation="h", marker_color=cmap.get(lv, "#555"),
            hovertemplate="%{y}<br>" + str(lv) + ": %{x:.2f} %<extra></extra>"))
    fig.update_xaxes(title=xtitle)
    style(fig, 160 + 16 * len(cols) * max(len(mean), 1))
    fig.update_layout(barmode="group", title=title, showlegend=True,
                      legend=dict(orientation="h", x=1, xanchor="right", y=1.0, yanchor="bottom"),
                      title_x=0, margin=dict(t=90))
    return fig


def layout(data):
    fr = fractions(data)
    f0 = fr[0]["value"] if fr else None
    full = data.has("prot_proteins") and data.has("prot_abundance")
    gen = data.q("SELECT taxon, avg(rel) AS r FROM prot_taxonomy WHERE rank='genus' GROUP BY 1 ORDER BY 2 DESC "
                 "LIMIT 400")["taxon"].tolist() if data.has("prot_taxonomy") else []
    man = data.manifest.get("proteomics", {})
    return html.Div([
        html.H2("Proteomics"),
        html.P("Metaproteomics of the same samples, imported with 'imp-explorer proteomics': proteins identified "
               "against the protein catalogue of the IMP3 assemblies, quantified by their peak area (AUC), per "
               "fraction. The protein profiles per taxon and KEGG term are on the Taxonomy and Functions pages "
               "next to mg and mt (omics %s), and on the Cohort page as layers 'Metaproteome …' once "
               "'imp-explorer cohort' has run after the import." % ", ".join(o["value"] for o in fr),
               className="lead"),
        html.Div([control("fraction", dd("pr-frac", fr, f0)), control("library metric", dd("pr-metric", METRICS,
                                                                                         "n_proteins"))],
                 className="controls"),
        card("Libraries", graph("pr-qc"), note="per sample and fraction: proteins identified, peptides, total AUC "
             "(%d proteins in the catalogue)" % int(man.get("proteins", 0) or 0)),
        html.Div(style={"height": "16px"}),
        *([html.Div([control("database", dd("pr-db", dbs(data), "kegg_module")),
                     control("term", dd("pr-term", [], None, placeholder="type to search")),
                     control("genus", dd("pr-genus", gen, gen[0] if gen else None))], className="controls"),
           html.Div([card("Who makes it", graph("pr-who"), note="share of the term's protein abundance made by "
                          "each genus (genus of the protein's contig), mean per group"),
                     card("What the genus makes", graph("pr-what"), note="its 15 main terms of the chosen database: "
                          "share of the genus's protein abundance with a term of that database, mean per group "
                          "(a protein with several terms counts for each)")], className="grid g2"),
           html.Div(style={"height": "16px"}),
           card("Proteins of the term", html.Div(id="pr-link", className="linkbox"),
                grid("pr-grid", pd.DataFrame(), height=360, fit=True),
                note="the term's proteins in this fraction: taxon, product, samples detected, mean AUC")]
          if full else [html.P("The protein catalogue and abundances were not imported (--no-proteins): taxon x "
                               "function views are not available.", className="hint")]),
    ])


def register(app, data):
    @app.callback(Output("pr-qc", "figure"), Input("pr-frac", "value"), Input("pr-metric", "value"),
                  Input(GLOBAL[0], "value"), Input(GLOBAL[1], "value"), Input(GLOBAL[2], "value"))
    def _qc(omic, metric, group, fcol, fval):
        if not omic or not metric:
            return empty_fig()
        ids, g, cmap = groups(data, group, fcol, fval, omic)
        d = data.q("SELECT * FROM prot_samples WHERE omic=?", [omic])
        if metric not in d.columns:
            return empty_fig("%s not in the import" % metric)
        v = d.set_index("sample_id")[metric].astype(float).reindex(ids)
        return box_by_group(v, g, cmap, dict((m["value"], m["label"]) for m in METRICS)[metric],
                            title=omic_label(data, omic), height=380)

    if not (data.has("prot_proteins") and data.has("prot_abundance")):
        return

    @app.callback(Output("pr-term", "options"), Output("pr-term", "value"), Input("pr-db", "value"),
                  Input("pr-frac", "value"))
    def _terms(db, omic):
        if not db or not omic:
            return [], None
        d = data.q("SELECT feature, avg(rel) AS r FROM prot_functions WHERE db=? AND omic=? GROUP BY 1 "
                   "ORDER BY 2 DESC LIMIT 3000", [db, omic])
        opts = list(d["feature"])
        return opts, opts[0] if opts else None

    @app.callback(Output("pr-who", "figure"), Output("pr-grid", "rowData"), Output("pr-grid", "columnDefs"),
                  Output("pr-link", "children"), Input("pr-term", "value"), Input("pr-db", "value"),
                  Input("pr-frac", "value"), Input(GLOBAL[0], "value"), Input(GLOBAL[1], "value"),
                  Input(GLOBAL[2], "value"))
    def _who(term, db, omic, group, fcol, fval):
        if not term or not db or not omic:
            return empty_fig("choose a term"), [], [], ""
        ids, g, cmap = groups(data, group, fcol, fval, omic)
        M = term_genera(data, omic, db, term, ids)
        fig = group_bars(M, g, 10, "%s · %s" % (term, omic), "% of the term's protein abundance")
        t = data.q("SELECT p.protein_id, p.taxon, p.rank, p.product, p.%s AS terms, count(DISTINCT a.sample_id) AS "
                   "samples, avg(a.auc) AS mean_auc FROM prot_abundance a JOIN (SELECT * FROM prot_proteins WHERE "
                   "list_contains(string_split(%s, ','), ?)) p USING (protein_id) WHERE a.omic = ? AND a.sample_id IN "
                   "(SELECT UNNEST(?::VARCHAR[])) GROUP BY 1, 2, 3, 4, 5 ORDER BY 6 DESC, 7 DESC LIMIT 500"
                   % (db, db), [term, omic, list(ids)])
        return fig, t.to_dict("records"), column_defs(t, fit=True), anchor(db, term)

    @app.callback(Output("pr-what", "figure"), Input("pr-genus", "value"), Input("pr-db", "value"),
                  Input("pr-frac", "value"), Input(GLOBAL[0], "value"), Input(GLOBAL[1], "value"),
                  Input(GLOBAL[2], "value"))
    def _what(genus, db, omic, group, fcol, fval):
        if not genus or not db or not omic:
            return empty_fig("choose a genus")
        ids, g, cmap = groups(data, group, fcol, fval, omic)
        M = genus_terms(data, omic, db, genus, ids)
        return term_bars(M, g, cmap, 15, "%s · %s" % (genus, omic), "% of the genus's annotated protein abundance")


PAGE = Page(path="/proteomics", name="Proteomics", order=46, layout=layout, register=register,
            requires=("prot_samples",))
