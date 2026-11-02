"""sRNA: small-RNA results of sRNA-IMP imported with 'imp-explorer srna' (optional layer): library QC, RNA classes,
composition of any feature set (human miRNAs, taxa of the rRNA reads, Rfam, tRNA, functions on the IMP3 assembly) and
one feature next to the metagenome and metatranscriptome of the same samples."""
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from dash import Input, Output, html
from dash.exceptions import PreventUpdate

from ..links import anchor
from ..ui import GLOBAL, PALETTE, Page, card, colors, control, dd, empty_fig, graph, style
from ._common import (HEAT_COLS, HEAT_ROWS, HEAT_VALUES, box_by_group, clustered_heatmap, matrix, selection,
                      stacked)
from .functions import annot_options, annot_series, default_annot

SET_LABEL = {"rna_class": "RNA classes", "mirna": "human miRNAs", "genus": "genera (rRNA reads)",
             "family": "families (rRNA reads)", "species": "species (rRNA reads)", "rfam_type": "Rfam types",
             "rfam_family": "Rfam families", "trna_aa": "tRNA amino acids", "fn_kegg_ko": "KEGG KO (reads on CDS)",
             "fn_kegg_module": "KEGG modules (reads on CDS)", "fn_kegg_pathway": "KEGG pathways (reads on CDS)",
             "fn_pfam": "Pfam (reads on CDS)", "fn_cog": "COG (reads on CDS)"}
SET_ORDER = list(SET_LABEL)
# sRNA feature set -> (IMP3 table, filters) for the comparison with the metagenome / metatranscriptome
IMP3_TWIN = {"genus": ("taxonomy", "genus"), "family": ("taxonomy", "family"), "species": ("taxonomy", "species"),
             "fn_kegg_module": ("functions", "kegg_module"), "fn_kegg_pathway": ("functions", "kegg_pathway"),
             "fn_kegg_ko": ("functions", "kegg_ko")}


def fractions(data):
    fr = data.manifest.get("srna", {}).get("fractions", {})
    have = data.distinct("srna_features", "omic")
    return [{"label": fr.get(f, f), "value": f} for f in have]


def sets(data, frac):
    have = data.distinct("srna_features", '"set"', "omic = ?", [frac]) if frac else []
    have = sorted(have, key=lambda s: SET_ORDER.index(s) if s in SET_ORDER else 99)
    return [{"label": SET_LABEL.get(s, s), "value": s} for s in have]


def groups(data, group, fcol, fval, frac):
    """samples of this sRNA fraction with their group: the global 'group by' if the release has it, else the group
    of the sRNA cohort analysis"""
    s, g, cmap = selection(data, group, fcol, fval)
    ids = set(data.q("SELECT DISTINCT sample_id FROM srna_samples WHERE omic=?", [frac])["sample_id"].astype(str))
    if g.nunique() <= 1:
        sg = data.q("SELECT sample_id, \"group\" FROM srna_samples WHERE omic=?", [frac])
        g = pd.Series(sg["group"].astype(str).values, index=sg["sample_id"].astype(str))
        cmap = colors(g.unique())
    keep = [x for x in s["sample_id"].astype(str) if x in ids] if len(s) else sorted(ids)
    if not keep:                                  # sRNA samples not in the release: use them anyway
        keep = sorted(ids)
    g = g.reindex(keep).fillna("NA")
    return keep, g, cmap


def load(data, frac, st, ids):
    d = data.q("SELECT sample_id, feature, rel FROM srna_features WHERE omic=? AND \"set\"=? "
               "AND sample_id IN (SELECT UNNEST(?::VARCHAR[]))", [frac, st, list(ids)])
    return matrix(d, value="rel", feature="feature")


def layout(data):
    fr = fractions(data)
    f0 = fr[0]["value"] if fr else None
    st = sets(data, f0)
    s0 = next((x["value"] for x in st if x["value"] == "genus"), st[0]["value"] if st else None)
    mets = data.distinct("srna_qc", "metric")
    m0 = next((m for m in ("human_fraction", "nonhuman_rRNA_fraction", "trimmed_reads") if m in mets),
              mets[0] if mets else None)
    return html.Div([
        html.H2("Small RNA"),
        html.P("Results of the sRNA-IMP pipeline for the same samples, imported from its cohort analysis: every "
               "library of each fraction (e.g. extracellular and intracellular sRNA), its QC, the RNA classes of "
               "the non-human reads, human miRNAs, taxa of the rRNA reads (Kraken2/Bracken), Rfam families, tRNAs "
               "and the reads on the genes of the IMP3 assembly. Values are shares within the library and set.",
               className="lead"),
        html.Div([control("fraction", dd("sr-frac", fr, f0)),
                  control("feature set", dd("sr-set", st, s0)),
                  control("features shown", dd("sr-top", [10, 15, 20, 30, 40], 15)),
                  control("view", dd("sr-view", [{"label": "stacked bars", "value": "bars"},
                                                 {"label": "heatmap", "value": "heat"}], "bars"))],
                 className="controls"),
        card("Composition",
             html.Div([control("rows", dd("sr-hrows", HEAT_ROWS, "cluster")),
                       control("samples", dd("sr-hcols", HEAT_COLS, "within")),
                       control("values", dd("sr-hval", HEAT_VALUES, "clr")),
                       control("colour bar", dd("sr-hann", annot_options(data), default_annot(data)))],
                      id="sr-hctl", className="controls", style={"display": "none"}),
             graph("sr-comp"), note="libraries of the fraction, sorted by group; groups from the 'group by' control, "
             "or from the sRNA cohort analysis when the release has no group column"),
        html.Div(style={"height": "16px"}),
        html.Div([card("Library QC", html.Div([dd("sr-metric", mets, m0), graph("sr-qc")])),
                  card("RNA classes by group", graph("sr-classes"),
                       note="mean share of the non-human reads per RNA class")], className="grid g2"),
        html.Div(style={"height": "16px"}),
        html.Div([control("feature", dd("sr-feat", [], None, placeholder="type to search")),
                  html.Div(id="sr-link", className="linkbox")], className="controls"),
        html.Div([card("Feature by group", graph("sr-box")),
                  card("sRNA vs metagenome and metatranscriptome", graph("sr-twin"),
                       note="the same taxon or function in the IMP3 profiles of the same samples (Kraken2 for taxa, "
                            "featureCounts for KEGG); log scale, 0 at the axis minimum")], className="grid g2"),
    ])


def register(app, data):
    @app.callback(Output("sr-set", "options"), Output("sr-set", "value"), Input("sr-frac", "value"),
                  Input("sr-set", "value"))
    def _sets(frac, cur):
        opts = sets(data, frac)
        vals = [o["value"] for o in opts]
        return opts, cur if cur in vals else (next((v for v in vals if v == "genus"), vals[0] if vals else None))

    @app.callback(Output("sr-hctl", "style"), Input("sr-view", "value"))
    def _hctl(view):
        return {"display": "flex"} if view == "heat" else {"display": "none"}

    @app.callback(Output("sr-comp", "figure"), Output("sr-feat", "options"), Output("sr-feat", "value"),
                  Input("sr-frac", "value"), Input("sr-set", "value"), Input("sr-top", "value"),
                  Input("sr-view", "value"), Input("sr-hrows", "value"), Input("sr-hcols", "value"),
                  Input("sr-hval", "value"), Input("sr-hann", "value"),
                  Input(GLOBAL[0], "value"), Input(GLOBAL[1], "value"), Input(GLOBAL[2], "value"))
    def _comp(frac, st, top, view, hrows, hcols, hval, hann, group, fcol, fval):
        if not frac or not st:
            return empty_fig(), [], None
        ids, g, cmap = groups(data, group, fcol, fval, frac)
        M = load(data, frac, st, ids)
        if M.empty:
            return empty_fig(), [], None
        if view == "heat":
            P = M.loc[:, (M > 0).mean() >= 0.1]
            s = data.subset(fcol, fval)
            fig = clustered_heatmap(P, g, cmap, rows=hrows, cols=hcols, values=hval, top=int(top or 15),
                                    annot=annot_series(data, hann, s) if hann and hann != group else None)
        else:
            fig = stacked(M, g, top=int(top or 15), ytitle="share of the %s reads" % SET_LABEL.get(st, st))
        order = list(M.mean().sort_values(ascending=False).index[:3000])
        return fig, [{"label": x, "value": x} for x in order], (order[0] if order else None)

    @app.callback(Output("sr-qc", "figure"), Input("sr-frac", "value"), Input("sr-metric", "value"),
                  Input(GLOBAL[0], "value"), Input(GLOBAL[1], "value"), Input(GLOBAL[2], "value"))
    def _qc(frac, metric, group, fcol, fval):
        if not frac or not metric:
            return empty_fig()
        ids, g, cmap = groups(data, group, fcol, fval, frac)
        d = data.q("SELECT sample_id, value FROM srna_qc WHERE omic=? AND metric=?", [frac, metric])
        v = d.set_index("sample_id")["value"].reindex(ids)
        return box_by_group(v, g, cmap, metric, title="%s · %s" % (metric, frac), height=380)

    @app.callback(Output("sr-classes", "figure"), Input("sr-frac", "value"),
                  Input(GLOBAL[0], "value"), Input(GLOBAL[1], "value"), Input(GLOBAL[2], "value"))
    def _classes(frac, group, fcol, fval):
        if not frac:
            return empty_fig()
        ids, g, cmap = groups(data, group, fcol, fval, frac)
        M = load(data, frac, "rna_class", ids)
        if M.empty:
            return empty_fig("no RNA classes for this fraction")
        mean = M.groupby(g.reindex(M.index).values).mean()
        cols = list(M.mean().sort_values(ascending=False).index)
        fig = go.Figure()
        pal = ["#8c564b", "#1f77b4", "#e377c2", "#2ca02c", "#98df8a", "#ff7f0e", "#d62728", "#ff9896", "#c49c94",
               "#9467bd", "#c5b0d5", "#17becf", "#9edae5", "#bcbd22", "#dbdb8d", "#393b79", "#6b6ecf", "#9c9ede",
               "#637939", "#7f7f7f", "#dddddd"]
        for i, c in enumerate(cols):
            fig.add_trace(go.Bar(x=list(mean.index), y=mean[c].values, name=str(c), marker_color=pal[i % len(pal)],
                                 hovertemplate="%{x}<br>" + str(c) + ": %{y:.1%}<extra></extra>"))
        fig.update_layout(barmode="stack")
        fig.update_yaxes(tickformat=".0%", title="mean share of non-human reads")
        style(fig, 380)
        fig.update_layout(legend=dict(orientation="v", x=1.01, y=1, font=dict(size=10)), margin=dict(r=160))
        return fig

    @app.callback(Output("sr-link", "children"), Input("sr-feat", "value"), Input("sr-set", "value"))
    def _link(feat, st):
        if not feat or not st:
            return None
        db = "gtdb" if st in ("genus", "family", "species") else (st[3:] if st.startswith("fn_kegg") else None)
        return anchor(db, feat) if db else None

    @app.callback(Output("sr-box", "figure"), Output("sr-twin", "figure"), Input("sr-feat", "value"),
                  Input("sr-frac", "value"), Input("sr-set", "value"),
                  Input(GLOBAL[0], "value"), Input(GLOBAL[1], "value"), Input(GLOBAL[2], "value"))
    def _feature(feat, frac, st, group, fcol, fval):
        if not feat or not frac or not st:
            return empty_fig("choose a feature"), empty_fig("choose a feature")
        ids, g, cmap = groups(data, group, fcol, fval, frac)
        d = data.q("SELECT sample_id, rel FROM srna_features WHERE omic=? AND \"set\"=? AND feature=?",
                   [frac, st, feat])
        v = d.set_index("sample_id")["rel"].reindex(ids).fillna(0) * 100
        pos = v[v > 0]
        log = len(pos) > 3 and (pos.max() / max(pos.min(), 1e-12) > 30 or (v == 0).mean() > 0.2)
        box = box_by_group(v, g, cmap, "% of the set" + (" (0 at the bottom)" if log else ""),
                           title="%s · %s · detected in %d of %d" % (feat, frac, int((v > 0).sum()), len(v)), log=log,
                           pseudo=(pos.min() / 2 if log else None), height=400)
        twin = IMP3_TWIN.get(st)
        if not twin or not data.has(twin[0]):
            return box, empty_fig("no IMP3 counterpart for this feature set")
        if twin[0] == "taxonomy":
            t = data.q("SELECT sample_id, omic, rel_abundance AS v FROM taxonomy WHERE source='kraken2_bracken' "
                       "AND rank=? AND taxon=?", [twin[1], feat])
        else:
            t = data.q("SELECT f.sample_id, f.omic, f.reads / t.tot AS v FROM functions f JOIN (SELECT sample_id, omic, "
                       "sum(reads) AS tot FROM functions WHERE db=? GROUP BY 1, 2) t USING (sample_id, omic) "
                       "WHERE f.db=? AND f.feature=?", [twin[1], twin[1], feat])
        if t.empty:
            return box, empty_fig("%s not in the IMP3 profiles" % feat)
        p = t.pivot_table(index="sample_id", columns="omic", values="v", aggfunc="sum").reindex(ids).fillna(0) * 100
        fig = go.Figure()
        gg = g.reindex(p.index).fillna("NA")
        allv = pd.concat([v] + [p[c] for c in p.columns])
        lo = max(allv[allv > 0].min() / 2, 1e-6) if (allv > 0).any() else 1e-6
        for j, omic in enumerate([c for c in ("mg", "mt") if c in p.columns]):
            for lv in sorted(gg.unique()):
                idx = gg.index[gg == lv]
                fig.add_trace(go.Scatter(x=p.loc[idx, omic].clip(lower=lo), y=v.reindex(idx).clip(lower=lo),
                                         mode="markers", name="%s · %s" % (omic, lv), text=list(idx),
                                         marker=dict(size=7, color=cmap.get(lv, "#555"),
                                                     symbol="circle" if omic == "mg" else "triangle-up",
                                                     line=dict(width=0.5, color="#fff")),
                                         hovertemplate="%{text}<br>" + omic + " %{x:.3g} %<br>sRNA %{y:.3g} %"
                                                       "<extra></extra>"))
        rho = {o: (np.corrcoef(p[o].rank(), v.reindex(p.index).rank())[0, 1] if len(p) > 3 else np.nan)
               for o in ("mg", "mt") if o in p.columns}
        fig.update_xaxes(type="log", title="% in the IMP3 profile (circle mg, triangle mt)")
        fig.update_yaxes(type="log", title="% in the sRNA library")
        fig.update_layout(title="%s · Spearman ρ %s" % (feat, ", ".join("%s %.2f" % (k, r) for k, r in rho.items())))
        return box, style(fig, 400)


PAGE = Page(path="/srna", name="sRNA", order=45, layout=layout, register=register, requires=("srna_features",))
