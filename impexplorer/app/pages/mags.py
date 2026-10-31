"""MAGs: the bins of all samples (CheckM2 quality, GTDB-Tk taxonomy) and how much each taxon is present (mg) and
active (mt) across samples."""
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from dash import Input, Output, html

from ...model import GTDB_RANKS
from ..links import anchor
from ..ui import GLOBAL, Page, card, control, dd, empty_fig, graph, style
from ._common import box_by_group, grid, selection

QCOL = {"high": "#2a6f97", "medium": "#edae49", "low": "#9aa5b1", "unknown": "#9aa5b1"}


def bins_of(data, samples, tiers):
    b = data.q("SELECT * FROM bins WHERE sample_id IN (SELECT UNNEST(?::VARCHAR[]))", [list(samples["sample_id"])])
    return b[b["quality"].isin(tiers)] if tiers else b


def layout(data):
    return html.Div([
        html.H2("MAGs"),
        html.P("Metagenome-assembled genomes of every sample (IMP3 MetaWrap refinement of several binners, CheckM2 "
               "quality, GTDB-Tk taxonomy). Abundance = share of the sample's mapped metagenome reads on the bin's "
               "contigs; activity = the same for the metatranscriptome.", className="lead"),
        html.Div([control("quality", dd("mg-tier", ["high", "medium", "low"], ["high", "medium"], multi=True)),
                  control("taxonomic rank", dd("mg-rank", list(GTDB_RANKS[1:]), "genus"))], className="controls"),
        html.Div([card("Quality", graph("mg-qual"), note="CheckM2 completeness vs contamination, MIMAG tiers"),
                  card("Taxonomy of the MAGs", graph("mg-sun"), note="number of MAGs per phylum and the chosen rank "
                       "(top 40)")], className="grid g2"),
        html.Div(style={"height": "16px"}),
        html.Div([control("taxon (chosen rank)", dd("mg-taxon", [], None, placeholder="type to search")),
                  html.Div(id="mg-link", className="linkbox")], className="controls"),
        html.Div([card("Abundance and activity by group", graph("mg-box"),
                       note="per sample: summed share of mapped reads of all MAGs of this taxon"),
                  card("Activity vs abundance per MAG", graph("mg-ratio"),
                       note="log2(mt share / mg share) of each MAG of this taxon; > 0 = more transcribed than its "
                            "genome copy number")], className="grid g2"),
        html.Div(style={"height": "16px"}),
        card("MAG catalogue", html.Div(id="mg-table"), note="all bins of the selected samples and tiers"),
    ])


def register(app, data):
    @app.callback(Output("mg-qual", "figure"), Output("mg-sun", "figure"), Output("mg-table", "children"),
                  Output("mg-taxon", "options"), Output("mg-taxon", "value"),
                  Input("mg-tier", "value"), Input("mg-rank", "value"),
                  Input(GLOBAL[0], "value"), Input(GLOBAL[1], "value"), Input(GLOBAL[2], "value"))
    def _main(tiers, rank, group, fcol, fval):
        s, g, cmap = selection(data, group, fcol, fval)
        b = bins_of(data, s, tiers or [])
        if b.empty:
            return empty_fig(), empty_fig(), None, [], None
        q = go.Figure()
        for tier in ("high", "medium", "low"):
            x = b[b["quality"] == tier]
            if len(x):
                q.add_trace(go.Scattergl(x=x["completeness"], y=x["contamination"], mode="markers", name=tier,
                                         marker=dict(size=4, color=QCOL[tier], opacity=0.55),
                                         text=x["bin_id"] + "<br>" + x["genus"].fillna(""),
                                         hovertemplate="%{text}<br>%{x:.1f} % complete, %{y:.1f} % contamination"
                                                       "<extra></extra>"))
        q.update_xaxes(title="completeness (%)")
        q.update_yaxes(title="contamination (%)")
        style(q, 420)
        rank = rank or "genus"
        c = b.assign(phylum=b["phylum"].replace("", "unclassified"), t=b[rank].replace("", "unclassified"))
        top = c["t"].value_counts().index[:40]
        c = c[c["t"].isin(top)]
        agg = c.groupby(["phylum", "t"]).size().reset_index(name="n")
        if rank == "phylum":
            sun = go.Figure(go.Bar(x=agg["n"], y=agg["t"], orientation="h", marker_color="#2a6f97"))
            sun.update_yaxes(autorange="reversed")
        else:
            ids = list(agg["phylum"].unique()) + list(agg["phylum"] + "/" + agg["t"])
            labels = list(agg["phylum"].unique()) + list(agg["t"])
            parents = [""] * agg["phylum"].nunique() + list(agg["phylum"])
            values = list(agg.groupby("phylum")["n"].sum().reindex(agg["phylum"].unique())) + list(agg["n"])
            sun = go.Figure(go.Sunburst(ids=ids, labels=labels, parents=parents, values=values,
                                        branchvalues="total", insidetextorientation="radial",
                                        hovertemplate="%{label}: %{value} MAGs<extra></extra>"))
        style(sun, 420, legend=False)
        sun.update_layout(margin=dict(l=10, r=10, t=20, b=10))
        cols = ["bin_id", "sample_id", "quality", "completeness", "contamination", "genome_size", "contig_n50",
                "phylum", "family", "genus", "species", "closest_genome_reference", "closest_genome_ani"]
        tbl = grid("mg-grid", b[[x for x in cols if x in b.columns]], height=420, page_size=100)
        counts = b[rank].replace("", np.nan).dropna().value_counts()
        opts = [{"label": "%s (%d MAGs)" % (k, v), "value": k} for k, v in counts.items()]
        return q, sun, tbl, opts, (counts.index[0] if len(counts) else None)

    @app.callback(Output("mg-link", "children"), Input("mg-taxon", "value"))
    def _link(taxon):
        return anchor("gtdb", taxon) if taxon else None

    @app.callback(Output("mg-box", "figure"), Output("mg-ratio", "figure"), Input("mg-taxon", "value"),
                  Input("mg-rank", "value"), Input("mg-tier", "value"),
                  Input(GLOBAL[0], "value"), Input(GLOBAL[1], "value"), Input(GLOBAL[2], "value"))
    def _taxon(taxon, rank, tiers, group, fcol, fval):
        if not taxon:
            return empty_fig("choose a taxon"), empty_fig("choose a taxon")
        s, g, cmap = selection(data, group, fcol, fval)
        b = bins_of(data, s, tiers or [])
        sel = b[b[rank or "genus"] == taxon]
        if sel.empty or not data.has("bin_abundance"):
            return empty_fig(), empty_fig()
        a = data.q("SELECT sample_id, bin_id, omic, share_mapped FROM bin_abundance WHERE bin_id IN "
                   "(SELECT UNNEST(?::VARCHAR[]))", [list(sel["bin_id"])])
        per = a.groupby(["sample_id", "omic"])["share_mapped"].sum().unstack("omic")
        per = per.reindex(s["sample_id"]).fillna(0)       # samples without such a MAG: 0
        fig = go.Figure()
        gl = sorted(g.unique())
        for omic in ("mg", "mt"):
            if omic not in per.columns:
                continue
            for lv in gl:
                v = per.loc[g.reindex(per.index).values == lv, omic] * 100
                fig.add_trace(go.Box(y=v.values, x=[lv] * len(v), name=omic, legendgroup=omic,
                                     showlegend=(lv == gl[0]), offsetgroup=omic, boxpoints="all", jitter=0.4,
                                     pointpos=0, marker=dict(size=4, color=cmap.get(lv), opacity=0.7),
                                     line=dict(color="#1f2933" if omic == "mg" else "#d1495b", width=1.3),
                                     fillcolor="rgba(0,0,0,0)", text=list(v.index),
                                     hovertemplate="%{text}<br>%{y:.3g} %<extra>" + omic + "</extra>"))
        fig.update_layout(boxmode="group", title="%s · present in %d of %d samples" % (
            taxon, sel["sample_id"].nunique(), len(s)))
        fig.update_yaxes(title="% of mapped reads (mg black, mt red)")
        style(fig, 420)
        p = a.pivot_table(index="bin_id", columns="omic", values="share_mapped")
        if not {"mg", "mt"} <= set(p.columns):
            return fig, empty_fig("needs mg and mt")
        p = p[(p["mg"] > 0)]
        r = np.log2((p["mt"] + 1e-7) / p["mg"])
        smp = pd.Series([i.split(":")[0] for i in r.index], index=r.index)
        r.index = smp.values
        return fig, box_by_group(r, g, cmap, "log2(mt / mg) per MAG", title=taxon, height=420)


PAGE = Page(path="/mags", name="MAGs", order=30, layout=layout, register=register, requires=("bins",))
