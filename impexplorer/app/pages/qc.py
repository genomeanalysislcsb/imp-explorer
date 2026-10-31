"""QC: any QC metric by group, and two metrics against each other (outliers, depth effects)."""
import numpy as np
import plotly.graph_objects as go
from dash import Input, Output, html

from ..ui import GLOBAL, Page, card, control, dd, graph, style, empty_fig
from ._common import box_by_group, selection

CURATED = [("mg", "raw:pairs", "mg read pairs (raw)"), ("mg", "trimmed:pairs", "mg read pairs (trimmed)"),
           ("mt", "raw:pairs", "mt read pairs (raw)"), ("mt", "trimmed:pairs", "mt read pairs (trimmed)"),
           ("mg", "trimmed:q30_pct", "mg Q30 % (trimmed R1)"), ("mt", "trimmed:q30_pct", "mt Q30 % (trimmed R1)"),
           ("mg", "mapped_fraction", "mg reads mapped to the assembly"),
           ("mt", "mapped_fraction", "mt reads mapped to the assembly"),
           ("mg", "nonpareil_diversity", "Nonpareil diversity (mg)"),
           ("mgmt", "contigs", "assembly: contigs"), ("mgmt", "n50", "assembly: N50"),
           ("mgmt", "assembly_bp", "assembly: total length"), ("mgmt", "imp:total_bins", "bins (IMP3 summary)")]


def metric_options(data):
    have = data.q("SELECT DISTINCT omic, metric FROM qc")
    hs = set(zip(have["omic"], have["metric"]))
    opts = [{"label": lab, "value": "%s|%s" % (o, m)} for o, m, lab in CURATED if (o, m) in hs]
    cur = {o["value"] for o in opts}
    for o, m in sorted(hs):
        v = "%s|%s" % (o, m)
        if v not in cur:
            opts.append({"label": "%s · %s" % (o, m), "value": v})
    return opts


def values(data, key, samples):
    o, m = key.split("|", 1)
    d = data.q("SELECT sample_id, value FROM qc WHERE omic = ? AND metric = ?", [o, m])
    d = d[d["sample_id"].isin(samples["sample_id"])]
    return d.set_index("sample_id")["value"]


def label(opts, key):
    return next((o["label"] for o in opts if o["value"] == key), key)


def layout(data):
    opts = metric_options(data)
    v0 = opts[0]["value"] if opts else None
    return html.Div([
        html.H2("Quality control"),
        html.P("Sequencing depth, read quality, mapping to the sample's own assembly, diversity and assembly "
               "statistics. Curated metrics come first; every IMP3 summary statistic is in the list too.",
               className="lead"),
        html.Div([
            card("Metric by group", html.Div([control("metric", dd("qc-metric", opts, v0)),
                                               control("scale", dd("qc-log", ["linear", "log"], "linear"))],
                                              className="controls"), graph("qc-box")),
            card("Metric vs metric", html.Div([control("x", dd("qc-x", opts, v0)),
                                               control("y", dd("qc-y", opts, opts[6]["value"] if len(opts) > 6 else v0))],
                                              className="controls"), graph("qc-xy"),
                 note="one point per sample; hover for the id"),
        ], className="grid g2"),
    ])


def register(app, data):
    @app.callback(Output("qc-box", "figure"), Input("qc-metric", "value"), Input("qc-log", "value"),
                  Input(GLOBAL[0], "value"), Input(GLOBAL[1], "value"), Input(GLOBAL[2], "value"))
    def _box(key, scale, group, fcol, fval):
        if not key:
            return empty_fig()
        s, g, cmap = selection(data, group, fcol, fval)
        return box_by_group(values(data, key, s), g, cmap, label(metric_options(data), key), log=scale == "log",
                            height=420)

    @app.callback(Output("qc-xy", "figure"), Input("qc-x", "value"), Input("qc-y", "value"),
                  Input(GLOBAL[0], "value"), Input(GLOBAL[1], "value"), Input(GLOBAL[2], "value"))
    def _xy(kx, ky, group, fcol, fval):
        if not kx or not ky:
            return empty_fig()
        s, g, cmap = selection(data, group, fcol, fval)
        opts = metric_options(data)
        x, y = values(data, kx, s), values(data, ky, s)
        idx = x.index.intersection(y.index)
        if not len(idx):
            return empty_fig()
        fig = go.Figure()
        gg = g.reindex(idx).fillna("NA")
        for lv in sorted(gg.unique()):
            i = idx[gg.values == lv]
            fig.add_trace(go.Scatter(x=x[i], y=y[i], mode="markers", name=lv, text=list(i),
                                     marker=dict(size=8, color=cmap.get(lv), opacity=0.8,
                                                 line=dict(width=0.5, color="#fff")),
                                     hovertemplate="%{text}<br>x %{x:.4g}<br>y %{y:.4g}<extra></extra>"))
        if len(idx) > 3:
            r = np.corrcoef(x[idx].rank(), y[idx].rank())[0, 1]
            fig.update_layout(title="Spearman ρ = %.2f (n = %d)" % (r, len(idx)))
        fig.update_xaxes(title=label(opts, kx))
        fig.update_yaxes(title=label(opts, ky))
        return style(fig, 420)


PAGE = Page(path="/qc", name="QC", order=10, layout=layout, register=register, requires=("qc",))
