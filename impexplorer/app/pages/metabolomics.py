"""Metabolomics: assays imported with 'imp-explorer metabolomics' (optional layer): profiles of each assay, one
metabolite by group, and the metabolite next to the microbiome (taxa, functions) of the same samples."""
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from dash import Input, Output, html

from ..ui import GLOBAL, Page, card, control, dd, empty_fig, graph, style
from ._common import HEAT_COLS, HEAT_ROWS, box_by_group, clustered_heatmap, column_defs, grid, selection
from .cohort import metab_units
from .functions import annot_options, annot_series, default_annot

ANNOT = [{"label": "all features", "value": "all"}, {"label": "identified only", "value": "known"}]
WEAK = ("unknown", "name from mass")          # annotation levels left out by 'identified only'
PROFILES = [("kraken.genus.mg", "Kraken2 genera, metagenome"), ("kraken.genus.mt", "Kraken2 genera, metatranscriptome"),
            ("mag.genus.mg", "MAG genera (GTDB), metagenome"), ("kegg_module.mg", "KEGG modules, metagenome"),
            ("kegg_module.mt", "KEGG modules, metatranscriptome"),
            ("prot.genus.mp_extra", "genera, metaproteome extracellular"),
            ("prot.genus.mp_intra", "genera, metaproteome intracellular"),
            ("prot.kegg_module.mp_extra", "KEGG modules, metaproteome extracellular"),
            ("prot.kegg_module.mp_intra", "KEGG modules, metaproteome intracellular")]


def assays(data):
    a = data.q("SELECT * FROM metab_assays") if data.has("metab_assays") else pd.DataFrame(columns=["assay", "label"])
    return [{"label": "%s (%s)" % (r.get("label") or r["assay"], r["assay"]), "value": r["assay"]}
            for r in a.to_dict("records")]


def features(data, assay, annot="all"):
    f = data.q("SELECT * FROM metab_features WHERE assay=?", [assay])
    if annot == "known" and "annotation" in f.columns:
        a = f["annotation"].fillna("").astype(str).str.lower()
        f = f[~a.str.startswith(WEAK)]
    return f


def values(data, assay, ids, feats=None):
    d = data.q("SELECT sample_id, feature_id, value FROM metab_values WHERE assay=? AND "
               "sample_id IN (SELECT UNNEST(?::VARCHAR[]))", [assay, list(ids)])
    if feats is not None:
        d = d[d["feature_id"].isin(set(feats))]
    return d.pivot_table(index="sample_id", columns="feature_id", values="value", aggfunc="first")


def groups(data, group, fcol, fval, assay):
    s, g, cmap = selection(data, group, fcol, fval)
    have = set(data.q("SELECT DISTINCT sample_id FROM metab_values WHERE assay=?", [assay])["sample_id"].astype(str))
    keep = [x for x in s["sample_id"].astype(str) if x in have]
    return keep, g.reindex(keep).fillna("NA"), cmap, s


def profile_matrix(data, key, ids):
    """samples x features (relative abundance) of an IMP3 profile, features in >= 10 % of the samples"""
    p = key.split(".")
    if p[0] == "prot":
        if p[1] == "genus" and data.has("prot_taxonomy"):
            d = data.q("SELECT sample_id, taxon AS f, rel AS v FROM prot_taxonomy WHERE rank=? AND omic=? AND "
                       "sample_id IN (SELECT UNNEST(?::VARCHAR[]))", [p[1], p[2], list(ids)])
        elif data.has("prot_functions"):
            d = data.q("SELECT sample_id, feature AS f, rel AS v FROM prot_functions WHERE db=? AND omic=? AND "
                       "sample_id IN (SELECT UNNEST(?::VARCHAR[]))", [p[1], p[2], list(ids)])
        else:
            return pd.DataFrame()
    elif p[0] == "kraken" and data.has("taxonomy"):
        d = data.q("SELECT sample_id, taxon AS f, rel_abundance AS v FROM taxonomy WHERE source='kraken2_bracken' "
                   "AND rank=? AND omic=? AND sample_id IN (SELECT UNNEST(?::VARCHAR[]))", [p[1], p[2], list(ids)])
    elif p[0] == "mag" and data.has("bin_abundance"):
        d = data.q("SELECT a.sample_id, b.\"%s\" AS f, sum(a.share_mapped) AS v FROM bin_abundance a JOIN bins b "
                   "USING (bin_id) WHERE a.omic=? AND a.sample_id IN (SELECT UNNEST(?::VARCHAR[])) GROUP BY 1, 2"
                   % p[1], [p[2], list(ids)])
        d = d[d["f"].fillna("") != ""]
    elif data.has("functions"):
        d = data.q("SELECT f.sample_id, f.feature AS f, f.reads / t.tot AS v FROM functions f JOIN (SELECT sample_id, "
                   "omic, sum(reads) AS tot FROM functions WHERE db=? GROUP BY 1, 2) t USING (sample_id, omic) "
                   "WHERE f.db=? AND f.omic=? AND f.sample_id IN (SELECT UNNEST(?::VARCHAR[]))",
                   [p[0], p[0], p[1], list(ids)])
    else:
        return pd.DataFrame()
    if d.empty:
        return pd.DataFrame()
    M = d.pivot_table(index="sample_id", columns="f", values="v", aggfunc="sum", fill_value=0.0)
    return M.loc[:, (M > 0).mean() >= 0.1]


def spearman_all(v, M):
    """Spearman rho of v (Series) with every column of M (same index)"""
    rv = v.rank().values
    R = M.rank().values
    rv = (rv - rv.mean()) / (rv.std() or 1)
    sd = R.std(axis=0)
    sd[sd == 0] = np.nan
    Rz = (R - R.mean(axis=0)) / sd
    return pd.Series((Rz * rv[:, None]).mean(axis=0), index=M.columns)


def layout(data):
    opts = assays(data)
    a0 = opts[0]["value"] if opts else None
    prof = [{"label": lab, "value": k} for k, lab in PROFILES
            if (k.startswith("kraken") and data.has("taxonomy")) or (k.startswith("mag") and data.has("bins"))
            or (k.startswith("kegg") and data.has("functions"))
            or (k.startswith("prot.") and data.has("prot_samples") and k.split(".")[2] in data.distinct("prot_samples", "omic"))]
    man = data.manifest.get("metabolomics", {}).get("assays", {})
    rows = [html.Li("%s: %s · %s · %s, acquired %s; %s%s" % (
        a, m.get("label") or a, m.get("platform") or "", m.get("matrix") or "", m.get("acquisition_date") or "?",
        m.get("value_type") or "values", (" (%s)" % m["units"]) if m.get("units") else ""))
        for a, m in man.items()]
    return html.Div([
        html.H2("Metabolomics"),
        html.P("Metabolomics measurements of the same samples, imported with 'imp-explorer metabolomics': one assay "
               "per platform or panel. Group comparisons and adjusted models of the assays are on the Cohort page "
               "(layers 'Metabolomics: …') once 'imp-explorer cohort' has been run after the import.",
               className="lead"),
        html.Ul(rows, className="hint") if rows else None,
        html.Div([control("assay", dd("mb-assay", opts, a0)),
                  control("features", dd("mb-annot", ANNOT, "all")),
                  control("features shown", dd("mb-top", [15, 25, 40, 60, 100], 40)),
                  control("rows", dd("mb-hrows", HEAT_ROWS, "cluster")),
                  control("samples", dd("mb-hcols", HEAT_COLS, "within")),
                  control("colour bar", dd("mb-hann", annot_options(data), default_annot(data)))],
                 className="controls"),
        card("Profiles", graph("mb-heat"), note="z-scores of log10 values of the most variable features; samples "
             "sorted by group ('group by' in the top bar)"),
        html.Div(style={"height": "16px"}),
        html.Div([control("feature", dd("mb-feat", [], None, placeholder="type to search")),
                  control("compare with", dd("mb-prof", prof, prof[0]["value"] if prof else None))],
                 className="controls"),
        html.Div([card("Feature by group", graph("mb-box")),
                  card("Feature and the microbiome", graph("mb-cor"),
                       note="Spearman correlation with every taxon / function of the IMP3 profile present in "
                            ">= 10 % of the samples; the 20 strongest. Exploratory, not adjusted for anything")],
                 className="grid g2"),
        html.Div(style={"height": "16px"}),
        card("Features", grid("mb-grid", pd.DataFrame(), height=380, fit=True),
             note="all features of the assay: annotation, detection and median per group"),
    ])


def register(app, data):
    @app.callback(Output("mb-heat", "figure"), Output("mb-feat", "options"), Output("mb-feat", "value"),
                  Output("mb-grid", "rowData"), Output("mb-grid", "columnDefs"),
                  Input("mb-assay", "value"), Input("mb-annot", "value"), Input("mb-top", "value"),
                  Input("mb-hrows", "value"), Input("mb-hcols", "value"), Input("mb-hann", "value"),
                  Input(GLOBAL[0], "value"), Input(GLOBAL[1], "value"), Input(GLOBAL[2], "value"))
    def _assay(assay, annot, top, hrows, hcols, hann, group, fcol, fval):
        if not assay:
            return empty_fig(), [], None, [], []
        ids, g, cmap, s = groups(data, group, fcol, fval, assay)
        f = features(data, assay, annot)
        M = values(data, assay, ids, f["feature_id"])
        if M.empty:
            return empty_fig("no samples of this assay in the selection"), [], None, [], []
        names = dict(zip(f["feature_id"], f["name"].astype(str)))
        P = M.fillna(0)
        P = P.loc[:, ((P > 0).mean() >= 0.1) & (P.std() > 0)]
        fig = clustered_heatmap(P, g, cmap, rows=hrows, cols=hcols, values="log", top=int(top or 40),
                                annot=annot_series(data, hann, s) if hann and hann != group else None,
                                ylabels=names) if P.shape[1] >= 2 else empty_fig("fewer than 2 features to show")
        # table: annotation, detection, median per group
        t = f.drop(columns=[c for c in ("assay", "source_name") if c in f.columns]).set_index("feature_id")
        t = t.reindex(M.columns)
        t["detected_pct"] = ((M > 0).mean() * 100).round(1)
        for lv in sorted(g.unique()):
            t["median " + lv] = M.loc[g.index[g == lv]].median()
        t = t.reset_index()
        keep = [c for c in ("feature_id", "name", "annotation", "class", "detected_pct") if c in t.columns]
        keep += [c for c in t.columns if c.startswith("median ")]
        keep += [c for c in ("formula", "mz", "rt_min", "ri", "internal_standard", "lab_notes") if c in t.columns]
        t = t[keep]
        t = t[[c for c in t.columns if c in ("feature_id", "name") or t[c].notna().any() and
               (t[c].astype(str).str.strip() != "").any()]]           # columns empty for this assay
        order = list(M.median().sort_values(ascending=False).index)
        opts = [{"label": names.get(x, x) if names.get(x, x) == x else "%s (%s)" % (names.get(x, x), x), "value": x}
                for x in order]
        return fig, opts, order[0] if order else None, t.to_dict("records"), column_defs(t, fit=True)

    @app.callback(Output("mb-box", "figure"), Output("mb-cor", "figure"), Input("mb-feat", "value"),
                  Input("mb-assay", "value"), Input("mb-prof", "value"),
                  Input(GLOBAL[0], "value"), Input(GLOBAL[1], "value"), Input(GLOBAL[2], "value"))
    def _feature(feat, assay, prof, group, fcol, fval):
        if not feat or not assay:
            return empty_fig("choose a feature"), empty_fig("choose a feature")
        ids, g, cmap, _ = groups(data, group, fcol, fval, assay)
        d = data.q("SELECT sample_id, value FROM metab_values WHERE assay=? AND feature_id=?", [assay, feat])
        v = d.set_index("sample_id")["value"].astype(float).reindex(ids).dropna()
        nm = data.q("SELECT name FROM metab_features WHERE assay=? AND feature_id=?", [assay, feat])
        name = nm["name"].iloc[0] if len(nm) else feat
        pos = v[v > 0]
        log = len(pos) > 3 and (v >= 0).all() and (pos.max() / max(pos.min(), 1e-12) > 30 or (v == 0).mean() > 0.2)
        box = box_by_group(v, g, cmap, metab_units(data, assay) + (", 0 at the bottom line" if log else ""),
                           title="%s · detected in %d of %d" % (name, int((v > 0).sum()), len(v)), log=log,
                           pseudo=(pos.min() / 2 if log else None), height=420)
        if not prof:
            return box, empty_fig("no IMP3 profiles in the release")
        M = profile_matrix(data, prof, list(v.index))
        M = M.reindex(v.index).dropna() if not M.empty else M
        if M.empty or len(M) < 8:
            return box, empty_fig("too few samples with both")
        r = spearman_all(v.loc[M.index], M).dropna()
        r = r.reindex(r.abs().sort_values(ascending=False).index[:20]).sort_values()
        fig = go.Figure(go.Bar(x=r.values, y=[str(i) for i in r.index], orientation="h",
                               marker_color=["#2a6f97" if x > 0 else "#d1495b" for x in r.values],
                               hovertemplate="%{y}<br>ρ %{x:.2f}<extra></extra>"))
        fig.update_xaxes(title="Spearman ρ (n = %d samples)" % len(M), range=[-1, 1])
        fig.update_yaxes(tickfont=dict(size=10))
        fig.update_layout(title="%s vs %s" % (name, dict(PROFILES).get(prof, prof)), margin=dict(l=220))
        return box, style(fig, 420, legend=False)


PAGE = Page(path="/metabolomics", name="Metabolomics", order=47, layout=layout, register=register,
            requires=("metab_values",))
