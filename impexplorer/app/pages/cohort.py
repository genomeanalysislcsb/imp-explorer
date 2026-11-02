"""Cohort: results of 'imp-explorer cohort' (IMP3 layers) and, if imported with 'imp-explorer srna', of the sRNA-IMP
cohort analysis (sRNA layers): group comparisons within blocks, adjusted models, beta diversity, QC and covariates,
browsable per layer, comparison, model and term. Shown when the release has cohort results of either."""
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from dash import Input, Output, State, ctx, dcc, html
from dash.exceptions import PreventUpdate

from ..links import md
from ..ui import PALETTE, Page, card, control, dd, empty_fig, graph, style
from ._common import box_by_group, grid


def settings(data):
    """settings of the IMP3 cohort analysis, else of the imported sRNA cohort analysis"""
    return data.manifest.get("cohort") or data.manifest.get("srna", {})


def srna(data):
    return data.manifest.get("srna", {})


def is_srna(layer):
    return str(layer).startswith("srna_")


def layer_label(data, key):
    if key == "metrics":
        return "IMP3 library metrics"
    for st in (data.manifest.get("cohort", {}), srna(data)):
        lab = st.get("layers", {}).get(key, {}).get("label")
        if lab:
            return lab
    return key


def layer_db(key):
    p = key.split(".")
    if is_srna(key):
        s = p[1] if len(p) > 1 else ""
        if s in ("genus", "family", "species"):
            return "gtdb"
        return s[3:] if s.startswith("fn_kegg_") or s == "fn_cog" else None
    if p[0] in ("kraken", "mag", "motus") or (p[0] == "ratio" and p[1] == "mag"):
        return "gtdb" if p[0] != "motus" else None
    return p[1] if p[0] == "ratio" else p[0]


def ordered(data, table, col):
    """distinct values in the order of the cohort settings (layers) or of first appearance (comparisons)"""
    vals = data.distinct(table, col)
    if col == "layer":
        order = ["metrics"] + list(data.manifest.get("cohort", {}).get("layers", {})) + list(srna(data).get("layers", {}))
        return [v for v in order if v in vals] + [v for v in vals if v not in order]
    d = data.q("SELECT DISTINCT %s AS v, kind FROM %s" % (col, table)) if col != "layer" else pd.DataFrame()
    rank = {"within": 0, "pooled": 0, "block": 1, "groups": 2, "partial": 3, "kruskal": 4}
    d["r"] = d["kind"].map(rank).fillna(5)
    return list(d.sort_values(["r", "v"])["v"].drop_duplicates())


def gcolors(data):
    gs = settings(data).get("groups") or srna(data).get("groups") or []
    return {g: PALETTE[i % len(PALETTE)] for i, g in enumerate(gs)}


def cells(data):
    """sample -> 'group · block', colour map of these cells (colour of the group); IMP3 cohort samples, plus the
    sRNA samples"""
    parts = []
    if data.has("cohort_samples"):
        parts.append(data.q("SELECT sample_id, \"group\", block FROM cohort_samples"))
    if data.has("srna_samples"):
        parts.append(data.q("SELECT DISTINCT sample_id, \"group\", block FROM srna_samples"))
    s = pd.concat(parts).drop_duplicates("sample_id") if parts else pd.DataFrame(columns=["sample_id", "group",
                                                                                          "block"])
    st = settings(data)
    lab = s["group"].astype(str) + ((" · " + s["block"].astype(str)) if st.get("block") else "")
    gc = gcolors(data)
    cmap = {l: gc.get(g, "#555") for l, g in zip(lab, s["group"].astype(str))}
    return pd.Series(lab.values, index=s["sample_id"].astype(str)), cmap


def feature_values(data, layer, feature):
    """per-sample value of one feature of a layer, from the release tables; (Series, axis title)"""
    p = layer.split(".")
    if is_srna(layer):
        frac, st = p[0], p[1]
        if st == "metrics":
            d = data.q("SELECT sample_id, value AS v FROM srna_qc WHERE omic=? AND metric=?", [frac, feature])
            return d.set_index("sample_id")["v"].astype(float), "%s (%s)" % (feature, frac)
        d = data.q("SELECT sample_id, rel AS v FROM srna_features WHERE omic=? AND \"set\"=? AND feature=?",
                   [frac, st, feature])
        return d.set_index("sample_id")["v"], "share of the %s reads (%s)" % (
            {"mirna": "human miRNA", "rna_class": "non-human"}.get(st, st.replace("fn_", "")), frac)
    if layer == "metrics":
        d = data.q("SELECT * FROM cohort_samples")
        if feature not in d.columns:
            return pd.Series(dtype=float), feature
        return d.set_index("sample_id")[feature].astype(float), settings(data).get("metrics", {}).get(feature, feature)

    def tax(src, rank, omic):
        d = data.q("SELECT sample_id, rel_abundance AS v FROM taxonomy WHERE source=? AND rank=? AND omic=? AND "
                   "taxon=?", [src, rank, omic, feature])
        return d.set_index("sample_id")["v"]

    def mag(rank, omic):
        d = data.q("SELECT a.sample_id, sum(a.share_mapped) AS v FROM bin_abundance a JOIN bins b USING (bin_id) "
                   "WHERE a.omic=? AND b.\"%s\"=? GROUP BY 1" % rank, [omic, feature])
        return d.set_index("sample_id")["v"]

    def fun(db, omic):
        d = data.q("SELECT sample_id, sum(CASE WHEN feature=? THEN reads ELSE 0 END) / sum(reads) AS v FROM functions "
                   "WHERE db=? AND omic=? GROUP BY 1", [feature, db, omic])
        return d.set_index("sample_id")["v"]

    def get(parts):
        if parts[0] in ("kraken", "motus"):
            return tax("kraken2_bracken" if parts[0] == "kraken" else "motus", parts[1], parts[2])
        if parts[0] == "mag":
            return mag(parts[1], parts[2])
        return fun(parts[0], parts[1])
    if p[0] == "ratio":
        base = p[1:]
        a, b = get(base + ["mg"]), get(base + ["mt"])
        idx = a.index.intersection(b.index)
        with np.errstate(divide="ignore", invalid="ignore"):
            r = np.log2((b.loc[idx] + 1e-7) / (a.loc[idx] + 1e-7))
        return r, "log2(mt share / mg share)"
    v = get(p)
    unit = {"mag": "share of mapped reads on its MAGs"}.get(p[0], "relative abundance")
    return v, "%s (%s)" % (unit, p[-1])


def feature_box(data, layer, feature, subtitle=""):
    if not layer or feature is None:
        return empty_fig("click a row of the table")
    v, ytitle = feature_values(data, layer, feature)
    lab, cmap = cells(data)
    v = v.reindex(lab.index)
    plain = layer.startswith(("metrics", "ratio")) or layer.endswith(".metrics")
    if not plain:
        v = v.fillna(0)
    pos = v[v > 0]
    log = not plain and len(pos) > 3 and (
        pos.max() / max(pos.min(), 1e-12) > 30 or (v == 0).mean() > 0.2)
    pseudo = pos.min() / 2 if log else None
    name = ytitle if layer == "metrics" else feature
    if is_srna(layer):
        v = v.dropna() if plain else v   # samples without this sRNA library are not 0
    return box_by_group(v, lab, cmap, ytitle + (", 0 at the bottom line" if log else ""), log=log, pseudo=pseudo,
                        title="%s · %s%s" % (name, layer_label(data, layer), subtitle), height=420)


def volcano(d, x, xtitle, q, title):
    if d.empty:
        return empty_fig("no results for this selection")
    y = -np.log10(d["p"].clip(lower=1e-300))
    sig = d["q"] < q
    fig = go.Figure()
    for flag, name, col in ((False, "q ≥ %g" % q, "#b8c2cc"), (True, "q < %g" % q, "#d1495b")):
        m = sig == flag
        fig.add_trace(go.Scattergl(x=d.loc[m, x], y=y[m], mode="markers", name=name,
                                   text=(d.loc[m, "feature"].astype(str) + " " + d.loc[m, "name"].fillna("").astype(str)),
                                   customdata=d.loc[m, "q"],
                                   marker=dict(size=7 if flag else 5, color=col, opacity=0.85 if flag else 0.5),
                                   hovertemplate="%{text}<br>" + xtitle + " %{x:.3g}<br>p %{y:.2f} (-log10)"
                                                 "<br>q %{customdata:.3g}<extra></extra>"))
    fig.add_vline(x=0, line=dict(color="#9aa5b1", width=1))
    fig.update_xaxes(title=xtitle)
    fig.update_yaxes(title="-log10 p")
    fig.update_layout(title=title)
    return style(fig, 420)


def table(id_, d, cols, layer):
    cols = [c for c in cols if c in d.columns and d[c].notna().any() and not (d[c].astype(str) == "").all()]
    d = d[cols].copy()
    if "name" in d.columns and not d["name"].fillna("").astype(str).str.strip().replace("nan", "").any():
        d = d.drop(columns=["name"])
    db = layer_db(layer) if layer else None
    if db and len(d) and "feature" in d.columns:
        d.insert(d.columns.get_loc("feature") + 1, "link", [md(db, f) for f in d["feature"]])
    return grid(id_, d.head(2000), height=420, fit=True)


def tab_tests(data):
    lays = ordered(data, "all_tests", "layer")
    comps = ordered(data, "all_tests", "comparison")
    c0 = next((c for c in comps if c != "all groups"), comps[0] if comps else None)
    return html.Div([
        html.Div([control("layer", dd("co-t-layer", [{"label": layer_label(data, k), "value": k} for k in lays],
                                      lays[0] if lays else None)),
                  control("comparison", dd("co-t-comp", comps, c0))], className="controls"),
        html.Div([card("Volcano", graph("co-t-volcano"), note="effect = P(first group > second) − 0.5 (van Elteren "
                       "within blocks / Mann-Whitney); all groups: Kruskal-Wallis H"),
                  card("Features", html.Div(grid("co-t-grid", pd.DataFrame(columns=["feature"]), height=420, fit=True), id="co-t-table"), note="click a row to plot the feature")],
                 className="grid g2"),
        html.Div(style={"height": "16px"}),
        card("Feature by group and block", graph("co-t-box")),
    ])


def tab_models(data):
    ov = data.q("SELECT model, formula, max(note) AS note, max(CASE WHEN layer LIKE 'srna_%' THEN 1 ELSE 0 END) "
                "AS has_srna, max(CASE WHEN layer LIKE 'srna_%' THEN 0 ELSE 1 END) AS has_imp "
                "FROM all_models_overview GROUP BY 1, 2")
    ov["k"] = ov["model"].str.extract(r"(\d+)")[0].astype(float)
    ov = ov.sort_values(["k", "has_srna"])
    src = lambda r: "IMP3 + sRNA" if r.has_srna and r.has_imp else ("sRNA" if r.has_srna else "IMP3")
    opts = [{"label": "%s: %s%s  [%s]" % (r.model, r.formula, " — " + r.note if r.note not in (None, "", "nan")
                                          else "", src(r)),
             "value": "%s|%s" % (r.model, r.formula)} for r in ov.itertuples()]
    return html.Div([
        html.Div([control("model", dd("co-m-model", opts, opts[0]["value"] if opts else None)),
                  control("layer", dd("co-m-layer", [], None)),
                  control("term", dd("co-m-term", [], None))], className="controls"),
        html.Div(id="co-m-info", className="hint"),
        html.Div([card("Volcano", graph("co-m-volcano"), note="OLS on CLR (ratio layers: log2 mt/mg; metrics: own "
                       "scale), HC3 robust SE; numeric terms per SD, categorical vs the reference level"),
                  card("Features", html.Div(grid("co-m-grid", pd.DataFrame(columns=["feature"]), height=420, fit=True), id="co-m-table"), note="BH-FDR per model, layer and term; click a row")],
                 className="grid g2"),
        html.Div(style={"height": "16px"}),
        card("Feature by group and block", graph("co-m-box"), note="raw values of the release (relative "
             "abundance, share of reads on MAGs, mt/mg ratio); the model itself uses CLR and the covariates"),
    ])


def tab_beta(data):
    lays = ordered(data, "all_ordination", "layer")
    return html.Div([
        html.Div([control("layer", dd("co-b-layer", [{"label": layer_label(data, k), "value": k} for k in lays],
                                      lays[0] if lays else None)),
                  control("method", dd("co-b-method", ["PCoA", "NMDS"], "PCoA"))], className="controls"),
        html.Div([card("Ordination (Bray-Curtis)", graph("co-b-ord")),
                  html.Div([card("PERMANOVA", html.Div(id="co-b-perm")),
                            card("Dispersion", html.Div(id="co-b-disp"), note="distance to the own centroid; "
                                 "ratio > 1: first group more variable")])], className="grid g2"),
    ])


def tab_qc(data):
    mets = data.manifest.get("cohort", {}).get("metrics", {})
    return html.Div([
        card("Covariates by group", html.Div(id="co-q-cov"), note="Kruskal-Wallis / χ² overall and per comparison; "
             "imbalanced covariates belong in the models"),
        html.Div(style={"height": "16px"}),
        html.Div([card("Metric tests", html.Div(id="co-q-tests"), note="click a row to plot the metric"),
                  card("Metric by group and block", html.Div([dd("co-q-metric", [{"label": v, "value": k}
                                                                               for k, v in mets.items()],
                                                                 next(iter(mets), None)), graph("co-q-box")]))],
                 className="grid g2"),
        html.Div(style={"height": "16px"}),
        card("QC flags", html.Div(id="co-q-flags"), note="FAIL: left out; WARN: kept (mt layers without samples "
             "lacking mt reads)"),
    ])


def layout(data):
    st, sr = data.manifest.get("cohort", {}), srna(data)
    parts = []
    if st:
        parts.append("IMP3 layers from imp-explorer cohort (%s): %d of %d samples, groups %s (reference %s)%s, "
                     "%d models; full report: cohort/cohort.report.html of the release." % (
                         st.get("created", "?")[:10], st.get("kept", 0), st.get("samples", 0),
                         ", ".join(st.get("groups", [])), st.get("reference", "–"),
                         ", block %s" % st["block"] if st.get("block") else "", len(st.get("models", []))))
    if sr:
        parts.append("sRNA layers (%s) from the sRNA-IMP cohort analysis, as in its report (%d samples, %d in this "
                     "release, %d models)." % (", ".join(sr.get("fractions", {}).values()), sr.get("samples", 0),
                                                 sr.get("samples_in_release", 0), len({m["model"] for m in sr.get("models", [])})))
    lead = " ".join(parts) + (" The global 'group by' and filter controls do not apply here; the analyses are "
                              "fixed. Models with the same number and formula are shown together.")
    return html.Div([
        html.H2("Cohort"), html.P(lead, className="lead"),
        dcc.Tabs(id="co-tabs", value="tests", children=[
            dcc.Tab(label="Group comparisons", value="tests", children=tab_tests(data)),
            dcc.Tab(label="Adjusted models", value="models", children=tab_models(data)),
            dcc.Tab(label="Beta diversity", value="beta", children=tab_beta(data)),
            dcc.Tab(label="QC and covariates", value="qc", children=tab_qc(data)),
        ], className="tabs"),
    ])


def register(app, data):
    q = settings(data).get("sig_q", 0.1)

    # ---------------- group comparisons
    @app.callback(Output("co-t-volcano", "figure"), Output("co-t-table", "children"),
                  Input("co-t-layer", "value"), Input("co-t-comp", "value"))
    def _tests(layer, comp):
        if not layer or not comp:
            return empty_fig(), None
        d = data.q("SELECT * FROM all_tests WHERE layer=? AND comparison=? ORDER BY p", [layer, comp])
        kw = comp == "all groups"
        d["x"] = d["effect"] if kw else d["effect"] - 0.5
        fig = volcano(d, "x", "H" if kw else "P(a > b) − 0.5", q, "%s · %s" % (layer_label(data, layer), comp))
        cols = ["feature", "name", "medians", "means", "prevalence", "p", "q"] if kw else \
            ["feature", "name", "n_a", "n_b", "median_a", "median_b", "mean_a", "mean_b", "effect", "p", "q"]
        return fig, table("co-t-grid", d, cols, layer)

    @app.callback(Output("co-t-box", "figure"), Input("co-t-grid", "selectedRows"), Input("co-t-layer", "value"),
                  Input("co-t-comp", "value"))
    def _tbox(rows, layer, comp):
        """the clicked feature; until a row is clicked (and after a new layer / comparison) the top hit"""
        if rows and ctx.triggered_id == "co-t-grid":
            return feature_box(data, layer, rows[0]["feature"])
        if not layer or not comp:
            return empty_fig("choose a layer and a comparison")
        top = data.q("SELECT feature FROM all_tests WHERE layer=? AND comparison=? ORDER BY p NULLS LAST LIMIT 1",
                     [layer, comp])
        if top.empty:
            return empty_fig("no features for this selection")
        return feature_box(data, layer, top["feature"].iloc[0], subtitle=" · top hit (click a row for another)")

    # ---------------- models
    @app.callback(Output("co-m-layer", "options"), Output("co-m-layer", "value"), Input("co-m-model", "value"))
    def _mlayers(key):
        if not key:
            return [], None
        model, formula = key.split("|", 1)
        ov = data.q("SELECT layer, n, features, features_q_focus FROM all_models_overview WHERE model=? "
                    "AND formula=? AND n > 0", [model, formula])
        opts = [{"label": "%s (n=%d, %d hits)" % (layer_label(data, r.layer), r.n, r.features_q_focus),
                 "value": r.layer} for r in ov.itertuples()]
        best = ov.sort_values("features_q_focus", ascending=False)["layer"].iloc[0] if len(ov) else None
        return opts, best

    @app.callback(Output("co-m-term", "options"), Output("co-m-term", "value"), Output("co-m-info", "children"),
                  Input("co-m-model", "value"), Input("co-m-layer", "value"))
    def _mterms(key, layer):
        if not key or not layer:
            return [], None, None
        model = key.split("|", 1)[0]
        t = data.q("SELECT term, bool_or(focus) AS focus FROM all_models WHERE model=? AND layer=? GROUP BY 1",
                   [model, layer])
        t = t.sort_values("focus", ascending=False)
        opts = [{"label": ("★ " if f else "") + x, "value": x} for x, f in zip(t["term"], t["focus"])]
        info = data.q("SELECT info FROM all_models_overview WHERE model=? AND layer=?", [model, layer])
        txt = info["info"].iloc[0] if len(info) and info["info"].iloc[0] not in (None, "nan", "") else None
        return opts, (t["term"].iloc[0] if len(t) else None), ("design: " + txt) if txt else None

    @app.callback(Output("co-m-volcano", "figure"), Output("co-m-table", "children"),
                  Input("co-m-model", "value"), Input("co-m-layer", "value"), Input("co-m-term", "value"))
    def _models(key, layer, term):
        if not (key and layer and term):
            return empty_fig(), None
        model = key.split("|", 1)[0]
        d = data.q("SELECT * FROM all_models WHERE model=? AND layer=? AND term=? ORDER BY p", [model, layer, term])
        fig = volcano(d, "estimate", "estimate", q, "%s · %s · %s" % (model, layer_label(data, layer), term))
        return fig, table("co-m-grid", d, ["feature", "name", "n", "estimate", "se", "p", "q"], layer)

    @app.callback(Output("co-m-box", "figure"), Input("co-m-grid", "selectedRows"), Input("co-m-model", "value"),
                  Input("co-m-layer", "value"), Input("co-m-term", "value"))
    def _mbox(rows, key, layer, term):
        if rows and ctx.triggered_id == "co-m-grid":
            return feature_box(data, layer, rows[0]["feature"])
        if not (key and layer and term):
            return empty_fig("choose a model, layer and term")
        top = data.q("SELECT feature FROM all_models WHERE model=? AND layer=? AND term=? ORDER BY p NULLS LAST "
                     "LIMIT 1", [key.split("|", 1)[0], layer, term])
        if top.empty:
            return empty_fig("no features for this selection")
        return feature_box(data, layer, top["feature"].iloc[0], subtitle=" · top hit (click a row for another)")

    # ---------------- beta diversity
    @app.callback(Output("co-b-ord", "figure"), Output("co-b-perm", "children"), Output("co-b-disp", "children"),
                  Input("co-b-layer", "value"), Input("co-b-method", "value"))
    def _beta(layer, method):
        if not layer:
            return empty_fig(), None, None
        o = data.q("SELECT * FROM all_ordination WHERE layer=? AND method=?", [layer, method])
        if o.empty:                                    # sRNA-IMP writes NMDS only, imp-explorer cohort PCoA (+ NMDS)
            o = data.q("SELECT * FROM all_ordination WHERE layer=?", [layer])
            method = o["method"].iloc[0] if len(o) else method
        gc = gcolors(data)
        syms = ["circle", "triangle-up", "square", "diamond", "cross", "x"]
        bl = sorted(o["block"].astype(str).unique()) if len(o) else []
        fig = go.Figure()
        for (g, b), x in o.groupby(["group", "block"]):
            fig.add_trace(go.Scatter(x=x["axis1"], y=x["axis2"], mode="markers", name="%s · %s" % (g, b),
                                     text=x["sample_id"], marker=dict(size=8, color=gc.get(str(g), "#555"),
                                                                      symbol=syms[bl.index(str(b)) % len(syms)],
                                                                      line=dict(width=0.5, color="#fff")),
                                     hovertemplate="%{text}<extra>" + "%s · %s" % (g, b) + "</extra>"))
        if len(o) and method == "PCoA":
            fig.update_xaxes(title="PCo1 (%.0f %%)" % (100 * o["explained1"].iloc[0]))
            fig.update_yaxes(title="PCo2 (%.0f %%)" % (100 * o["explained2"].iloc[0]))
            fig.update_layout(title=layer_label(data, layer))
        elif len(o):
            fig.update_layout(title="%s · stress %.3f" % (layer_label(data, layer), o["explained1"].iloc[0]))
        else:
            return empty_fig("no %s for this layer (cohort --nmds)" % method), None, None
        pm = data.q("SELECT comparison, n, F, R2, p FROM all_permanova WHERE layer=? ORDER BY CASE kind WHEN 'within' "
                    "THEN 0 WHEN 'pooled' THEN 0 WHEN 'block' THEN 1 WHEN 'groups' THEN 2 ELSE 3 END", [layer])
        if pm.empty and is_srna(layer):                  # sRNA-IMP: matched PERMANOVA on genus profiles
            pm = data.q("SELECT comparison, n, F, R2, p FROM all_permanova WHERE layer=?",
                        [layer.split(".")[0] + ".genus"])
        ds = data.q("SELECT comparison, n_a, n_b, median_dist_a, median_dist_b, ratio, p FROM all_dispersion "
                    "WHERE layer=?", [layer])
        return style(fig, 460), grid("co-b-pgrid", pm, height=230, fit=True), grid("co-b-dgrid", ds, height=170,
                                                                                    fit=True)

    # ---------------- QC and covariates
    @app.callback(Output("co-q-cov", "children"), Output("co-q-tests", "children"), Output("co-q-flags", "children"),
                  Input("co-tabs", "value"))
    def _qc(tab):
        if tab != "qc":
            raise PreventUpdate
        cov = data.q("SELECT * FROM cohort_covariates") if data.has("cohort_covariates") else (
            data.q("SELECT * FROM srna_covariates") if data.has("srna_covariates") else pd.DataFrame())
        mt = data.q("SELECT feature, name, comparison, n_a, n_b, median_a, median_b, effect, p, q FROM "
                    "cohort_metric_tests ORDER BY p") if data.has("cohort_metric_tests") else pd.DataFrame()
        fl = data.q("SELECT sample_id, \"group\", block, level, excluded, reasons FROM cohort_qc WHERE level <> 'PASS'") \
            if data.has("cohort_qc") else pd.DataFrame()
        if data.has("srna_samples"):
            sf = data.q("SELECT sample_id, \"group\", block, qc AS level, in_statistics, type AS fraction, rules, details "
                        "FROM srna_samples WHERE qc <> 'PASS'")
            fl = pd.concat([fl, sf], ignore_index=True)
        return (grid("co-q-cgrid", cov, height=260, fit=True), grid("co-q-tgrid", mt, height=380, fit=True),
                grid("co-q-fgrid", fl, height=260, fit=True))

    @app.callback(Output("co-q-metric", "value"), Input("co-q-tgrid", "selectedRows"), prevent_initial_call=True)
    def _qpick(rows):
        if not rows:
            raise PreventUpdate
        return rows[0]["feature"]

    @app.callback(Output("co-q-box", "figure"), Input("co-q-metric", "value"), Input("co-tabs", "value"))
    def _qbox(metric, tab):
        if tab != "qc" or not metric:
            raise PreventUpdate
        return feature_box(data, "metrics", metric)


PAGE = Page(path="/cohort", name="Cohort", order=70, layout=layout, register=register, requires=("all_tests",))
