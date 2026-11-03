"""Overview: what is in the release, samples with their metadata and key numbers."""
import pandas as pd
import plotly.graph_objects as go
from dash import Input, Output, html

from ..ui import GLOBAL, Page, card, colors, fmt_int, graph, stat_tile, style
from ._common import grid, selection

KEY_QC = [("mg", "trimmed:pairs", "mg read pairs (trimmed)"), ("mt", "trimmed:pairs", "mt read pairs (trimmed)"),
          ("mg", "mapped_fraction", "mg mapped to assembly"), ("mt", "mapped_fraction", "mt mapped to assembly"),
          ("mgmt", "n50", "assembly N50"), ("mg", "nonpareil_diversity", "Nonpareil diversity (mg)")]


def srna_tiles(data):
    """second row of key numbers when sRNA-IMP results are imported (imp-explorer srna)"""
    if not data.has("srna_samples"):
        return None
    man = data.manifest.get("srna", {})
    fr = man.get("fractions", {})
    lib = data.q("SELECT omic, count(DISTINCT sample_id) AS n FROM srna_samples GROUP BY 1 ORDER BY 1")
    ids = set(data.q("SELECT DISTINCT sample_id FROM srna_samples")["sample_id"].astype(str))
    inrel = len(ids & set(data.samples["sample_id"].astype(str)))
    tiles = [stat_tile(fmt_int(len(ids)), "sRNA samples (%s in this release)" % fmt_int(inrel)),
             stat_tile(" / ".join(fmt_int(n) for n in lib["n"]),
                       "libraries: " + " / ".join(fr.get(o, o) for o in lib["omic"]))]
    if data.has("srna_features"):
        nf = data.q('SELECT "set", count(DISTINCT feature) AS n FROM srna_features GROUP BY 1')
        nf = dict(zip(nf["set"], nf["n"]))
        for st, lab in (("mirna", "human miRNAs"), ("genus", "genera (sRNA rRNA reads)")):
            if st in nf:
                tiles.append(stat_tile(fmt_int(nf[st]), lab))
    if data.has("srna_tests"):
        q = man.get("sig_q", 0.1)
        n = data.q("SELECT count(DISTINCT layer || '|' || feature) AS n FROM srna_tests WHERE q < ? "
                   "AND layer NOT LIKE '%.metrics'", [q])["n"].iloc[0]
        tiles.append(stat_tile(fmt_int(n), "sRNA features with q < %g in a group comparison" % q))
    return html.Div(tiles, className="tiles")


def metab_tiles(data):
    """key numbers of the metabolomics assays (imp-explorer metabolomics)"""
    if not data.has("metab_assays"):
        return None
    a = data.q("SELECT * FROM metab_assays")
    ids = set(data.q("SELECT DISTINCT sample_id FROM metab_values")["sample_id"].astype(str))
    inrel = len(ids & set(data.samples["sample_id"].astype(str)))
    nf = int(pd.to_numeric(a.get("n_features"), errors="coerce").fillna(0).sum()) if "n_features" in a else 0
    tiles = [stat_tile(fmt_int(len(ids)), "metabolomics samples (%s in this release)" % fmt_int(inrel)),
             stat_tile(fmt_int(len(a)), "assays: " + ", ".join(a["assay"].astype(str))),
             stat_tile(fmt_int(nf), "features measured")]
    if data.has("metab_features") and "annotation" in data.q("SELECT * FROM metab_features LIMIT 1").columns:
        n = data.q("SELECT count(*) AS n FROM metab_features WHERE NOT (lower(coalesce(annotation, '')) LIKE "
                   "'unknown%' OR lower(coalesce(annotation, '')) LIKE 'name from mass%')")["n"].iloc[0]
        tiles.append(stat_tile(fmt_int(n), "identified features"))
    return html.Div(tiles, className="tiles")


def sample_table(data):
    q = data.q("SELECT sample_id, omic, metric, value FROM qc") if data.has("qc") else None
    s = data.samples.drop(columns=[c for c in ("imp_dir",) if c in data.samples.columns]).copy()
    if q is not None and len(q):
        for omic, metric, label in KEY_QC:
            v = q[(q["omic"] == omic) & (q["metric"] == metric)].set_index("sample_id")["value"]
            s[label] = s["sample_id"].map(v)
    if data.has("bins"):
        b = data.q("SELECT sample_id, quality, count(*) AS n FROM bins GROUP BY 1, 2")
        for tier in ("high", "medium"):
            s["MAGs %s" % tier] = s["sample_id"].map(b[b["quality"] == tier].set_index("sample_id")["n"]).fillna(0)
    if data.has("srna_samples"):
        fr = data.manifest.get("srna", {}).get("fractions", {})
        lib = data.q("SELECT DISTINCT sample_id, omic FROM srna_samples")
        for o in sorted(lib["omic"].unique()):
            have = set(lib.loc[lib["omic"] == o, "sample_id"].astype(str))
            s[fr.get(o, o)] = s["sample_id"].astype(str).map(lambda x: "yes" if x in have else "")
    if data.has("metab_values"):
        m = data.q("SELECT DISTINCT sample_id, assay FROM metab_values")
        for assay in sorted(m["assay"].unique()):
            have = set(m.loc[m["assay"] == assay, "sample_id"].astype(str))
            s["metab " + assay] = s["sample_id"].astype(str).map(lambda x: "yes" if x in have else "")
    return s


def layout(data):
    n = len(data.samples)
    nb = data.q("SELECT quality, count(*) AS n FROM bins GROUP BY 1") if data.has("bins") else None
    tiers = dict(zip(nb["quality"], nb["n"])) if nb is not None else {}
    ng = data.q("SELECT count(DISTINCT taxon) AS n FROM taxonomy WHERE rank='genus'")["n"].iloc[0] \
        if data.has("taxonomy") else 0
    nk = data.q("SELECT count(DISTINCT feature) AS n FROM functions WHERE db='kegg_ko'")["n"].iloc[0] \
        if data.has("functions") else 0
    man = data.manifest
    tiles = html.Div([
        stat_tile(fmt_int(n), "samples"),
        stat_tile("%s / %s" % (fmt_int(data.samples.get("has_mg", []).sum()), fmt_int(data.samples.get("has_mt", []).sum())),
                  "with metagenome / metatranscriptome"),
        stat_tile(fmt_int(sum(tiers.values())), "MAGs (bins)"),
        stat_tile("%s / %s" % (fmt_int(tiers.get("high", 0)), fmt_int(tiers.get("medium", 0))), "high / medium quality"),
        stat_tile(fmt_int(ng), "genera (read-based)"),
        stat_tile(fmt_int(nk), "KEGG orthologs"),
    ], className="tiles")
    srna = srna_tiles(data)
    metab = metab_tiles(data)
    return html.Div([
        html.H2("Overview"),
        html.P("Release %s of %s, exported %s from %s. Use 'group by' and 'filter' in the top bar; every page "
               "follows them." % (man.get("release", "?"), data.study, (man.get("created") or "?")[:10],
                                   man.get("sources", {}).get("pipeline", "IMP3")), className="lead"),
        tiles,
        *([html.Div("Small RNA (sRNA-IMP)", className="tile-label", style={"margin": "4px 0 6px"}), srna]
          if srna is not None else []),
        *([html.Div("Metabolomics", className="tile-label", style={"margin": "4px 0 6px"}), metab]
          if metab is not None else []),
        html.Div([card("Samples per group", graph("ov-groups")),
                  card("MAGs per sample", graph("ov-mags"), note="high (≥ 90 % complete, ≤ 5 % contamination) and "
                       "medium quality (≥ 50 %, ≤ 10 %) bins of each sample, MIMAG tiers")],
                 className="grid g2"),
        html.Div(style={"height": "16px"}),
        card("Samples", html.Div(id="ov-table"), note="metadata and key numbers; sort and filter in the column headers"),
    ])


def register(app, data):
    @app.callback(Output("ov-groups", "figure"), Output("ov-mags", "figure"), Output("ov-table", "children"),
                  Input(GLOBAL[0], "value"), Input(GLOBAL[1], "value"), Input(GLOBAL[2], "value"))
    def _update(group, fcol, fval):
        s, g, cmap = selection(data, group, fcol, fval)
        counts = g.value_counts().sort_index()
        f1 = go.Figure(go.Bar(x=[str(i) for i in counts.index], y=counts.values,
                              marker_color=[cmap.get(str(i)) for i in counts.index], text=counts.values,
                              textposition="outside"))
        f1.update_yaxes(title="samples", range=[0, max(counts.values) * 1.15 if len(counts) else 1])
        style(f1, 300, legend=False)
        f2 = go.Figure()
        if data.has("bins"):
            b = data.q("SELECT sample_id, quality, count(*) AS n FROM bins GROUP BY 1, 2")
            b = b[b["sample_id"].isin(s["sample_id"])]
            piv = b.pivot_table(index="sample_id", columns="quality", values="n", fill_value=0)
            piv["g"] = g.reindex(piv.index).fillna("NA")
            piv = piv.sort_values(["g"] + [c for c in ("high",) if c in piv.columns], ascending=[True, False])
            for tier, col in (("high", "#2a6f97"), ("medium", "#8fb8d0"), ("low", "#d9dde2")):
                if tier in piv.columns:
                    f2.add_trace(go.Bar(x=list(piv.index), y=piv[tier], name=tier, marker_color=col))
            f2.update_layout(barmode="stack")
            f2.update_xaxes(showticklabels=False, title="samples (sorted by group, then high-quality MAGs)")
            f2.update_yaxes(title="MAGs")
        style(f2, 300)
        t = sample_table(data)
        t = t[t["sample_id"].isin(s["sample_id"])]
        return f1, f2, grid("ov-grid", t, height=420)


PAGE = Page(path="/", name="Overview", order=0, layout=layout, register=register, requires=("samples",))
