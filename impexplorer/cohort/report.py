"""Cohort report: PDF (matplotlib) and HTML with the same sections, figures and tables."""
import base64
import html
import io
import math
import os
import textwrap

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                       # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

from .data import METRIC_LABEL                        # noqa: E402

PALETTE = ["#2a6f97", "#d1495b", "#edae49", "#00798c", "#8f5aa8", "#66a182", "#e07a5f", "#3d405b"]
MARKERS = ["o", "^", "s", "D", "v", "P"]
TAXA_PAL = ["#2a6f97", "#d1495b", "#edae49", "#00798c", "#8f5aa8", "#66a182", "#e07a5f", "#3d405b", "#81b29a",
            "#f2cc8f", "#6d597a", "#b56576", "#457b9d", "#9c6644"]
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 9, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.titleweight": "bold"})
A4 = (8.27, 11.69)
A4L = (11.69, 8.27)

EXPLAIN = {
    "design":
        "Samples with IMP3 output by group and block (e.g. site), and how many were used after QC. In a nested design "
        "(a group sampled at one block only) group and block effects cannot be separated: every group is compared with "
        "the reference group within the blocks they share, and the reference group split by block shows how large "
        "block (site) differences are without disease.",
    "findings":
        "Automatic list of results that pass the thresholds, by area. It is a guide to the report pages, not a result "
        "list: many tests are run; a group difference that is as large as the block contrast within the reference "
        "group, or that goes with a covariate that differs between the groups, is likely confounded.",
    "covariates":
        "Balance of the covariates between the groups: overall (Kruskal-Wallis / χ²) and for each comparison. A "
        "covariate that differs between the groups compared (e.g. sex, age, stool consistency, medication) can explain "
        "apparent group differences and belongs in the adjusted models.",
    "qc":
        "QC rules per sample: FAIL without metagenome reads or below the minimum metagenome depth (left out of all "
        "statistics); WARN without (enough) metatranscriptome reads (left out of the mt layers) or as robust outlier "
        "(|robust z| above the threshold) in a key metric. Boxes by group and block: differences in depth, host "
        "share, rRNA share or mapping can reflect sampling and extraction rather than biology.",
    "alpha":
        "Within-sample diversity: Shannon index and richness (taxa with reads) of the read-based profiles. Richness "
        "depends on sequencing depth; compare it with the depth panels of the QC page.",
    "metric_tests":
        "Per-sample metrics compared between groups: each group vs the reference within shared blocks (van Elteren, "
        "stratified Wilcoxon; effect = P(group > reference)), the block contrast within the reference group, and "
        "Kruskal-Wallis across all groups. BH-FDR across metrics within each comparison.",
    "composition":
        "Mean relative abundance of the most abundant taxa per group and block (read-based Kraken2/Bracken, and reads "
        "on MAGs summed per GTDB genus), metagenome (present) and metatranscriptome (active).",
    "beta":
        "Bray-Curtis PCoA of the relative abundances (colour: group, marker: block). PERMANOVA: all groups with "
        "permutations within blocks; group after block and block after group (partial R²); every comparison within its "
        "blocks; the block contrast within the reference group. Metabolomics layers: Euclidean distances of the "
        "autoscaled log2 values instead of Bray-Curtis.",
    "nmds":
        "Non-metric MDS (Kruskal) of the Bray-Curtis distances, best of a PCoA and five random starts; stress below "
        "0.2 is a usable 2-D picture. Samples far from all others are left out of the NMDS (they collapse the rest).",
    "features":
        "Per feature with at least 10 reads in at least the prevalence share of the samples (CLR-free, on relative "
        "abundances): each group vs the reference within shared blocks (van Elteren), the block contrast within the "
        "reference group, Kruskal-Wallis across groups. BH-FDR per layer and comparison. Table: features with q below "
        "the threshold (top 25 by p per comparison).",
    "models":
        "Covariate-adjusted linear models per feature: OLS of the CLR-transformed counts (pseudocount 0.5) on the "
        "model terms, HC3 robust standard errors; numeric terms per SD, categorical terms against the reference level "
        "(group: the reference group). The first term of a formula is its term of interest. 'terms | column=value' "
        "fits the model on these samples only. Ratio layers are modelled directly (log2 mt/mg), metrics on their "
        "scale. BH-FDR per model, layer and term. Fits with fewer than 20 samples, < 3 samples in a level or a sample "
        "with leverage > 0.95 are skipped.",
    "mags":
        "MAGs (medium and high quality, MetaWrap-refined bins with GTDB-Tk taxonomy) per sample and the share of "
        "samples in which a MAG of a genus was recovered, by group. Recovery depends on depth and on assembly; it is "
        "not the same as presence. P values: van Elteren on recovery (yes / no) within blocks.",
    "ratio":
        "Expression relative to gene copies: CLR(mt) - CLR(mg) per feature in log2 units, for features tested in both "
        "omics; > 0 = more transcribed than the average feature relative to its abundance in the metagenome. Left: "
        "median ratio of the most variable features per group. Table: group comparisons (q below the threshold).",
    "dispersion":
        "Heterogeneity: distance of every sample to the centroid of its own side of a comparison (PCoA of the "
        "Bray-Curtis distances, all axes); ratio > 1 = the first group is more variable. Mann-Whitney test, BH-FDR. "
        "Higher dispersion with a similar centroid means patients diverge from controls in different directions.",
    "methods":
        "Settings of this run and the notes of the analysis.",
}


def esc(x):
    return html.escape(str(x))


def fmtv(v):
    if isinstance(v, (float, np.floating)):
        if math.isnan(v):
            return "–"
        a = abs(v)
        return ("%.3g" % v) if (a >= 1e-3 or a == 0) and a < 1e5 else ("%.2e" % v)
    if isinstance(v, (bool, np.bool_)):
        return "yes" if v else ""
    return str(v)


class Report:
    def __init__(self, C, A, res, st, out):
        self.C, self.A, self.res, self.st, self.out = C, A, res, st, out
        self.sections = []           # (key, title, [items]) items: ("fig", fig) | ("table", caption, df) | ("text", s)
        self.q = st.get("sig_q", 0.1)
        self.gcol = {g: PALETTE[i % len(PALETTE)] for i, g in enumerate(C.groups or ["all"])}
        blocks = sorted({C.block_of(s) for s in C.samples.index})
        self.bmark = {b: MARKERS[i % len(MARKERS)] for i, b in enumerate(blocks)}

    def add(self, key, title, *items):
        self.sections.append((key, title, list(items)))

    def label(self, layer):
        return "library metrics" if layer == "metrics" else self.C.layers[layer].label if layer in self.C.layers \
            else layer

    # ------------------------------------------------------------------ plot helpers
    def cells(self, samples):
        C = self.C
        out = []
        for g in (C.groups or ["all"]):
            for b in sorted({C.block_of(s) for s in samples if C.group(s) == g}):
                ids = [s for s in samples if C.group(s) == g and C.block_of(s) == b]
                if ids:
                    out.append((g, b, ids))
        return out

    def box(self, ax, values, title, ylab=None, log=False):
        v = values.dropna()
        if log:
            v = v[v > 0]
        cells = self.cells(list(v.index))
        data = [v.loc[ids].values for g, b, ids in cells]
        if not data:
            ax.set_visible(False)
            return
        bp = ax.boxplot(data, widths=0.6, showfliers=False, patch_artist=True)
        for i, (g, b, ids) in enumerate(cells):
            bp["boxes"][i].set(facecolor="none", edgecolor=self.gcol.get(g, "#555"))
            bp["medians"][i].set(color=self.gcol.get(g, "#555"))
            x = np.random.default_rng(i).normal(i + 1, 0.07, len(ids))
            ax.scatter(x, v.loc[ids].values, s=6, color=self.gcol.get(g, "#555"), alpha=0.6,
                       marker=self.bmark.get(b, "o"), linewidths=0)
        ax.set_xticks(range(1, len(cells) + 1))
        ax.set_xticklabels(["%s\n%s" % (g, b) if self.C.block else g for g, b, _ in cells], fontsize=6.5)
        if log:
            ax.set_yscale("log")
        ax.set_title(textwrap.fill(title, 34), fontsize=7.5)
        if ylab:
            ax.set_ylabel(ylab)

    def box_grid(self, metrics, ncol=4):
        M = self.res["sample_metrics"].set_index("sample_id")
        metrics = [m for m in metrics if m in M.columns]
        if not metrics:
            return None
        nrow = int(math.ceil(len(metrics) / ncol))
        fig, axs = plt.subplots(nrow, ncol, figsize=(A4[0], 2.3 * nrow + 0.4))
        axs = np.atleast_1d(axs).ravel()
        for ax, m in zip(axs, metrics):
            self.box(ax, M[m].astype(float), METRIC_LABEL.get(m, m))
        for ax in axs[len(metrics):]:
            ax.set_visible(False)
        fig.tight_layout()
        return fig

    def legend_fig(self, fig):
        hs = [plt.Line2D([], [], color=c, marker="o", ls="", label=g) for g, c in self.gcol.items()]
        if self.C.block:
            hs += [plt.Line2D([], [], color="#555", marker=m, ls="", label=b) for b, m in self.bmark.items()]
        fig.legend(handles=hs, loc="lower center", ncol=min(len(hs), 8), frameon=False, fontsize=7)

    # ------------------------------------------------------------------ sections
    def s_design(self):
        d = self.res["design"].rename(columns={"block": self.C.block or "block"})
        qs = self.res["qc_summary"]
        txt = "%s · %d samples with IMP3 output, %d used. Groups: %s (reference %s). Block: %s." % (
            self.st["name"], len(self.C.samples), len(self.C.kept), ", ".join(self.C.groups) or "–",
            self.C.ref or "–", self.C.block or "–")
        self.add("design", "Samples and design", ("text", txt), ("table", "Samples per group and block", d),
                 ("table", "QC by group (exclusion vs group: χ² p = %s)" % fmtv(qs.attrs.get("p_excluded", np.nan)),
                  qs))

    def s_findings(self):
        f = self.res["findings"]
        df = pd.DataFrame(f, columns=["area", "finding"]) if f else pd.DataFrame({"finding": ["no result passes"]})
        self.add("findings", "Key findings", ("table", "", df))

    def s_covariates(self):
        cv = self.res["covariates"]
        if len(cv):
            self.add("covariates", "Covariates by group", ("table", "", cv))

    def s_qc(self):
        items = []
        fig = self.box_grid(["mg_pairs_M", "mt_pairs_M", "mg_removed_fraction", "mt_rrna_fraction",
                             "mg_mapped_fraction", "mt_mapped_fraction", "mg_nonpareil", "assembly_Mb",
                             "assembly_n50_kb", "mags_hq", "mg_binned_fraction", "mt_binned_fraction"])
        if fig:
            self.legend_fig(fig)
            fig.subplots_adjust(bottom=0.08)
            items.append(("fig", fig))
        fl = self.res["qc_flags"]
        fl = fl[fl["level"] != "PASS"][["sample_id", "group", "block", "level", "excluded", "reasons"]]
        items.append(("table", "Samples with WARN or FAIL (%d)" % len(fl), fl))
        self.add("qc", "QC", *items)

    def s_alpha(self):
        ms = [m for m in self.res["sample_metrics"].columns if m.startswith(("shannon_", "richness_"))]
        fig = self.box_grid(ms, ncol=4)
        if fig:
            self.legend_fig(fig)
            fig.subplots_adjust(bottom=0.1)
            self.add("alpha", "Alpha diversity", ("fig", fig))

    def s_metric_tests(self):
        mt = self.res["metrics_tests"]
        if not len(mt):
            return
        t = mt[mt["kind"] != "kruskal"].copy()
        piv = t.pivot_table(index="feature", columns="comparison", values="p", aggfunc="first")
        eff = t.pivot_table(index="feature", columns="comparison", values="effect", aggfunc="first")
        qq = t.pivot_table(index="feature", columns="comparison", values="q", aggfunc="first")
        kw = mt[mt["kind"] == "kruskal"].set_index("feature")
        rows = []
        for f in piv.index:
            r = {"metric": METRIC_LABEL.get(f, f)}
            if f in kw.index:
                r["all groups p"] = kw.at[f, "p"]
            for c in piv.columns:
                r[c] = "%s (%s)%s" % (fmtv(eff.at[f, c]), fmtv(piv.at[f, c]), " *" if qq.at[f, c] < self.q else "")
            rows.append(r)
        self.add("metric_tests", "Metrics by group",
                 ("table", "cells: effect P(group > reference) (p); * q < %g" % self.q, pd.DataFrame(rows)))

    def s_composition(self):
        lays = [k for k in ("kraken.genus.mg", "kraken.genus.mt", "mag.genus.mg", "mag.genus.mt") if k in self.C.layers]
        if not lays:
            return
        fig, axs = plt.subplots(len(lays), 1, figsize=(A4[0], 2.5 * len(lays) + 0.5))
        axs = np.atleast_1d(axs)
        for ax, k in zip(axs, lays):
            R = self.C.layers[k].rel()
            cells = self.cells(list(R.index))
            top = list(R.mean().sort_values(ascending=False).index[:12])
            mean = pd.DataFrame([R.loc[ids, top].mean() for g, b, ids in cells])
            bottom = np.zeros(len(cells))
            for i, t in enumerate(top):
                ax.bar(range(len(cells)), mean[t].values, bottom=bottom, color=TAXA_PAL[i % len(TAXA_PAL)], label=t,
                       width=0.75)
                bottom += mean[t].values
            ax.bar(range(len(cells)), 1 - bottom, bottom=bottom, color="#d9dde2", label="other", width=0.75)
            ax.set_xticks(range(len(cells)))
            ax.set_xticklabels(["%s · %s (n=%d)" % (g, b, len(ids)) if self.C.block else "%s (n=%d)" % (g, len(ids))
                                for g, b, ids in cells], fontsize=7)
            ax.set_ylim(0, 1)
            ax.set_ylabel("mean relative abundance")
            ax.set_title(self.C.layers[k].label, loc="left")
            ax.legend(fontsize=6, bbox_to_anchor=(1.01, 1), loc="upper left", frameon=False)
        fig.tight_layout()
        self.add("composition", "Taxonomic composition", ("fig", fig))

    def scatter_grid(self, method, layers, title_fn):
        o = self.res["ordination"]
        o = o[(o["method"] == method) & o["layer"].isin(layers)]
        lays = [k for k in layers if k in set(o["layer"])]
        if not lays:
            return None
        ncol = 2
        nrow = int(math.ceil(len(lays) / ncol))
        fig, axs = plt.subplots(nrow, ncol, figsize=(A4[0], 3.3 * nrow + 0.5))
        axs = np.atleast_1d(axs).ravel()
        for ax, k in zip(axs, lays):
            d = o[o["layer"] == k]
            for (g, b), x in d.groupby(["group", "block"]):
                ax.scatter(x["axis1"], x["axis2"], s=12, color=self.gcol.get(g, "#555"), marker=self.bmark.get(b, "o"),
                           alpha=0.75, linewidths=0)
            e1, e2 = d["explained1"].iloc[0], d["explained2"].iloc[0]
            if method == "PCoA":
                ax.set_xlabel("PCo1 (%.0f %%)" % (100 * e1))
                ax.set_ylabel("PCo2 (%.0f %%)" % (100 * e2))
            else:
                ax.set_xlabel("NMDS1")
                ax.set_ylabel("NMDS2")
            ax.set_title(textwrap.fill(title_fn(k, d), 60), fontsize=7.5)
        for ax in axs[len(lays):]:
            ax.set_visible(False)
        self.legend_fig(fig)
        fig.tight_layout(rect=(0, 0.04, 1, 1))
        return fig

    def s_beta(self):
        pm = self.res["permanova"]
        main = [k for k in ("kraken.genus.mg", "kraken.genus.mt", "motus.species.mg", "motus.species.mt",
                            "mag.genus.mg", "mag.genus.mt", "kegg_pathway.mg", "kegg_pathway.mt")
                if k in self.C.layers]

        def title(k, d):
            x = pm[(pm["layer"] == k) & (pm["kind"].isin(["within", "block", "pooled"]))]
            sub = "; ".join("%s R² %.3f p %s" % (r.comparison, r.R2, fmtv(r.p)) for r in x.itertuples())
            return "%s\n%s" % (self.C.layers[k].label, sub)
        items = []
        fig = self.scatter_grid("PCoA", main, title)
        if fig:
            items.append(("fig", fig))
        t = pm.copy()
        t["layer"] = t["layer"].map(self.label)
        items.append(("table", "PERMANOVA (Bray-Curtis; metabolomics: Euclidean of autoscaled log2 values; %d permutations)" % self.st["permutations"], t))
        self.add("beta", "Beta diversity", *items)
        o = self.res["ordination"]
        nm = [k for k in self.C.layers if ((o["method"] == "NMDS") & (o["layer"] == k)).any()] if len(o) else []
        if nm:
            fig = self.scatter_grid("NMDS", nm, lambda k, d: "%s · stress %.3f" % (self.C.layers[k].label,
                                                                                    d["explained1"].iloc[0]))
            self.add("nmds", "NMDS", ("fig", fig))

    def s_features(self):
        ft = self.res["feature_tests"]
        if not len(ft):
            return
        t = ft[ft["kind"] != "kruskal"]
        summ = t.assign(sig=t["q"] < self.q).groupby(["layer", "comparison"]).agg(
            tested=("feature", "size"), sig=("sig", "sum")).reset_index()
        piv = summ.pivot_table(index="layer", columns="comparison", values="sig", aggfunc="first").fillna(0).astype(int)
        nt = summ.groupby("layer")["tested"].max()
        piv.insert(0, "tested", nt)
        piv.index = [self.label(k) for k in piv.index]
        items = [("table", "features with q < %g per layer and comparison" % self.q, piv.reset_index().rename(
            columns={"index": "layer"}))]
        for comp, d in t.groupby("comparison", sort=False):
            s = d[d["q"] < self.q].sort_values("p").groupby("layer").head(25)
            if not len(s):
                continue
            s = s.assign(layer=s["layer"].map(self.label),
                         direction=np.where(s["effect"] > 0.5, "higher", "lower"))
            items.append(("table", comp, s[["layer", "feature", "name", "n_a", "n_b", "median_a", "median_b",
                                            "effect", "direction", "p", "q"]]))
        self.add("features", "Group comparisons per feature", *items)

    def s_models(self):
        mo, ov = self.res["models"], self.res["models_overview"]
        specs = self.res["model_specs"]
        if not specs:
            return
        piv = ov.pivot_table(index="layer", columns="model", values="features_q_focus", aggfunc="first")
        piv = piv.reindex([c for c in ["metrics"] + list(self.C.layers) if c in piv.index])
        piv = piv[[s["model"] for s in specs if s["model"] in piv.columns]].fillna(0).astype(int)
        piv.index = [self.label(c) for c in piv.index]
        piv.index.name = "layer"
        head = pd.DataFrame([{"model": s["model"], "formula": s["formula"], "question": s["note"],
                              "n (metrics)": int(ov[(ov["model"] == s["model"]) & (ov["layer"] == "metrics")]["n"].max()
                                                 if len(ov[(ov["model"] == s["model"]) & (ov["layer"] == "metrics")])
                                                 else 0)} for s in specs])
        items = [("table", "models", head),
                 ("table", "features with the term of interest at q < %g (rows: layers, columns: models)" % self.q,
                  piv.reset_index())]
        info = ov[ov["info"].fillna("").astype(str).str.len() > 0][["model", "layer", "info"]].drop_duplicates(
            ["model", "info"])
        if len(info):
            items.append(("table", "design notes (dropped collinear columns, missing terms)", info))
        # volcano per model (focus terms, all layers)
        foc = mo[mo["focus"]] if len(mo) else mo
        vfoc = foc[foc["layer"] != "metrics"] if len(foc) else foc        # metrics are on their own scales
        ms = [s["model"] for s in specs if len(vfoc) and (vfoc["model"] == s["model"]).any()]
        if ms:
            ncol = 2
            nrow = int(math.ceil(len(ms) / ncol))
            fig, axs = plt.subplots(nrow, ncol, figsize=(A4[0], 3.0 * nrow + 0.4))
            axs = np.atleast_1d(axs).ravel()
            kinds = {"tax": "#2a6f97", "mag": "#00798c", "fn": "#edae49", "ratio": "#d1495b"}

            def kind(l):
                return "metrics" if l == "metrics" else "ratio" if l.startswith("ratio") else "mag" if l.startswith(
                    "mag") else "tax" if l.startswith(("kraken", "motus")) else "fn"
            for ax, m in zip(axs, ms):
                d = vfoc[(vfoc["model"] == m) & vfoc["p"].notna()]
                for kd, x in d.groupby(d["layer"].map(kind)):
                    sig = x["q"] < self.q
                    ax.scatter(x["estimate"], -np.log10(x["p"].clip(lower=1e-300)), s=np.where(sig, 10, 3),
                               color=kinds[kd], alpha=np.where(sig, 0.9, 0.25).tolist(), linewidths=0, label=kd)
                spec = next(s for s in specs if s["model"] == m)
                ax.set_title(textwrap.fill("%s: %s%s" % (m, spec["formula"], " — " + spec["note"] if spec["note"]
                                                        else ""), 55), fontsize=7.5)
                ax.set_xlabel("estimate (CLR or log2 ratio; per SD / vs reference)")
                ax.set_ylabel("-log10 p")
                ax.axvline(0, color="#9aa5b1", lw=0.6)
            for ax in axs[len(ms):]:
                ax.set_visible(False)
            hs = [plt.Line2D([], [], color=c, marker="o", ls="", label=k) for k, c in kinds.items()]
            fig.legend(handles=hs, loc="lower center", ncol=4, frameon=False, fontsize=7)
            fig.tight_layout(rect=(0, 0.04, 1, 1))
            items.append(("fig", fig))
        for s in specs:
            d = foc[(foc["model"] == s["model"]) & (foc["q"] < self.q)].sort_values("p") if len(foc) else foc
            if len(d):
                d = d.groupby("layer").head(15)
                d = d.assign(layer=d["layer"].map(self.label))
                items.append(("table", "%s %s — term of interest, q < %g" % (s["model"], s["formula"], self.q),
                              d[["layer", "feature", "name", "n", "term", "estimate", "se", "p", "q"]]))
        self.add("models", "Adjusted models", *items)

    def s_mags(self):
        mp = self.res.get("mag_prevalence")
        items = []
        fig = self.box_grid(["mags_hq", "mags_mq", "mg_binned_fraction", "mt_binned_fraction"], ncol=4)
        if fig:
            self.legend_fig(fig)
            fig.subplots_adjust(bottom=0.2)
            items.append(("fig", fig))
        if mp is not None and len(mp):
            top = mp.head(30)
            gcols = [c for c in mp.columns if c.startswith("prev ")]
            fig, ax = plt.subplots(figsize=(A4[0], 6.5))
            y = np.arange(len(top))
            h = 0.8 / max(len(gcols), 1)
            for i, c in enumerate(gcols):
                ax.barh(y + i * h, top[c].values * 100, height=h, color=self.gcol.get(c[5:], "#555"), label=c[5:])
            ax.set_yticks(y + 0.4 - h / 2)
            ax.set_yticklabels(["%s (%s)" % (g, p) for g, p in zip(top["genus"], top["phylum"])], fontsize=6.5)
            ax.invert_yaxis()
            ax.set_xlabel("% of samples with a MAG of the genus")
            ax.legend(frameon=False, fontsize=7)
            ax.set_title("MAG recovery of the 30 most often recovered genera", loc="left")
            fig.tight_layout()
            items.append(("fig", fig))
            qc = [c for c in mp.columns if c.startswith("q ")]
            sig = mp[(mp[qc] < self.q).any(axis=1)] if qc else mp.iloc[:0]
            items.append(("table", "genera whose MAG recovery differs (q < %g in a comparison)" % self.q,
                          sig.head(40)))
        if items:
            self.add("mags", "MAGs", *items)

    def s_ratio(self):
        lays = [k for k, L in self.C.layers.items() if L.kind == "ratio"]
        if not lays:
            return
        items = []
        fig, axs = plt.subplots(1, len(lays), figsize=(A4[0], 5.2), gridspec_kw={"wspace": 1.6})
        axs = np.atleast_1d(axs)
        for ax, k in zip(axs, lays):
            Y = self.C.layers[k].M
            top = Y.var().sort_values(ascending=False).index[:25]
            groups = self.C.groups or ["all"]
            med = pd.DataFrame({g: Y.loc[[s for s in Y.index if self.C.group(s) == g or g == "all"], top].median()
                                for g in groups})
            im = ax.imshow(med.values, aspect="auto", cmap="RdBu_r", vmin=-3, vmax=3)
            ax.set_yticks(range(len(top)))
            ax.set_yticklabels([str(t)[:24] for t in top], fontsize=5.5)
            ax.set_xticks(range(len(groups)))
            ax.set_xticklabels(groups, fontsize=7)
            ax.set_title(textwrap.fill(self.C.layers[k].label, 28), fontsize=7.5)
        fig.subplots_adjust(left=0.17, right=0.9, top=0.92, bottom=0.08)
        fig.colorbar(im, ax=list(axs), shrink=0.5, label="median log2 mt/mg (relative)")
        items.append(("fig", fig))
        ft = self.res["feature_tests"]
        if len(ft):
            t = ft[ft["layer"].isin(lays) & (ft["kind"] != "kruskal") & (ft["q"] < self.q)].sort_values("p")
            if len(t):
                t = t.assign(layer=t["layer"].map(self.label))
                items.append(("table", "group comparisons of the ratios (q < %g)" % self.q,
                              t[["layer", "feature", "name", "comparison", "median_a", "median_b", "effect", "p",
                                 "q"]].head(60)))
        self.add("ratio", "Functions and MAGs: metatranscriptome / metagenome", *items)

    def s_dispersion(self):
        ds = self.res["dispersion"]
        if len(ds):
            self.add("dispersion", "Heterogeneity (dispersion)", ("table", "", ds.assign(layer=ds["layer"].map(
                self.label))))

    def s_methods(self):
        st = self.st
        rows = [("release", os.path.abspath(self.C.release)), ("group column", st["group_col"]),
                ("groups", ", ".join(st["groups"])), ("reference", st["reference"]), ("block", st["block"] or "–"),
                ("covariates", ", ".join(st["covariates"])), ("QC drop", st["qc_drop"]),
                ("PERMANOVA permutations", st["permutations"]), ("prevalence filter", st["prevalence"]),
                ("FDR threshold", st["sig_q"]), ("settings file", st.get("config") or "–"),
                ("from the settings file", ", ".join(st.get("from_config", [])) or "–")]
        for s in st["models"]:
            rows.append((s["model"], s["formula"] + (" — " + s["note"] if s.get("note") else "")))
        for k, L in st["layers"].items():
            rows.append(("layer " + k, "%s: %d samples x %d features" % (L["label"], L["samples"], L["features"])))
        for n in st["notes"]:
            rows.append(("note", n))
        self.add("methods", "Settings and notes", ("table", "", pd.DataFrame(rows, columns=["setting", "value"])))

    def build(self):
        for fn in (self.s_design, self.s_findings, self.s_covariates, self.s_qc, self.s_alpha, self.s_metric_tests,
                   self.s_composition, self.s_beta, self.s_features, self.s_models, self.s_mags, self.s_ratio,
                   self.s_dispersion, self.s_methods):
            try:
                fn()
            except Exception as e:                       # a broken section must not cost the report
                import traceback
                traceback.print_exc()
                self.add("error", "Section %s failed" % fn.__name__, ("text", "%s: %s" % (type(e).__name__, e)))

    # ------------------------------------------------------------------ PDF (flow layout: items stacked on A4 pages)
    WRAP = {"finding": 120, "reasons": 70, "info": 80, "value": 95, "question": 30, "formula": 45, "name": 40,
            "columns": 60}

    def _cells(self, df):
        cols = list(map(str, df.columns))
        rows = []
        for r in df.itertuples(index=False):
            row = []
            for c, v in zip(cols, r):
                t = fmtv(v)
                w = self.WRAP.get(c)
                row.append(textwrap.fill(t, w) if w else textwrap.shorten(t, 48, placeholder="…"))
            rows.append(row)
        return cols, rows

    LINE = 0.0118              # one text line of a table cell, fraction of the A4 page height

    def _row_h(self, row):
        return self.LINE * max(c.count("\n") + 1 for c in row) + 0.004

    def new_page(self, pdf, title, explain):
        if self.page is not None:
            pdf.savefig(self.page)
            plt.close(self.page)
        self.page = plt.figure(figsize=A4)
        self.header(self.page, title, explain)
        n = len(textwrap.wrap(explain, 120)) + 1 if explain else 0
        self.y = 0.93 - 0.012 * n

    def put_table(self, pdf, title, caption, df):
        cols, rows = self._cells(df)
        if not rows:
            self.put_text(pdf, title, (caption + ": " if caption else "") + "none")
            return
        widths = np.array([max([len(c) * 0.8] + [max(len(x) for x in r[i].split("\n")) for r in rows])
                           for i, c in enumerate(cols)], float)
        widths = np.clip(widths, 3, 120)
        tot = widths.sum()
        fs = 6.8 if tot < 150 else 6.0 if tot < 190 else 5.2
        hdr = self.LINE * 2 + 0.006
        i = 0
        first = True
        while i < len(rows):
            need = (0.02 if caption else 0) + hdr + self._row_h(rows[i])
            if self.y - need < 0.04:
                self.new_page(pdf, title, None)
            if caption:
                self.page.text(0.04, self.y, caption + ("" if first else " (continued)"), fontsize=8, weight="bold",
                               va="top")
                self.y -= 0.02
            j, h = i, hdr
            while j < len(rows) and self.y - h - self._row_h(rows[j]) >= 0.04:
                h += self._row_h(rows[j])
                j += 1
            j = max(j, i + 1)
            ax = self.page.add_axes([0.04, self.y - h, 0.92, h])
            ax.axis("off")
            t = ax.table(cellText=rows[i:j], colLabels=[textwrap.fill(c, 16) for c in cols], loc="upper left",
                         colWidths=list(widths / tot), cellLoc="left", bbox=[0, 0, 1, 1])
            t.auto_set_font_size(False)
            t.set_fontsize(fs)
            heights = [hdr] + [self._row_h(r) for r in rows[i:j]]
            for (r, c), cell in t.get_celld().items():
                cell.set_edgecolor("#d9dde2")
                cell.set_height(heights[r] / h)
                cell.PAD = 0.02
                if r == 0:
                    cell.set_text_props(weight="bold")
                    cell.set_facecolor("#eef2f6")
            self.y -= h + 0.02
            i, first = j, False

    def put_fig(self, pdf, title, f):
        w, h = f.get_size_inches()
        hh = h / A4[1] * (0.94 * A4[0] / w)
        if self.y - hh < 0.03:
            self.new_page(pdf, title, None)
            hh = min(hh, self.y - 0.03)
        buf = io.BytesIO()
        f.savefig(buf, format="png", dpi=200)
        buf.seek(0)
        img = plt.imread(buf)
        ax = self.page.add_axes([0.03, self.y - hh, 0.94, hh])
        ax.imshow(img)
        ax.axis("off")
        self.y -= hh + 0.015

    def put_text(self, pdf, title, txt):
        lines = textwrap.wrap(txt, 125)
        if self.y - 0.012 * len(lines) < 0.04:
            self.new_page(pdf, title, None)
        self.page.text(0.04, self.y, "\n".join(lines), fontsize=8, va="top")
        self.y -= 0.0125 * len(lines) + 0.015

    def header(self, fig, title, explain, caption=None):
        fig.text(0.04, 0.975, title, fontsize=13, weight="bold", va="top")
        fig.text(0.96, 0.975, self.st["name"], fontsize=7, color="#7b8794", va="top", ha="right")
        if explain:
            fig.text(0.04, 0.95, textwrap.fill(explain, 120), fontsize=7, color="#3e4c59", va="top")

    def write_pdf(self, path):
        self.page = None
        with PdfPages(path) as pdf:
            fig = plt.figure(figsize=A4)
            fig.text(0.08, 0.8, self.st["name"], fontsize=20, weight="bold")
            fig.text(0.08, 0.76, "IMP3 cohort analysis · imp-explorer %s · %s" % (self.st["impexplorer"],
                                                                                  self.st["created"][:10]), fontsize=10)
            y = 0.68
            for i, (key, title, items) in enumerate(self.sections):
                fig.text(0.1, y, "%2d  %s" % (i + 1, title), fontsize=10)
                y -= 0.025
            pdf.savefig(fig)
            plt.close(fig)
            for key, title, items in self.sections:
                self.new_page(pdf, title, EXPLAIN.get(key))
                for it in items:
                    if it[0] == "fig":
                        self.put_fig(pdf, title, it[1])
                    elif it[0] == "table":
                        self.put_table(pdf, title, it[1], it[2])
                    else:
                        self.put_text(pdf, title, it[1])
            if self.page is not None:
                pdf.savefig(self.page)
                plt.close(self.page)

    # ------------------------------------------------------------------ HTML
    def write_html(self, path, max_rows=200):
        out = ["<!doctype html><html><head><meta charset='utf-8'><title>%s · cohort</title>" % esc(self.st["name"]),
               "<style>body{font-family:Inter,Helvetica,Arial,sans-serif;color:#1f2933;max-width:1180px;margin:0 auto;"
               "padding:20px 24px;font-size:14px}h1{font-size:24px}h2{border-bottom:2px solid #2a6f97;padding-bottom:4px;"
               "margin-top:40px}nav a{margin-right:12px;white-space:nowrap}p.ex{color:#3e4c59;max-width:950px}"
               ".tw{overflow-x:auto;max-height:560px;overflow-y:auto;border:1px solid #e4e7eb;margin:6px 0 18px}"
               "table{border-collapse:collapse;font-size:12px}th,td{padding:3px 8px;border-bottom:1px solid #e4e7eb;"
               "text-align:left;white-space:nowrap}th{background:#eef2f6;position:sticky;top:0;cursor:pointer}"
               "td.n{text-align:right;font-variant-numeric:tabular-nums}tr.sig td{font-weight:600}"
               "img{max-width:100%;height:auto}.cap{font-weight:600;margin-top:14px}</style>",
               "<script>function srt(th){var t=th.closest('table'),i=[...th.parentNode.children].indexOf(th),"
               "b=t.tBodies[0],r=[...b.rows],d=th.dataset.d=th.dataset.d==='1'?'-1':'1';r.sort(function(x,y){"
               "var a=x.cells[i].dataset.v??x.cells[i].innerText,c=y.cells[i].dataset.v??y.cells[i].innerText,"
               "p=parseFloat(a),q=parseFloat(c);return (isNaN(p)||isNaN(q)?a.localeCompare(c):p-q)*d});"
               "r.forEach(function(e){b.appendChild(e)})}</script></head><body>",
               "<h1>%s</h1><p>IMP3 cohort analysis · imp-explorer %s · %s</p><nav>" % (
                   esc(self.st["name"]), esc(self.st["impexplorer"]), esc(self.st["created"][:16]))]
        out += ["<a href='#%s'>%s</a>" % (esc(k), esc(t)) for k, t, _ in self.sections]
        out.append("</nav>")
        for key, title, items in self.sections:
            out.append("<h2 id='%s'>%s</h2>" % (esc(key), esc(title)))
            if EXPLAIN.get(key):
                out.append("<p class='ex'>%s</p>" % esc(EXPLAIN[key]))
            for it in items:
                if it[0] == "text":
                    out.append("<p>%s</p>" % esc(it[1]))
                elif it[0] == "fig":
                    buf = io.BytesIO()
                    it[1].savefig(buf, format="png", dpi=110)
                    out.append("<img src='data:image/png;base64,%s'>" % base64.b64encode(buf.getvalue()).decode())
                else:
                    cap, df = it[1], it[2]
                    if cap:
                        out.append("<div class='cap'>%s</div>" % esc(cap))
                    out.append(self.html_table(df, max_rows))
        out.append("</body></html>")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("\n".join(out))

    def html_table(self, df, max_rows):
        n = len(df)
        df = df.head(max_rows)
        qcol = "q" if "q" in df.columns else None
        h = ["<div class='tw'><table><thead><tr>"] + ["<th onclick='srt(this)'>%s</th>" % esc(c) for c in df.columns]
        h.append("</tr></thead><tbody>")
        for r in df.itertuples(index=False):
            d = dict(zip(df.columns, r))
            sig = qcol and isinstance(d[qcol], (float, np.floating)) and d[qcol] < self.q
            h.append("<tr class='sig'>" if sig else "<tr>")
            for v in r:
                if isinstance(v, (int, float, np.integer, np.floating)) and not isinstance(v, bool):
                    h.append("<td class='n' data-v='%s'>%s</td>" % (v, esc(fmtv(v))))
                else:
                    h.append("<td>%s</td>" % esc(fmtv(v)))
            h.append("</tr>")
        h.append("</tbody></table></div>")
        if n > max_rows:
            h.append("<p class='ex'>first %d of %d rows; all rows in the TSV files under stats/</p>" % (max_rows, n))
        return "".join(h)


def write_report(C, A, res, st, out):
    R = Report(C, A, res, st, out)
    R.build()
    R.write_pdf(os.path.join(out, "cohort.report.pdf"))
    R.write_html(os.path.join(out, "cohort.report.html"))
    plt.close("all")
    C.note("report: %s/cohort.report.pdf + .html (%d sections)" % (out, len(R.sections)))
