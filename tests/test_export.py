import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import make_fixture  # noqa: E402
from impexplorer import export  # noqa: E402
from impexplorer.readers import imp3  # noqa: E402
from impexplorer.release import read_manifest, read_table  # noqa: E402


def build(tmp):
    root = make_fixture.main(str(tmp), n=6)
    return root, os.path.join(str(tmp), "metadata.tsv")


def test_sample_readers(tmp_path):
    root, _ = build(tmp_path)
    b = os.path.join(root, "S01")
    q = imp3.qc(b, "S01")
    m = dict(zip(q["omic"] + ":" + q["metric"], q["value"]))
    assert m["mg:trimmed:pairs"] > 0 and 0.9 < m["mg:mapped_fraction"] <= 1
    assert m["mgmt:contigs"] == 400 and m["mgmt:n50"] > 0
    t = imp3.taxonomy(b, "S01")
    assert set(t["source"]) == {"kraken2_bracken", "motus"}
    g = t[(t["omic"] == "mg") & (t["rank"] == "genus")]
    assert abs(g["rel_abundance"].sum() - 1) < 1e-9 and "unclassified" not in set(g["taxon"])
    bins = imp3.bins(b, "S01")
    assert len(bins) >= 6 and bins["genus"].str.len().gt(0).all() and set(bins["quality"]) <= {"high", "medium", "low"}
    assert bins["bin"].str.startswith("bin.").all()          # padded cells stripped
    ab = imp3.bin_abundance(b, "S01")
    assert set(ab["omic"]) == {"mg", "mt"} and ab["share_mapped"].between(0, 1).all()
    f = imp3.functions(b, "S01")
    assert {"kegg_ko", "kegg_module"} == set(f["db"]) and (f["reads"] > 0).all()


def test_new_bracken_layout(tmp_path):
    root, _ = build(tmp_path)
    sid = "S02"                                                # odd samples use Analysis/taxonomy/bracken/<omic>/
    assert os.path.isdir(os.path.join(root, sid, "Analysis", "taxonomy", "bracken"))
    t = imp3.taxonomy(os.path.join(root, sid), sid)
    assert len(t[(t["source"] == "kraken2_bracken") & (t["rank"] == "species") & (t["omic"] == "mg")]) == 10
    mo = t[(t["source"] == "motus") & (t["omic"] == "mg")]          # hierarchical mOTUs: one profile per rank
    assert set(mo["rank"]) == {"genus", "species"} and "unassigned" not in set(mo["taxon"])
    assert mo[mo["rank"] == "species"]["taxid"].str.startswith("ref_mOTU").all()
    assert abs(mo.groupby("rank")["rel_abundance"].sum() - 1).max() < 1e-9


def test_missing_files_give_empty_tables(tmp_path):
    d = tmp_path / "empty"
    d.mkdir()
    r = imp3.read_sample(str(d), "X")
    assert all(len(r[k]) == 0 for k in ("qc", "taxonomy", "bins", "bin_abundance", "functions"))


def test_export_release(tmp_path):
    root, md = build(tmp_path)
    out = tmp_path / "release"
    rc = export.main(["--imp-root", root, "--metadata", md, "--out", str(out), "--study", "fixture",
                      "--threads", "2", "--format", "tsv"])
    assert rc == 0
    s = read_table(str(out), "samples")
    assert len(s) == 6 and s["has_mg"].all() and "group" in s.columns
    man = read_manifest(str(out))
    assert man["study"] == "fixture" and man["tables"]["bins"]["rows"] > 30
    tax = read_table(str(out), "taxonomy")
    assert tax["sample_id"].nunique() == 6
    ab = read_table(str(out), "bin_abundance")
    assert set(ab["bin_id"]) <= set(read_table(str(out), "bins")["bin_id"])


def test_app_builds(tmp_path):
    """the app starts on a fixture release and every page renders its layout (skipped without dash)"""
    try:
        import dash  # noqa: F401
        import duckdb  # noqa: F401
    except ImportError:
        return
    from impexplorer.app import create_app
    root, md = build(tmp_path)
    out = tmp_path / "release"
    export.main(["--imp-root", root, "--metadata", md, "--out", str(out), "--threads", "1"])
    app = create_app(str(out))
    assert [p.name for p in app._imp_pages] == ["Overview", "QC", "Taxonomy", "MAGs", "Functions", "Search"]
    for p in app._imp_pages:
        assert p.layout(app._imp_data) is not None
    assert "group" in app._imp_data.group_cols


# ----------------------------------------------------------------------------- sample ids and paths
def test_derive_ids():
    names = ["run_s01", "run_s02", "other"]
    assert export.derive_ids(names[:2], id_regex=r"run_(s\d+)", id_case="upper") == ["S01", "S02"]
    assert export.derive_ids(names, id_map={"other": "S03"}, id_regex=r"run_(s\d+)") == ["s01", "s02", "S03"]
    assert export.derive_ids(["a", "b"], id_case="keep") == ["a", "b"]
    assert export.derive_ids(["ab_1"], id_regex=r"\d+") == ["1"]                  # no group: whole match
    for bad in (dict(names=["x", "X"], id_case="upper"),                           # duplicates after translation
                dict(names=["run_s01", "zzz"], id_regex=r"run_(s\d+)")):          # regex does not match
        try:
            export.derive_ids(**bad)
        except ValueError:
            continue
        raise AssertionError("expected ValueError for %s" % bad)


def test_id_map_file(tmp_path):
    p = tmp_path / "map.tsv"
    p.write_text("folder\tsample_id\nrun_a\tA\nrun_b\tB\n")
    assert export.read_id_map(str(p)) == {"run_a": "A", "run_b": "B"}
    q = tmp_path / "map.csv"
    q.write_text("run_a,A\n")
    assert export.read_id_map(str(q)) == {"run_a": "A"}


def test_export_translated_ids_and_paths(tmp_path):
    """folders run_s01 ... + metadata ids S01 ...: --id-regex/--id-case translate, all samples get metadata"""
    make_fixture.main(str(tmp_path), n=4, folder_fmt="run_{sid_lower}")
    root, md = str(tmp_path / "imp"), str(tmp_path / "metadata.tsv")
    assert sorted(os.listdir(root)) == ["run_s01", "run_s02", "run_s03", "run_s04"]
    out = tmp_path / "rel"
    rc = export.main(["--imp-root", root, "--id-regex", r"run_(\w+)", "--id-case", "upper", "--metadata", md,
                      "--out", str(out), "--threads", "1", "--format", "tsv"])
    assert rc == 0
    s = read_table(str(out), "samples")
    assert list(s["sample_id"]) == ["S01", "S02", "S03", "S04"] and s["group"].notna().all()
    assert list(s["imp_dir"]) == ["run_s01", "run_s02", "run_s03", "run_s04"]       # relative to --imp-root
    assert read_manifest(str(out))["metadata_without_output"] == []
    assert read_table(str(out), "bins")["bin_id"].str.startswith("S0").all()
    # without translation the metadata does not match, and the manifest says so
    out2 = tmp_path / "rel2"
    export.main(["--imp-root", root, "--metadata", md, "--out", str(out2), "--threads", "1", "--format", "tsv",
                 "--paths", "none"])
    assert read_manifest(str(out2))["metadata_without_output"] == ["S01", "S02", "S03", "S04"]
    assert read_table(str(out2), "samples")["imp_dir"].fillna("").eq("").all()
    # an id map does the same as the regex
    m = tmp_path / "map.tsv"
    m.write_text("".join("run_s0%d\tS0%d\n" % (i, i) for i in range(1, 5)))
    out3 = tmp_path / "rel3"
    export.main(["--imp-root", root, "--id-map", str(m), "--metadata", md, "--out", str(out3), "--threads", "1",
                 "--format", "tsv", "--paths", "absolute"])
    s3 = read_table(str(out3), "samples")
    assert list(s3["sample_id"]) == ["S01", "S02", "S03", "S04"] and s3["imp_dir"].map(os.path.isabs).all()


def test_samples_sheet_relative_paths(tmp_path):
    root = make_fixture.main(str(tmp_path), n=3)
    sheet = tmp_path / "samples.tsv"
    sheet.write_text("sample_id\timp_dir\tgroup\n" + "".join("S0%d\timp/S0%d\tg%d\n" % (i, i, i % 2)
                                                              for i in range(1, 4)))
    out = tmp_path / "rel"
    assert export.main(["--samples", str(sheet), "--out", str(out), "--threads", "1", "--format", "tsv"]) == 0
    s = read_table(str(out), "samples")
    assert list(s["imp_dir"]) == ["imp/S01", "imp/S02", "imp/S03"] and s["has_mt"].all()


# ----------------------------------------------------------------------------- cohort analysis
def test_config_settings(tmp_path):
    """one settings file per study: 'key: value' (written by --write-config) or bash COHORT_* (sRNA-IMP config)"""
    from impexplorer.cohort import settings
    sh = tmp_path / "srna.config.sh"
    sh.write_text('THREADS=8\nCOHORT_GROUP_COL="${COHORT_GROUP_COL-group}"   # comment\n'
                  'COHORT_GROUP_ORDER="${COHORT_GROUP_ORDER-control,case_a,case_b}"\n'
                  'COHORT_BLOCK="${COHORT_BLOCK-${COHORT_MATCH_COL-site}}"\n'
                  'COHORT_PERMUTATIONS="${COHORT_PERMUTATIONS-99}"\nCOHORT_MIN_MIRNA_READS=50\n'
                  'COHORT_MODEL="${COHORT_MODEL-group+site, age+sex | group=case_a}"\n'
                  'COHORT_MODEL_NOTES="${COHORT_MODEL_NOTES-M2=age in case_a}"\n')
    a = settings.parse(["--data", str(tmp_path), "--config", str(sh), "--permutations", "19"])
    assert a.group_order == "control,case_a,case_b" and a.block == "site"
    assert a.model == "group+site, age+sex | group=case_a" and a.model_notes == "M2=age in case_a"
    assert a.permutations == 19                     # the command line wins over the file
    a = settings.parse(["--data", str(tmp_path), "--srna-config", str(sh)])      # old option name still works
    assert a.permutations == 99
    # template with the values given, read back unchanged
    t = tmp_path / "cohort.yaml"
    settings.parse(["--write-config", str(t), "--block", "centre", "--model", "group+centre, group+centre+age",
                    "--prevalence", "0.2"])
    b = settings.parse(["--data", str(tmp_path), "--config", str(t)])
    assert b.block == "centre" and b.model == "group+centre, group+centre+age" and b.prevalence == 0.2
    assert b.group_col == "group" and b.permutations == 999 and b.qc_drop == "FAIL"
    # hand-written file, mixed key styles
    y = tmp_path / "study.yaml"
    y.write_text("# study design\nGroup_Col: diagnosis\ncohort_reference: HC\nblock: 'site'   # quoted\n"
                 "model: diagnosis+site, diagnosis+site+sex  # two models\nsig-q: 0.05\n")
    c = settings.parse(["--data", str(tmp_path), "--config", str(y)])
    assert (c.group_col, c.reference, c.block, c.sig_q) == ("diagnosis", "HC", "site", 0.05)
    assert c.model == "diagnosis+site, diagnosis+site+sex"


def test_permanova_same_as_pairwise_loop():
    import numpy as np
    from impexplorer.cohort.stats import bray_curtis, permanova
    rng = np.random.default_rng(2)
    X = rng.random((24, 8))
    X[:12, 0] += 2                                  # a real difference
    D = bray_curtis(X)
    lab = np.array(["a"] * 12 + ["b"] * 12)
    F, r2, p = permanova(D, lab, 199)
    # reference: pairwise sums of squares (Anderson 2001)
    n, D2 = len(lab), D ** 2
    sst = D2[np.triu_indices(n, 1)].sum() / n
    ssw = sum(D2[np.ix_(lab == g, lab == g)][np.triu_indices(12, 1)].sum() / 12 for g in "ab")
    assert abs(F - (sst - ssw) / (ssw / (n - 2))) < 1e-9 and abs(r2 - (sst - ssw) / sst) < 1e-12 and p < 0.05


def test_ols_hc3_estimates():
    import numpy as np
    from impexplorer.cohort.stats import ols_hc3
    rng = np.random.default_rng(1)
    X = np.column_stack([np.ones(50), rng.normal(size=50)])
    Y = (X @ np.array([1.0, 2.0]))[:, None] + rng.normal(scale=0.1, size=(50, 1))
    B, SE, T, P, df = ols_hc3(X, Y)
    assert np.allclose(B[:, 0], np.linalg.lstsq(X, Y, rcond=None)[0][:, 0]) and df == 48 and P[1, 0] < 1e-10


def test_cohort_on_fixture(tmp_path):
    """imp-explorer cohort on a 30-sample release: tables, report, cohort tables in the release, app page"""
    from impexplorer import cohort
    root = make_fixture.main(str(tmp_path), n=30)
    out = tmp_path / "release"
    export.main(["--imp-root", root, "--metadata", str(tmp_path / "metadata.tsv"), "--out", str(out),
                 "--threads", "2"])
    try:
        import matplotlib  # noqa: F401
        report = []
    except ImportError:
        report = ["--no-report"]
    rc = cohort.main(["--data", str(out), "--block", "site", "--group-order", "control,case_a,case_b",
                      "--model", "group+site, group+site+sex+age, age+sex | group=case_a",
                      "--model-notes", "M3=age within case_a", "--permutations", "49", "--nmds", "kraken.genus.mg",
                      "--layers", "kraken.genus.mg,kraken.genus.mt,mag.genus.mg,mag.genus.mt,kegg_module.mg,"
                      "kegg_module.mt,ratio.kegg_module,ratio.mag.genus"] + report)
    assert rc == 0
    c = out / "cohort"
    mo = pd.read_csv(c / "stats" / "cohort.models.tsv", sep="\t")
    assert set(mo["model"]) >= {"M1", "M2"} and mo["focus"].any()
    assert {"group: case_a vs control", "group: case_b vs control"} <= set(mo["term"])
    ov = pd.read_csv(c / "stats" / "cohort.models_overview.tsv", sep="\t")
    assert (ov[ov["model"] == "M3"]["features"].fillna(0) == 0).all()      # 10 case_a samples: too few, skipped
    ft = pd.read_csv(c / "stats" / "cohort.feature_tests.tsv", sep="\t")
    assert {"case_a vs control (site_1, site_2)", "control: site_2 vs site_1"} <= set(ft["comparison"])
    pm = pd.read_csv(c / "stats" / "cohort.permanova.tsv", sep="\t")
    assert "groups | site (partial)" in set(pm["comparison"])
    assert (c / "tables" / "cohort.ratio.kegg_module.log2ratio.tsv").exists()
    if not report:
        assert (c / "cohort.report.pdf").stat().st_size > 10000 and (c / "cohort.report.html").exists()
    man = read_manifest(str(out))
    assert man["cohort"]["reference"] == "control" and man["cohort"]["models"][2]["note"] == "age within case_a"
    assert len(read_table(str(out), "cohort_models")) == len(mo)
    try:
        import dash  # noqa: F401
    except ImportError:
        return
    from impexplorer.app import create_app
    app = create_app(str(out))
    assert "Cohort" in [p.name for p in app._imp_pages]
    page = next(p for p in app._imp_pages if p.name == "Cohort")
    assert page.layout(app._imp_data) is not None


def test_clustered_heatmap_modes():
    try:
        import plotly  # noqa: F401
        import scipy  # noqa: F401
        import dash_ag_grid  # noqa: F401
    except ImportError:
        return
    import numpy as np
    from impexplorer.app.pages._common import clustered_heatmap
    rng = np.random.default_rng(0)
    M = pd.DataFrame(rng.random((24, 30)), index=["S%02d" % i for i in range(24)],
                     columns=["K%05d" % i for i in range(30)])
    M.iloc[:, 3] = 0                                    # a constant feature must not break the clustering
    g = pd.Series(["a", "b", "c"] * 8, index=M.index)
    site = {"site": pd.Series(["x"] * 12 + ["y"] * 12, index=M.index)}
    for rows in ("variance", "cluster"):
        for cols in ("top", "within", "all"):
            for values in ("log", "clr"):
                fig = clustered_heatmap(M, g, {"a": "#111", "b": "#222", "c": "#333"}, rows, cols, values, top=20,
                                        annot=site)
                heat = [t for t in fig.data if t.type == "heatmap"][-1]
                assert np.asarray(heat.z).shape == (20, 24)


# ----------------------------------------------------------------------------- sRNA layer
def make_srna_cohort(root, samples):
    """minimal Cohort folder as written by sRNA-IMP (srna_cohort.py), two fractions"""
    import numpy as np
    rng = np.random.default_rng(5)
    t, s = os.path.join(root, "tables"), os.path.join(root, "stats")
    os.makedirs(t)
    os.makedirs(s)
    groups = ["control", "case_a", "case_b"]
    g = {x: groups[i % 3] for i, x in enumerate(samples)}
    site = {x: "site_%d" % (1 + (i // 3) % 2) for i, x in enumerate(samples)}
    rows = []
    for ty in ("sRNA-Ext", "sRNA-Int"):
        for x in samples:
            rows.append([x, ty, g[x], "PASS", "yes", site[x], "F", 60, int(rng.integers(5e6, 9e6)),
                         float(rng.random() * 0.1)])
    pd.DataFrame(rows, columns=["sample", "type", "group", "qc", "in_statistics", "site", "sex", "age",
                                "raw_reads", "human_fraction"]).to_csv(os.path.join(t, "cohort.sample_stats.tsv"),
                                                                       sep="\t", index=False)
    pd.DataFrame({"sample": samples, "group": [g[x] for x in samples], "site": [site[x] for x in samples],
                  "sex": "F", "age": 60}).to_csv(os.path.join(t, "cohort.metadata.tsv"), sep="\t", index=False)
    pd.DataFrame([[x, ty, g[x], "PASS", "yes", "-", "-"] for ty in ("sRNA-Ext", "sRNA-Int") for x in samples],
                 columns=["sample", "type", "group", "qc", "in_statistics", "rules", "details"]).to_csv(
        os.path.join(t, "cohort.qc_flags.tsv"), sep="\t", index=False)
    feats = {"genus": ["Bacteroides", "Faecalibacterium", "Blautia_A"], "mirna": ["hsa-miR-21-5p", "hsa-miR-192-5p"],
             "rna_class": ["rRNA", "tRNA", "unknown"], "fnas_kegg_ko": ["K00001"]}
    for ty in ("sRNA-Ext", "sRNA-Int"):
        for st, fs in feats.items():
            m = pd.DataFrame(rng.integers(0, 500, (len(fs), len(samples))), index=fs, columns=samples)
            m.index.name = "feature"
            m.to_csv(os.path.join(t, "cohort.%s.%s.counts.tsv" % (ty, st)), sep="\t")
        pd.DataFrame({"feature": feats["genus"], "prevalence": 1.0, "H": [9.1, 1.0, 0.5], "p": [0.01, 0.6, 0.8],
                      "q": [0.03, 0.9, 0.9], "mean_rel_case_a": 0.2, "mean_rel_case_b": 0.3,
                      "mean_rel_control": 0.4}).to_csv(os.path.join(s, "cohort.%s.genus.by_group.tsv" % ty),
                                                       sep="\t", index=False)
        pd.DataFrame({"metric": ["raw_reads"], "n": [len(samples)], "H": [1.2], "p": [0.5], "q": [0.5],
                      "median_case_a": [6e6], "median_case_b": [7e6], "median_control": [6.5e6]}).to_csv(
            os.path.join(s, "cohort.%s.metrics_by_group.tsv" % ty), sep="\t", index=False)
        for tag, lab in (("case_a_vs_control", "case_a vs control (within site: site_1)"),
                         ("control_site_1_vs_control_site_2", "control site_1 vs control site_2 (site effect "
                                                              "within control)")):
            pd.DataFrame({"comparison": [lab], "n_a": [5], "n_b": [5], "pseudo_F": [1.5], "R2": [0.1],
                          "p": [0.2]}).to_csv(os.path.join(s, "cohort.%s.matched.%s.permanova.tsv" % (ty, tag)),
                                              sep="\t", index=False)
            pd.DataFrame({"feature": feats["genus"], "prevalence": 1.0, "mean_rel_a": 0.3, "mean_rel_b": 0.2,
                          "P(a>b)": [0.8, 0.5, 0.4], "p": [0.02, 0.9, 0.5], "q": [0.06, 0.9, 0.9]}).to_csv(
                os.path.join(s, "cohort.%s.matched.%s.genus.tsv" % (ty, tag)), sep="\t", index=False)
        pd.DataFrame({"variable": ["group", "site"], "n": len(samples), "levels": [3, 2], "pseudo_F": [1.1, 2.0],
                      "R2": [0.05, 0.04], "p": [0.3, 0.01]}).to_csv(
            os.path.join(s, "cohort.%s.genus.permanova.tsv" % ty), sep="\t", index=False)
        pd.DataFrame({"sample": samples, "group": [g[x] for x in samples], "NMDS1": rng.random(len(samples)),
                      "NMDS2": rng.random(len(samples)), "stress": 0.15, "permanova_R2": 0.05,
                      "permanova_p": 0.3}).to_csv(os.path.join(s, "cohort.%s.genus.nmds.tsv" % ty), sep="\t",
                                                  index=False)
        for m in ("M1", "M2"):
            pd.DataFrame({"feature": feats["genus"] * 2, "name": "", "n": len(samples),
                          "term": ["group: case_a vs control"] * 3 + ["site: site_2 vs site_1"] * 3,
                          "estimate": [1.0, 0.1, -0.2, 0.3, 0.0, 0.1], "se": 0.2, "t": 1.0,
                          "p": [0.001, 0.5, 0.4, 0.2, 0.9, 0.6], "q": [0.003, 0.7, 0.7, 0.6, 0.9, 0.9]}).to_csv(
                os.path.join(s, "cohort.%s.model_%s.genus.tsv" % (ty, m)), sep="\t", index=False)
    pd.DataFrame([["M1", "group+site", ty, len(samples), "-"] for ty in ("sRNA-Ext", "sRNA-Int")] +
                 [["M2", "group+site+age", ty, len(samples), "dose"] for ty in ("sRNA-Ext", "sRNA-Int")],
                 columns=["model", "formula", "type", "libraries", "notes"]).to_csv(
        os.path.join(s, "cohort.models.tsv"), sep="\t", index=False)
    pd.DataFrame({"type": ["sRNA-Ext"], "set": ["genus"], "comparison": ["case_a vs control (site_1)"],
                  "group": ["case_a"], "block": ["site_1"], "n_cases": [5], "n_controls": [5],
                  "median_dist_cases": [0.4], "median_dist_controls": [0.3], "ratio": [1.3], "p": [0.2],
                  "q": [0.2]}).to_csv(os.path.join(s, "cohort.heterogeneity.dispersion.tsv"), sep="\t", index=False)
    return root


def test_srna_layer(tmp_path):
    """imp-explorer srna: sRNA-IMP cohort results into the release, Cohort / sRNA / Search pages"""
    from impexplorer import srna
    root = make_fixture.main(str(tmp_path), n=12)
    out = tmp_path / "release"
    export.main(["--imp-root", root, "--metadata", str(tmp_path / "metadata.tsv"), "--out", str(out),
                 "--threads", "2"])
    samples = ["s%02d" % i for i in range(1, 13)]          # sRNA-IMP names in lower case
    make_srna_cohort(str(tmp_path / "Cohort"), samples)
    assert srna.main(["--data", str(out), "--srna-cohort", str(tmp_path / "Cohort"), "--id-case", "upper"]) == 0
    man = read_manifest(str(out))["srna"]
    assert man["samples_in_release"] == 12 and set(man["fractions"]) == {"srna_ext", "srna_int"}
    f = read_table(str(out), "srna_features")
    assert set(f["set"]) == {"genus", "mirna", "rna_class"}                  # fnas_* left out
    assert abs(f.groupby(["sample_id", "omic", "set"])["rel"].sum() - 1).max() < 1e-9
    t = read_table(str(out), "srna_tests")
    assert {"all groups", "case_a vs control (site_1)", "control: site_1 vs site_2"} <= set(t["comparison"])
    assert {"srna_ext.genus", "srna_ext.metrics"} <= set(t["layer"])
    mo = read_table(str(out), "srna_models")
    assert mo[mo["focus"].astype(str).isin(["True", "true"])]["term"].unique().tolist() == ["group: case_a vs control"]
    ov = read_table(str(out), "srna_models_overview")
    assert ov[(ov["model"] == "M1") & (ov["layer"] == "srna_ext.genus")]["features_q_focus"].iloc[0] == 1
    pm = read_table(str(out), "srna_permanova")
    assert "control: site_1 vs site_2" in set(pm["comparison"]) and "site (all samples)" in set(pm["comparison"])
    try:
        import dash  # noqa: F401
    except ImportError:
        return
    from impexplorer.app import create_app
    from impexplorer.app.pages import cohort, search
    app = create_app(str(out))
    names = [p.name for p in app._imp_pages]
    assert "sRNA" in names and "Cohort" in names            # Cohort page with the sRNA results alone
    d = app._imp_data
    for p in app._imp_pages:
        assert p.layout(d) is not None
    v, _ = cohort.feature_values(d, "srna_int.genus", "Bacteroides")
    assert len(v) == 12 and (v > 0).all()
    from impexplorer.app.pages import overview
    assert overview.srna_tiles(d) is not None and "sRNA-Ext" in overview.sample_table(d).columns
    hits = search.find(d, "miR-21")
    assert (hits["layer"] == "sRNA").any()
    ser = search.series_for(d, "taxonomy", "kraken2_bracken · mg · genus", "Bacteroides", d.samples)
    assert "sRNA-Ext" in ser and "sRNA-Int" in ser


def make_metab_tables(folder, meta):
    """generic metabolomics tables (metabolomics_export.py format), sample ids as 'EB01_s01_t0' to test the mapping;
    feature 'butyrate' higher in case_a"""
    import numpy as np
    os.makedirs(folder, exist_ok=True)
    rng = np.random.default_rng(3)
    pd.DataFrame([{"assay": "scfa", "label": "Short-chain fatty acids", "platform": "GC-MS", "matrix": "faecal water",
                   "acquisition_date": "2023-11-06", "value_type": "concentration", "units": "umol/g"},
                  {"assay": "lipids", "label": "Lipids", "platform": "LC-MS", "matrix": "faecal extract",
                   "acquisition_date": "2024-09-26", "value_type": "normalised area", "units": "relative"}]
                 ).to_parquet(os.path.join(folder, "metabolomics_assays.parquet"), index=False)
    feats = [("scfa", "acetate", "acetate", "identified (targeted)"), ("scfa", "butyrate", "butyrate", "identified (targeted)"),
             ("scfa", "propionate", "propionate", "identified (targeted)")] + \
            [("lipids", "L%02d@%.1f" % (i, 5 + i / 3), "PC(3%d:1)" % i, "identified" if i < 6 else "unknown")
             for i in range(10)]
    pd.DataFrame(feats, columns=["assay", "feature_id", "name", "annotation"]).to_csv(
        os.path.join(folder, "metabolomics_features.tsv"), sep="\t", index=False)
    rows = []
    for _, r in meta.iterrows():
        sid = "EB01_%s_t0" % r["sample_id"].lower()
        for a, f, _, _ in feats:
            if a == "lipids" and r["sample_id"] == "S03":
                continue                                       # one sample not measured in one assay
            v = rng.lognormal(3 if a == "scfa" else 0, 0.4)
            if f == "butyrate" and r["group"] == "case_a":
                v *= 4
            if f == "L09@8.0" and rng.random() < 0.5:
                v = 0.0                                        # half not detected
            rows.append((sid, a, f, v))
    pd.DataFrame(rows, columns=["sample_id", "assay", "feature_id", "value"]).to_csv(
        os.path.join(folder, "metabolomics_values.tsv.gz"), sep="\t", index=False)


def test_metabolomics_layer(tmp_path):
    """imp-explorer metabolomics: generic tables into the release; cohort layers metab.*; app pages"""
    from impexplorer import cohort, metabolomics
    root = make_fixture.main(str(tmp_path), n=30)
    out = tmp_path / "release"
    export.main(["--imp-root", root, "--metadata", str(tmp_path / "metadata.tsv"), "--out", str(out),
                 "--threads", "2"])
    meta = pd.read_csv(tmp_path / "metadata.tsv", sep="\t")
    make_metab_tables(str(tmp_path / "metab"), meta)
    assert metabolomics.main(["--data", str(out), "--metab", str(tmp_path / "metab"),
                              "--id-regex", r"_(s\d+)_", "--id-case", "upper"]) == 0
    man = read_manifest(str(out))["metabolomics"]
    assert man["samples_in_release"] == 30 and set(man["assays"]) == {"scfa", "lipids"}
    assert man["assays"]["lipids"]["n_samples"] == 29
    v = read_table(str(out), "metab_values")
    assert set(v["sample_id"]) == set(meta["sample_id"]) and v["source_id"].str.startswith("EB01_").all()
    rc = cohort.main(["--data", str(out), "--block", "site", "--group-order", "control,case_a,case_b",
                      "--model", "group+site", "--permutations", "49", "--no-report",
                      "--layers", "kraken.genus.mg,metab"])
    assert rc == 0
    t = read_table(str(out), "cohort_tests")
    assert {"metab.scfa", "metab.lipids"} <= set(t["layer"])
    b = t[(t["layer"] == "metab.scfa") & (t["feature"] == "butyrate") & t["comparison"].str.startswith("case_a vs")]
    assert len(b) and b["effect"].iloc[0] > 0.8 and b["p"].iloc[0] < 0.01
    assert set(t[t["layer"] == "metab.lipids"]["name"]) >= {"PC(30:1)", "PC(39:1)"}
    mo = read_table(str(out), "cohort_models")
    m = mo[(mo["layer"] == "metab.scfa") & (mo["feature"] == "butyrate") & (mo["term"] == "group: case_a vs control")]
    assert len(m) and 1.3 < m["estimate"].iloc[0] < 2.7                # log2(4) = 2
    pm = read_table(str(out), "cohort_permanova")
    assert "metab.lipids" in set(pm["layer"])
    st = read_manifest(str(out))["cohort"]["layers"]
    assert st["metab.scfa"]["label"] == "Metabolomics: Short-chain fatty acids" and st["metab.lipids"]["samples"] == 29
    try:
        import dash  # noqa: F401
    except ImportError:
        return
    from impexplorer.app import create_app
    from impexplorer.app.pages import cohort as cpage, overview, search
    from impexplorer.app.pages import metabolomics as mpage
    app = create_app(str(out))
    assert "Metabolomics" in [p.name for p in app._imp_pages]
    d = app._imp_data
    for p in app._imp_pages:
        assert p.layout(d) is not None
    x, unit = cpage.feature_values(d, "metab.lipids", "L01@5.3")
    assert len(x) == 29 and "normalised area" in unit
    assert cpage.feature_box(d, "metab.scfa", "butyrate") is not None
    hits = search.find(d, "butyr")
    assert (hits["layer"] == "metabolomics").any()
    h = hits[hits["layer"] == "metabolomics"].iloc[0]
    ser = search.series_for(d, h["layer"], h["kind"], h["name"], d.samples, h["fid"])
    assert len(ser["scfa"]) == 30
    assert overview.metab_tiles(d) is not None and "metab scfa" in overview.sample_table(d).columns
    assert len(mpage.features(d, "lipids", "known")) == 6
    M = mpage.profile_matrix(d, "kraken.genus.mg", list(d.samples["sample_id"]))
    r = mpage.spearman_all(M.iloc[:, 0], M)
    assert abs(r.iloc[0] - 1) < 1e-9


def make_prot_tables(folder, meta):
    """generic proteomics tables (proteomics_export.py format): 40 proteins of 4 genera with KO / module; fractions
    extra (S05 missing) and intra; Akkermansia proteins higher in case_a"""
    import numpy as np
    os.makedirs(folder, exist_ok=True)
    rng = np.random.default_rng(5)
    gen = [("Bacteroides", "Bacteroidota"), ("Blautia_A", "Bacillota_A"), ("Akkermansia", "Verrucomicrobiota"),
           ("Faecalibacterium", "Bacillota_A")]
    prots = []
    for i in range(40):
        g, ph = gen[i % 4]
        prots.append({"protein_id": "expobiome_imp_%d" % i, "gene_id": "G%d" % i, "contig_id": "c%d" % i,
                      "assembly_sample": "S01", "taxid": "", "rank": "genus", "taxon": g,
                      "lineage": "d_Bacteria;p_%s;g_%s" % (ph, g), "product": "protein %d" % i, "gene": "",
                      "kegg_ko": "K%05d" % (i % 10) + (",K09999" if i == 3 else ""), "kegg_module": "M%05d" % (i % 3),
                      "kegg_pathway": "ko00010", "enzyme_ec": "", "cog": "", "pfam": "", "cazy": "", "go": "",
                      "tcdb": ""})
    P = pd.DataFrame(prots)
    P.to_csv(os.path.join(folder, "proteomics_proteins.tsv.gz"), sep="\t", index=False)
    rows = []
    for _, r in meta.iterrows():
        for fr in ("extra", "intra"):
            if fr == "extra" and r["sample_id"] == "S05":
                continue
            for _, q in P.iterrows():
                if rng.random() < 0.2:
                    continue
                v = rng.lognormal(16, 0.5) * (5 if q["taxon"] == "Akkermansia" and r["group"] == "case_a" else 1)
                rows.append((r["sample_id"].lower(), fr, q["protein_id"], v, int(rng.integers(1, 6))))
    A = pd.DataFrame(rows, columns=["sample_id", "fraction", "protein_id", "auc", "n_peptides"])
    A.to_csv(os.path.join(folder, "proteomics_abundance.tsv.gz"), sep="\t", index=False)
    tot = A.groupby(["sample_id", "fraction"])["auc"].sum() * 1.1          # 10 % of the AUC unannotated
    S = A.groupby(["sample_id", "fraction"]).agg(n_proteins=("protein_id", "size"), n_peptides=("n_peptides", "sum"))
    S["total_auc"] = tot
    S.reset_index().to_csv(os.path.join(folder, "proteomics_samples.tsv"), sep="\t", index=False)
    J = A.merge(P, on="protein_id")
    T = []
    for rank, col in (("genus", "taxon"), ("phylum", None)):
        x = J.assign(t=J["lineage"].str.extract(r"p_([^;]+)")[0] if col is None else J[col])
        x = x.groupby(["sample_id", "fraction", "t"], as_index=False).agg(auc=("auc", "sum"), n_proteins=("auc", "size"))
        T.append(x.rename(columns={"t": "taxon"}).assign(rank=rank))
    T = pd.concat(T)
    T["rel"] = T["auc"] / T.set_index(["sample_id", "fraction"]).index.map(tot).values
    T.to_csv(os.path.join(folder, "proteomics_taxonomy.tsv.gz"), sep="\t", index=False)
    F = []
    for db in ("kegg_ko", "kegg_module"):
        x = J.assign(feature=J[db].str.split(",")).explode("feature")
        x = x.groupby(["sample_id", "fraction", "feature"], as_index=False).agg(auc=("auc", "sum"),
                                                                               n_proteins=("auc", "size"))
        F.append(x.assign(db=db))
    F = pd.concat(F)
    F["rel"] = F["auc"] / F.set_index(["sample_id", "fraction"]).index.map(tot).values
    F.to_csv(os.path.join(folder, "proteomics_functions.tsv.gz"), sep="\t", index=False)


def test_proteomics_layer(tmp_path):
    """imp-explorer proteomics: protein profiles next to mg / mt, cohort layers prot.*, Proteomics page"""
    from impexplorer import cohort, proteomics
    root = make_fixture.main(str(tmp_path), n=30)
    out = tmp_path / "release"
    export.main(["--imp-root", root, "--metadata", str(tmp_path / "metadata.tsv"), "--out", str(out),
                 "--threads", "2"])
    meta = pd.read_csv(tmp_path / "metadata.tsv", sep="\t")
    make_prot_tables(str(tmp_path / "prot"), meta)
    assert proteomics.main(["--data", str(out), "--prot", str(tmp_path / "prot"), "--id-case", "upper"]) == 0
    man = read_manifest(str(out))["proteomics"]
    assert man["samples_in_release"] == 30 and set(man["fractions"]) == {"mp_extra", "mp_intra"}
    s = read_table(str(out), "prot_samples")
    assert len(s) == 59 and set(s["omic"]) == {"mp_extra", "mp_intra"}
    rc = cohort.main(["--data", str(out), "--block", "site", "--group-order", "control,case_a,case_b",
                      "--model", "group+site", "--permutations", "49", "--no-report",
                      "--layers", "kraken.genus.mg,prot"])
    assert rc == 0
    t = read_table(str(out), "cohort_tests")
    assert {"prot.genus.extra", "prot.genus.intra", "prot.kegg_module.extra", "prot.kegg_ko.intra"} <= set(t["layer"])
    a = t[(t["layer"] == "prot.genus.intra") & (t["feature"] == "Akkermansia") & t["comparison"].str.startswith("case_a vs")]
    assert len(a) and a["effect"].iloc[0] > 0.8
    st = read_manifest(str(out))["cohort"]["layers"]
    assert st["prot.genus.extra"]["samples"] == 29 and st["prot.genus.extra"]["label"] == "Metaproteome extracellular, genus"
    try:
        import dash  # noqa: F401
    except ImportError:
        return
    from impexplorer.app import create_app
    from impexplorer.app.pages import cohort as cpage, overview, proteomics as ppage, search, taxonomy
    app = create_app(str(out))
    assert "Proteomics" in [p.name for p in app._imp_pages]
    d = app._imp_data
    for p in app._imp_pages:
        assert p.layout(d) is not None
    # next to mg / mt in the shared views
    assert {"mg", "mt", "mp_extra", "mp_intra"} <= set(d.distinct("functions", "omic"))
    assert "metaproteomics" in set(d.distinct("taxonomy", "source"))
    ao = taxonomy.across_omics(d, "genus", "Akkermansia", d.samples)
    assert len(ao["proteins mp_extra"]) == 29 and len(ao["Kraken2 mg"]) == 30
    v, unit = cpage.feature_values(d, "prot.genus.extra", "Akkermansia")
    assert len(v) == 29 and "protein abundance" in unit
    hits = search.find(d, "akkermans")
    assert "taxonomy" in set(hits["layer"]) and any("mp_" in k for k in hits["kind"])
    M = ppage.term_genera(d, "mp_intra", "kegg_module", "M00000", list(d.samples["sample_id"]))
    assert abs(M.sum(axis=1) - 1).max() < 1e-9 and set(M.columns) <= {"Bacteroides", "Blautia_A", "Akkermansia",
                                                                       "Faecalibacterium"}
    G = ppage.genus_terms(d, "mp_intra", "kegg_ko", "Bacteroides", list(d.samples["sample_id"]))
    assert len(G) and (G.max() <= 1 + 1e-9).all().all()
    assert overview.prot_tiles(d) is not None and "proteins mp_extra" in overview.sample_table(d).columns


def test_pfam_cog_functions(tmp_path):
    """--function-dbs with pfam and cog: Pfam accessions only (no names), COG ids, links, cohort default layers"""
    from impexplorer import cohort
    from impexplorer.app.links import url
    root = make_fixture.main(str(tmp_path), n=12)
    out = tmp_path / "release"
    export.main(["--imp-root", root, "--metadata", str(tmp_path / "metadata.tsv"), "--out", str(out), "--threads", "1",
                 "--function-dbs", "kegg_ko,kegg_module,pfam,cog"])
    f = read_table(str(out), "functions")
    assert set(f[f["db"] == "pfam"]["feature"]) == {"PF00001", "PF00002", "PF00004", "PF00008", "PF00016"}
    assert set(f[f["db"] == "cog"]["feature"]) == {"COG0001", "COG0002", "COG0005", "COG0012"}
    assert url("pfam", "PF00001").startswith("https://www.ebi.ac.uk/interpro/") and "COG0005" in url("cog", "COG0005")
    names = tmp_path / "cog.def.tab"
    names.write_text("COG0001\tH\tGlutamate-1-semialdehyde aminotransferase\themL\n")
    assert cohort.main(["--data", str(out), "--permutations", "19", "--no-report", "--term-names", str(names),
                        "--model", "group"]) == 0
    t = read_table(str(out), "cohort_tests")
    assert {"pfam.mg", "pfam.mt", "cog.mg", "cog.mt"} <= set(t["layer"])
    assert t[t["feature"] == "COG0001"]["name"].iloc[0] == "Glutamate-1-semialdehyde aminotransferase"


def test_nmds_names_from_srna_config(tmp_path):
    """COHORT_NMDS of an sRNA-IMP config names sRNA sets: translated to the metagenome layers"""
    from impexplorer.cohort import settings
    cfg = tmp_path / "srna.config.sh"
    cfg.write_text('COHORT_NMDS="${COHORT_NMDS-genus,fn_kegg_pathway,mirna,rfam_family}"\n'
                   'COHORT_NMDS_STARTS="${COHORT_NMDS_STARTS-10}"\n')
    a = settings.parse(["--data", str(tmp_path), "--config", str(cfg)])
    assert a.nmds == "kraken.genus.mg,kegg_pathway.mg" and a.nmds_starts == 10
