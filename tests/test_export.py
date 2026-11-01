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
