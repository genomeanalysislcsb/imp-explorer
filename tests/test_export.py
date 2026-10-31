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
    assert set(mo["rank"]) == {"genus", "motu"} and "unassigned" not in set(mo["taxon"])
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
