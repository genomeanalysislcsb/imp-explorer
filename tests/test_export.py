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
    b = os.path.join(root, "BHC01")
    q = imp3.qc(b, "BHC01")
    m = dict(zip(q["omic"] + ":" + q["metric"], q["value"]))
    assert m["mg:trimmed:pairs"] > 0 and 0.9 < m["mg:mapped_fraction"] <= 1
    assert m["mgmt:contigs"] == 400 and m["mgmt:n50"] > 0
    t = imp3.taxonomy(b, "BHC01")
    assert set(t["source"]) == {"kraken2_bracken", "motus"}
    g = t[(t["omic"] == "mg") & (t["rank"] == "genus")]
    assert abs(g["rel_abundance"].sum() - 1) < 1e-9 and "unclassified" not in set(g["taxon"])
    bins = imp3.bins(b, "BHC01")
    assert len(bins) >= 6 and bins["genus"].str.len().gt(0).all() and set(bins["quality"]) <= {"high", "medium", "low"}
    assert bins["bin"].str.startswith("bin.").all()          # padded cells stripped
    ab = imp3.bin_abundance(b, "BHC01")
    assert set(ab["omic"]) == {"mg", "mt"} and ab["share_mapped"].between(0, 1).all()
    f = imp3.functions(b, "BHC01")
    assert {"kegg_ko", "kegg_module"} == set(f["db"]) and (f["reads"] > 0).all()


def test_new_bracken_layout(tmp_path):
    root, _ = build(tmp_path)
    sid = "KPC02"                                              # odd samples use Analysis/taxonomy/bracken/<omic>/
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
