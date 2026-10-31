# Development

```text
impexplorer/
├── model.py            release tables (columns, keys, docs), MIMAG tiers
├── readers/imp3.py     IMP3 sample folder -> tables
├── export.py           imp-explorer export
├── release.py          write / read tables and the manifest
├── cli.py              imp-explorer {export,run}
└── app/
    ├── __init__.py     create_app(), page registry, routing, global controls
    ├── data.py         DuckDB views over the release, samples and metadata
    ├── stats.py        Bray-Curtis, PCoA, PERMANOVA, group tests with BH-FDR
    ├── links.py        links to KEGG, CAZy, GTDB
    ├── ui.py           Page, palette, figure style, small components
    ├── assets/style.css
    └── pages/          overview, qc, taxonomy, mags, functions, search (+ _common helpers)
tests/
├── make_fixture.py     synthetic IMP3 sample folders in the real file formats
├── test_export.py
└── run_tests.py        runs the tests without pytest
```

## Tests

```bash
python tests/run_tests.py
python tests/make_fixture.py /tmp/demo 30          # a demo study of 30 samples
imp-explorer export --imp-root /tmp/demo/imp --metadata /tmp/demo/metadata.tsv --out /tmp/demo/release
imp-explorer run --data /tmp/demo/release --debug
```

The fixture writes the IMP3 file formats exactly as IMP3 3.0 does (both Bracken layouts, plain and hierarchical
mOTUs, space-padded bin tables), so reader changes are tested against them. The documentation screenshots are
taken from this demo study.

## Adding an IMP3 output

1. Add the columns to `TABLES` in `model.py` (or a new table).
2. Read the file in `readers/imp3.py`; return an empty frame when the file is missing.
3. Write the format into `tests/make_fixture.py` and test it in `tests/test_export.py`.
4. Show it on a page.
