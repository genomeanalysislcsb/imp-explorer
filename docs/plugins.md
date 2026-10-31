# Plug-ins and study portals

The IMP explorer is the core of study portals. A portal adds its own data as extra release tables and its own pages
as a plug-in; the core stays unchanged and keeps working for studies with IMP3 data only.

## A page

A page is an `impexplorer.app.ui.Page`:

| field | meaning |
|---|---|
| `path`, `name`, `order` | URL, menu label, position in the menu (core: 0, 10, … 90) |
| `layout(data)` | returns the Dash component tree of the page |
| `register(app, data)` | registers the page's callbacks |
| `requires` | tables that must exist and be non-empty; otherwise the page is hidden |

`data` (`impexplorer.app.data.Data`) gives SQL access to all tables (`data.q(sql, params)`, DuckDB), the samples and
metadata (`data.samples`, `data.group_cols`, `data.meta_cols`) and helpers for the global controls
(`impexplorer.app.pages._common.selection`).

```python
# my_portal/pages.py
from dash import Input, Output, html
from impexplorer.app.ui import GLOBAL, Page, card, graph
from impexplorer.app.pages._common import box_by_group, selection


def layout(data):
    return html.Div([html.H2("Small RNA"), card("miRNA-like reads by group", graph("srna-box"))])


def register(app, data):
    @app.callback(Output("srna-box", "figure"),
                  Input(GLOBAL[0], "value"), Input(GLOBAL[1], "value"), Input(GLOBAL[2], "value"))
    def _box(group, fcol, fval):
        s, g, cmap = selection(data, group, fcol, fval)
        d = data.q("SELECT sample_id, value FROM srna_qc WHERE metric = 'miRNA_like_fraction'")
        return box_by_group(d.set_index("sample_id")["value"], g, cmap, "share of non-human reads")


PAGES = [Page(path="/srna", name="sRNA", order=50, layout=layout, register=register, requires=("srna_qc",))]
```

```bash
imp-explorer run --data release_v0.1 --plugin my_portal.pages --table srna_qc
```

An installed package can also register pages through the entry-point group `impexplorer.pages`:

```toml
[project.entry-points."impexplorer.pages"]
my_portal = "my_portal.pages:PAGES"
```

## Extra tables

Write extra tables into the release directory as `<name>.parquet` (or `.tsv.gz`) with the shared keys
(`sample_id`, `omic`, a feature column). Pass them with `--table <name>`; they are then available as DuckDB views of
the same name.

## Example: a multi-omics study portal

A study with more data layers than IMP3 covers uses the IMP explorer for metagenomics and metatranscriptomics and
adds, as plug-ins: small-RNA results (e.g. from [sRNA-IMP](https://github.com/genomeanalysislcsb/sRNA-IMP)),
metaproteomics, metabolomics or host data as extra tables with the same keys, a cross-layer feature search, and
pages for the figures of its paper.
