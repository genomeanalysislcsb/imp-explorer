"""imp-explorer run: the Dash app on a release directory.

Pages are plug-ins (impexplorer.app.ui.Page). The core pages cover IMP3 (overview, QC, taxonomy, MAGs, functions,
search); a study portal adds its own pages with --plugin <module> (a module with a list PAGES) or through the
entry-point group 'impexplorer.pages'. A page only appears when the tables it needs are in the release.
"""
import argparse
import importlib
import sys

import dash
from dash import Input, Output, dcc, html

from .data import Data
from .ui import GLOBAL, Page, control, dd


def nav_id(path):
    return "nav-" + (path.strip("/").replace("/", "-") or "home")


def core_pages():
    from .pages import cohort, functions, mags, metabolomics, overview, proteomics, qc, search, srna, taxonomy
    return [overview.PAGE, qc.PAGE, taxonomy.PAGE, mags.PAGE, functions.PAGE, srna.PAGE, proteomics.PAGE,
            metabolomics.PAGE, search.PAGE, cohort.PAGE]


def plugin_pages(modules=()):
    pages = []
    for m in modules:
        pages += list(importlib.import_module(m).PAGES)
    try:
        from importlib.metadata import entry_points
        for ep in entry_points(group="impexplorer.pages"):
            pages += list(ep.load())
    except Exception:                                       # no entry points: fine
        pass
    return pages


def create_app(release, plugins=(), extra_tables=()):
    data = Data(release, extra_tables)
    pages = [p for p in core_pages() + plugin_pages(plugins) if isinstance(p, Page)]
    pages = sorted([p for p in pages if p.available(data)], key=lambda p: p.order)
    app = dash.Dash(__name__, suppress_callback_exceptions=True, title="%s · IMP explorer" % data.study,
                    update_title=None)
    gcols = data.group_cols
    default_group = next((c for c in ("group", "Group", "condition", "diagnosis") if c in gcols), gcols[0] if gcols else None)
    app.layout = html.Div([
        dcc.Location(id="url"),
        html.Div([
            html.Div([data.study, html.Span("IMP explorer")], className="brand"),
            html.Nav([dcc.Link(p.name, href=p.path, id=nav_id(p.path)) for p in pages], className="nav"),
            html.Div([
                control("group by", dd(GLOBAL[0], [{"label": "(none)", "value": ""}] + gcols, default_group or "")),
                control("filter", dd(GLOBAL[1], [{"label": "(all samples)", "value": ""}] + data.meta_cols, "")),
                control("values", dd(GLOBAL[2], [], [], multi=True, placeholder="all")),
            ], className="globals"),
        ], className="topbar"),
        html.Div(id="page", className="main"),
    ])

    @app.callback(Output(GLOBAL[2], "options"), Output(GLOBAL[2], "value"), Input(GLOBAL[1], "value"))
    def _filter_values(col):
        if not col:
            return [], []
        vals = sorted(data.samples[col].dropna().astype(str).unique())
        return [{"label": v, "value": v} for v in vals], []

    by_path = {p.path: p for p in pages}

    @app.callback(Output("page", "children"), *[Output(nav_id(p.path), "className") for p in pages],
                  Input("url", "pathname"))
    def _route(path):
        p = by_path.get(path or "/") or pages[0]
        return [p.layout(data)] + ["active" if q is p else "" for q in pages]

    for p in pages:
        p.register(app, data)
    app._imp_data, app._imp_pages = data, pages
    return app


def main(argv=None):
    ap = argparse.ArgumentParser(prog="imp-explorer run", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True, help="release directory written by 'imp-explorer export'")
    ap.add_argument("--plugin", action="append", default=[], help="module with extra PAGES (repeatable)")
    ap.add_argument("--table", action="append", default=[], help="extra release table a plugin needs (repeatable)")
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8050)
    ap.add_argument("--debug", action="store_true")
    a = ap.parse_args(argv)
    app = create_app(a.data, a.plugin, a.table)
    print("IMP explorer: http://%s:%d  (%s)" % (a.host, a.port, a.data), file=sys.stderr)
    app.run(host=a.host, port=a.port, debug=a.debug)
    return 0
