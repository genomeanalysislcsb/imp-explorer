"""Shared UI pieces: page registry, palette, figure style, small components."""
from dataclasses import dataclass, field
from typing import Callable, Tuple

import plotly.graph_objects as go
from dash import dcc, html

PALETTE = ["#2a6f97", "#d1495b", "#edae49", "#00798c", "#8f5aa8", "#66a182", "#e07a5f", "#3d405b",
           "#81b29a", "#f2cc8f", "#6d597a", "#b56576"]
OMIC_LABEL = {"mg": "metagenome (mg)", "mt": "metatranscriptome (mt)", "mgmt": "assembly"}
GLOBAL = ("g-group", "g-fcol", "g-fval")       # ids of the global controls (group by, filter column, values)


@dataclass
class Page:
    path: str
    name: str
    order: int
    layout: Callable                    # layout(data) -> component
    register: Callable                  # register(app, data) -> None (callbacks)
    requires: Tuple[str, ...] = field(default_factory=tuple)   # tables that must be non-empty

    def available(self, data):
        return all(data.has(t) for t in self.requires)


def colors(levels):
    return {str(g): PALETTE[i % len(PALETTE)] for i, g in enumerate(sorted(map(str, levels)))}


def style(fig, height=420, legend=True):
    fig.update_layout(template="simple_white", height=height, margin=dict(l=60, r=20, t=40, b=50),
                      font=dict(family="Inter, Helvetica, Arial, sans-serif", size=12, color="#1f2933"),
                      legend=dict(orientation="h", y=-0.18, x=0) if legend else dict(),
                      showlegend=legend, title=dict(font=dict(size=14), x=0, xanchor="left"),
                      paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)")
    return fig


def empty_fig(msg="no data for this selection", height=300):
    fig = go.Figure()
    fig.add_annotation(text=msg, showarrow=False, font=dict(size=13, color="#7b8794"), x=0.5, y=0.5,
                       xref="paper", yref="paper")
    fig.update_xaxes(visible=False)
    fig.update_yaxes(visible=False)
    return style(fig, height, legend=False)


def card(title, *children, note=None, cls=""):
    return html.Div([html.Div([html.H3(title), html.P(note, className="note") if note else None],
                              className="card-head"), *children], className="card " + cls)


def control(label, comp):
    return html.Div([html.Label(label), comp], className="control")


def dd(id_, options, value=None, multi=False, clearable=False, placeholder=None):
    opts = [o if isinstance(o, dict) else {"label": str(o), "value": o} for o in options]
    return dcc.Dropdown(id=id_, options=opts, value=value, multi=multi, clearable=clearable,
                        placeholder=placeholder)


def stat_tile(value, label):
    return html.Div([html.Div(value, className="tile-value"), html.Div(label, className="tile-label")],
                    className="tile")


def fmt_int(x):
    try:
        return "{:,}".format(int(round(float(x))))
    except (TypeError, ValueError):
        return "–"


def graph(id_, height=None):
    return dcc.Loading(dcc.Graph(id=id_, config={"displaylogo": False, "toImageButtonOptions": {"format": "svg"}},
                                 style={"height": height} if height else None), type="dot", color="#2a6f97")
