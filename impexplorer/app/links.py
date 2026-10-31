"""Links from feature ids to their reference pages (KEGG, CAZy, ExPASy, GTDB)."""
import re
from urllib.parse import quote

from dash import html


def url(db, feature):
    f = str(feature).strip()
    if db == "kegg_ko" and re.fullmatch(r"K\d{5}", f):
        return "https://www.genome.jp/dbget-bin/www_bget?ko:" + f
    if db == "kegg_pathway" and re.fullmatch(r"(ko|map)\d{5}", f):
        return "https://www.genome.jp/dbget-bin/www_bget?pathway+" + f
    if db == "kegg_module" and re.fullmatch(r"M\d{5}", f):
        return "https://www.kegg.jp/module/" + f
    if db == "enzyme_ec" and re.fullmatch(r"\d+\.\d+\.\d+\.[\dn-]+", f):
        return "https://www.genome.jp/dbget-bin/www_bget?ec:" + f
    if db == "cazy" and re.fullmatch(r"(GH|GT|PL|CE|AA|CBM)\d+", f):
        return "http://www.cazy.org/%s.html" % f
    if db == "gtdb" and f:
        return "https://gtdb.ecogenomic.org/searches?s=al&q=" + quote(f)
    return None


SITE = {"kegg_ko": "KEGG", "kegg_pathway": "KEGG", "kegg_module": "KEGG", "enzyme_ec": "KEGG", "cazy": "CAZy",
        "gtdb": "GTDB"}


def md(db, feature):
    """markdown link for an AG Grid column with cellRenderer 'markdown' ('' if none)"""
    u = url(db, feature)
    return "[%s ↗](%s)" % (SITE.get(db, "link"), u) if u else ""


def anchor(db, feature):
    u = url(db, feature)
    if not u:
        return None
    return html.A("%s in %s ↗" % (feature, SITE.get(db, "reference")), href=u, target="_blank",
                  className="extlink")
