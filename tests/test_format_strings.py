"""every string literal formatted with % in the package has as many placeholders as arguments and no '% ' (a
literal percent sign followed by a space, e.g. '% of ...', is read as the conversion '% o' and fails at run time
only, in an app callback)"""
import ast
import os
import re

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "impexplorer")
SPEC = re.compile(r"%(\([^)]*\))?([-+ #0]*)(\d+|\*)?(?:\.(\d+|\*))?([diouxXeEfFgGcrsa%])")


def placeholders(s):
    """number of arguments the format string takes, and the specs with a space flag"""
    n, spaced = 0, []
    for m in SPEC.finditer(s):
        key, flags, w, prec, t = m.groups()
        if t == "%":
            continue
        if " " in flags:
            spaced.append(m.group(0))
        n += (w == "*") + (prec == "*") + 1
    return n, spaced


def test_percent_format_strings():
    bad = []
    for d, _, files in os.walk(ROOT):
        for f in files:
            if not f.endswith(".py"):
                continue
            p = os.path.join(d, f)
            for n in ast.walk(ast.parse(open(p, encoding="utf-8").read())):
                if not (isinstance(n, ast.BinOp) and isinstance(n.op, ast.Mod) and isinstance(n.left, ast.Constant)
                        and isinstance(n.left.value, str)):
                    continue
                s = n.left.value
                where = "%s:%d %r" % (os.path.relpath(p, ROOT), n.lineno, s[:60])
                if "%(" in s:
                    continue                                    # mapping keys: checked where used
                k, spaced = placeholders(s)
                if spaced:
                    bad.append("%s: %s (a literal '%%' must be written '%%%%')" % (where, ", ".join(spaced)))
                if isinstance(n.right, ast.Tuple) and len(n.right.elts) != k:
                    bad.append("%s: %d placeholders, %d arguments" % (where, k, len(n.right.elts)))
                elif not isinstance(n.right, ast.Tuple) and k != 1 and not isinstance(n.right, (ast.Name,
                                                                                              ast.Call)):
                    bad.append("%s: %d placeholders, 1 argument" % (where, k))
    assert not bad, "\n".join(bad)
