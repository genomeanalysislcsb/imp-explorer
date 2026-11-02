"""Run the tests without pytest: python3 tests/run_tests.py (pytest works too)"""
import inspect
import os
import pathlib
import sys
import tempfile
import traceback

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import test_export as T  # noqa: E402


def main():
    fails = 0
    for name, fn in inspect.getmembers(T, inspect.isfunction):
        if not name.startswith("test_"):
            continue
        with tempfile.TemporaryDirectory() as d:
            try:
                fn(pathlib.Path(d)) if inspect.signature(fn).parameters else fn()
                print("ok   ", name)
            except Exception:
                fails += 1
                print("FAIL ", name)
                traceback.print_exc()
    return 1 if fails else 0


if __name__ == "__main__":   # guard: export workers are spawned and import this module
    sys.exit(main())
