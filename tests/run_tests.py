"""Run the tests without pytest: python3 tests/run_tests.py"""
import inspect, pathlib, sys, tempfile, traceback, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import test_export as T
fails = 0
for name, fn in inspect.getmembers(T, inspect.isfunction):
    if not name.startswith("test_"):
        continue
    with tempfile.TemporaryDirectory() as d:
        try:
            fn(pathlib.Path(d)); print("ok   ", name)
        except Exception:
            fails += 1; print("FAIL ", name); traceback.print_exc()
sys.exit(1 if fails else 0)
