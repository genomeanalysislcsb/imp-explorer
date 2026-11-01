"""imp-explorer command line: export | run | cohort"""
import sys


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in ("-h", "--help"):
        print("usage: imp-explorer {export,run,cohort} ...\n  export  IMP3 outputs -> release tables\n"
              "  cohort  cohort analysis of a release (tables, PDF/HTML report, app page)\n"
              "  run     start the explorer on a release directory")
        return 0
    cmd, rest = argv[0], argv[1:]
    if cmd == "export":
        from .export import main as m
    elif cmd == "cohort":
        from .cohort import main as m
    elif cmd == "run":
        from .app import main as m
    else:
        print("unknown command %s" % cmd, file=sys.stderr)
        return 1
    return m(rest)


if __name__ == "__main__":
    sys.exit(main())
