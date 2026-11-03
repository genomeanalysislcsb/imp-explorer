"""imp-explorer command line: export | cohort | srna | metabolomics | run"""
import sys


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in ("-h", "--help"):
        print("usage: imp-explorer {export,cohort,srna,metabolomics,run} ...\n  export  IMP3 outputs -> release tables\n"
              "  cohort  cohort analysis of a release (tables, PDF/HTML report, app page)\n"
              "  srna    add the results of an sRNA-IMP cohort analysis to a release (optional layer)\n"
              "  metabolomics  add metabolomics results (generic assay / feature / value tables) to a release\n"
              "  run     start the explorer on a release directory")
        return 0
    cmd, rest = argv[0], argv[1:]
    if cmd == "export":
        from .export import main as m
    elif cmd == "srna":
        from .srna import main as m
    elif cmd == "metabolomics":
        from .metabolomics import main as m
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
