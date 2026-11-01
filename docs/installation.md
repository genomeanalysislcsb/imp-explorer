# Installation

The export needs only Python ≥ 3.9 with pandas (pyarrow for Parquet output). The app needs Dash, Plotly, DuckDB,
PyArrow and SciPy; the cohort analysis numpy and pandas, plus matplotlib for its PDF/HTML report. One conda
environment covers all three.

## Conda / micromamba (recommended)

```bash
git clone https://github.com/genomeanalysislcsb/imp-explorer.git
cd imp-explorer
micromamba create -y -p /path/to/envs/imp_explorer -f environment.yml
micromamba activate /path/to/envs/imp_explorer
pip install --no-deps -e .
python tests/run_tests.py            # all tests ok
```

`environment.yml` pins Dash to 2.x (`dash>=2.17,<3`) and Plotly to 5.x (`plotly>=5.22,<6`); all packages come from
conda-forge.

!!! warning "Plotly 6 or 7 with Dash 2"
    Plotly ≥ 6 sends numbers to the browser as binary typed arrays, which the plotly.js shipped with Dash 2 does not
    fully read: figures can stay empty. The app avoids the known case (integer arrays), but an existing environment
    with a newer Plotly is safest pinned back: `micromamba install -n imp_explorer -c conda-forge "plotly>=5.22,<6"`.

### Updating

```bash
cd imp-explorer && git pull origin main       # with pip install -e . nothing else is needed
imp-explorer --help                           # lists export, cohort, run
```

## pip

```bash
pip install -e .              # export only (pandas, numpy)
pip install -e ".[parquet]"   # + pyarrow: Parquet output
pip install -e ".[app]"       # + the app (dash, dash-ag-grid, plotly, duckdb, pyarrow, scipy)
pip install -e ".[cohort]"    # + the cohort report (matplotlib, scipy)
```

## Offline machines

For a machine without internet access, download the wheels elsewhere and install from the folder:

```bash
pip download dash==2.18.2 dash-ag-grid==31.2.0 plotly==5.24.1 duckdb==1.1.3 pyarrow==18.1.0 scipy matplotlib \
    --dest wheels --only-binary=:all: --platform manylinux2014_x86_64 --python-version 3.12
pip install --no-index --find-links wheels dash dash-ag-grid plotly duckdb pyarrow scipy matplotlib
```
