# Installation

The export needs only Python ≥ 3.9 with pandas (pyarrow for Parquet output). The app needs Dash, Plotly, DuckDB,
PyArrow and SciPy. One conda environment covers both.

## Conda / micromamba (recommended)

```bash
git clone https://github.com/genomeanalysislcsb/imp-explorer.git
cd imp-explorer
micromamba create -y -p /path/to/envs/imp_explorer -f environment.yml
micromamba activate /path/to/envs/imp_explorer
pip install --no-deps -e .
python tests/run_tests.py            # 5 x ok
```

`environment.yml` pins Dash to 2.x (`dash>=2.17,<3`); all packages come from conda-forge.

## pip

```bash
pip install -e .              # export only (pandas, numpy)
pip install -e ".[parquet]"   # + pyarrow: Parquet output
pip install -e ".[app]"       # + the app (dash, dash-ag-grid, plotly, duckdb, pyarrow, scipy)
```

## Offline machines

For a machine without internet access, download the wheels elsewhere and install from the folder:

```bash
pip download dash==2.18.2 dash-ag-grid==31.2.0 plotly==5.24.1 duckdb==1.1.3 pyarrow==18.1.0 scipy \
    --dest wheels --only-binary=:all: --platform manylinux2014_x86_64 --python-version 3.12
pip install --no-index --find-links wheels dash dash-ag-grid plotly duckdb pyarrow scipy
```
