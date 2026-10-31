# Export IMP3 results

`imp-explorer export` reads the IMP3 output folder of every sample and writes one release directory. Run it where
the IMP3 results are (e.g. on the HPC); it needs only pandas.

## Samples from a root folder

Every sub-folder with IMP3 output (`Stats/` or `status/`) is one sample, named after the folder:

```bash
imp-explorer export --imp-root /work/projects/mystudy/IMP --id-case upper \
    --metadata metadata.tsv --study "My study" --out release_v0.1 --threads 8
```

`--id-case upper` turns folder names like `bhc01` into the ids used in the metadata (`BHC01`).

## Samples from a sheet

```bash
imp-explorer export --samples samples.tsv --out release_v0.1
```

`samples.tsv` (or `.csv`) has the columns `sample_id`, `imp_dir` and any metadata columns.

## Metadata

`--metadata` is a TSV or CSV with one row per sample; the id column is `sample_id` or, if missing, the first column.
All columns end up in the `samples` table. In the app, every column with 2–12 distinct values can be used to group
samples, and every column to filter them.

!!! warning "Check the design before comparing groups"
    If a group was sampled at one site only (e.g. all patients of one disease from one centre), group and site
    differences cannot be separated. Use the sample filter of the app to compare groups within one site.

## Options

| option | default | meaning |
|---|---|---|
| `--imp-root` / `--samples` | – | where the samples come from (one of the two) |
| `--metadata` | – | sample metadata joined by `sample_id` |
| `--id-case` | `keep` | `upper` / `lower`: change the case of folder-derived ids |
| `--out` | – | release directory to write |
| `--study` | `study` | study name shown in the app |
| `--function-dbs` | `kegg_ko,kegg_module,kegg_pathway,cazy,enzyme_ec` | Mantis databases to export (`go`, `pfam`, `cog`, `eggnog`, … also work but are large) |
| `--format` | `auto` | `parquet` (needs pyarrow), `tsv` (TSV.gz) or `auto` |
| `--threads` | 4 | samples read in parallel |

## Output

```text
release_v0.1/
├── manifest.json          study, date, versions, rows per table, failed samples
├── samples.parquet
├── qc.parquet
├── taxonomy.parquet
├── bins.parquet
├── bin_abundance.parquet
└── functions.parquet
```

A sample that cannot be read is listed under `failed_samples` in the manifest (exit code 2); the others are
exported. Missing files of a partly finished sample give empty rows, not an error. For 180 integrated
metagenome + metatranscriptome samples, the release is ~30 MB and the export takes a few minutes.

The tables are described in [Release tables](data-model.md), the IMP3 files behind them in
[IMP3 files read](imp3-files.md).
