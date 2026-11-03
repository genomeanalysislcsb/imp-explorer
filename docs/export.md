# Export IMP3 results

`imp-explorer export` reads the IMP3 output folder of every sample and writes one release directory. Run it where
the IMP3 results are (e.g. on the HPC); it needs only pandas.

## Samples from a root folder

Every sub-folder with IMP3 output (`Stats/` or `status/`) is one sample:

```bash
imp-explorer export --imp-root /path/to/IMP3_outputs --metadata /path/to/metadata.tsv \
    --study "My study" --out /path/to/release_v0.1 --threads 8
```

By default the folder name is the sample id.

## Samples from a sheet

```bash
imp-explorer export --samples /path/to/samples.tsv --out /path/to/release_v0.1
```

`samples.tsv` (or `.csv`) has the columns `sample_id`, `imp_dir` and any metadata columns. Relative `imp_dir`
paths are relative to the sheet's folder.

```text
sample_id   imp_dir                       group     site
S01         /path/to/IMP3_outputs/s01     control   site_1
S02         /path/to/IMP3_outputs/s02     case      site_1
```

## Sample ids that differ from the folder names

The ids of the IMP3 folders have to match the ids of the metadata. When they don't (e.g. folders `run_s01`,
`run_s02` and metadata ids `S01`, `S02`), translate the folder names with one or more of these options, applied in
this order:

| option | example | `run_s01` becomes |
|---|---|---|
| `--id-map FILE` | TSV/CSV with two columns: folder name, sample id (an optional header line is skipped) | the id in the file; folders not in the file keep their name |
| `--id-regex REGEX` | `--id-regex 'run_(\w+)'` | the first capture group (or the whole match): `s01` |
| `--id-case upper\|lower` | `--id-case upper` | `S01` (here after the regex) |

```bash
imp-explorer export --imp-root /path/to/IMP3_outputs --id-regex 'run_(\w+)' --id-case upper \
    --metadata /path/to/metadata.tsv --out /path/to/release_v0.1
```

The export stops if a folder does not match the regex, or if two folders end up with the same id. It reports how
many samples found their metadata and lists metadata ids without IMP3 output (also in `manifest.json` under
`metadata_without_output`), so a wrong translation shows at once. Alternatively, write a [sample sheet](#samples-from-a-sheet).

## Metadata

`--metadata` is a TSV or CSV with one row per sample; the id column is `sample_id` or, if missing, the first column.
All columns end up in the `samples` table. In the app, every column with 2–12 distinct values can be used to group
samples, and every column to filter them.

!!! warning "Check the design before comparing groups"
    If a group was sampled at one site only, group and site differences cannot be separated. Use the sample filter
    of the app to compare groups within one site.

## Paths in the release

The `samples` table keeps the IMP3 folder of each sample (`imp_dir`). By default it is stored relative to
`--imp-root` (or the sample sheet), so a release does not contain the folder structure of your HPC. Use
`--paths absolute` to keep full paths, or `--paths none` to leave them out, e.g. for a public release.

## Options

| option | default | meaning |
|---|---|---|
| `--imp-root` / `--samples` | – | where the samples come from (one of the two) |
| `--metadata` | – | sample metadata joined by `sample_id` |
| `--id-map`, `--id-regex`, `--id-case` | – / – / `keep` | translate folder names into sample ids (with `--imp-root`) |
| `--paths` | `relative` | `relative`, `absolute` or `none`: how `imp_dir` is stored |
| `--out` | – | release directory to write |
| `--study` | `study` | study name shown in the app |
| `--function-dbs` | `kegg_ko,kegg_module,kegg_pathway,cazy,enzyme_ec` | Mantis databases to export: `<omic>.mantis.<db>.counts.tsv` of each sample. `pfam` and `cog` are supported throughout (Pfam accessions / COG ids only, links to InterPro / NCBI, cohort layers `pfam.mg` … by default when exported); `go`, `eggnog`, … also work but are large |
| `--format` | `auto` | `parquet` (needs pyarrow), `tsv` (TSV.gz) or `auto` |
| `--threads` | 4 | samples read in parallel |

## Output

```text
/path/to/release_v0.1/
├── manifest.json          study, date, versions, rows per table, failed samples, unmatched metadata ids
├── samples.parquet
├── qc.parquet
├── taxonomy.parquet
├── bins.parquet
├── bin_abundance.parquet
└── functions.parquet
```

A sample that cannot be read is listed under `failed_samples` in the manifest (exit code 2); the others are
exported. Missing files of a partly finished sample give empty rows, not an error. For about 200 integrated
metagenome + metatranscriptome samples, the release is ~30 MB and the export takes a few minutes.

The tables are described in [Release tables](data-model.md), the IMP3 files behind them in
[IMP3 files read](imp3-files.md).
