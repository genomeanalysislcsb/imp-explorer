# Metabolomics layer

`imp-explorer metabolomics` adds metabolomics results of the same samples to a release. It reads three generic
tables, so any assay works: targeted panels, untargeted features, lipids, or any other per-sample measurement.

```bash
imp-explorer metabolomics --data /path/to/release_v0.1 --metab /path/to/metabolomics_export
# then, for group comparisons, models and PERMANOVA of the assays next to the IMP3 layers:
imp-explorer cohort --data /path/to/release_v0.1 --config study.cohort.yaml
```

## Input: three tables

Each table can be `.parquet`, `.tsv` or `.tsv.gz` in the `--metab` folder.

| table | required columns | further columns (kept, shown in the app) |
|---|---|---|
| `metabolomics_assays` | `assay` | `label`, `platform`, `matrix`, `acquisition_date`, `value_type`, `units`, `notes` |
| `metabolomics_features` | `assay`, `feature_id` | `name`, `annotation`, `class`, `formula`, `mz`, `rt_min`, `ri`, `internal_standard`, `lab_notes`, ... |
| `metabolomics_values` | `sample_id`, `assay`, `feature_id`, `value` | |

`metabolomics_samples` (`sample_id`, `assay`, `source_column`) is optional.

- An **assay** is one platform or panel, e.g. `scfa`, `bile_acids` or `lipids_targeted`.
- A **feature** is a measured compound or an unknown (m/z / RT / RI).
- A **value** is a concentration, a normalised intensity or a ratio. Values should be ≥ 0, with 0 meaning not
  detected. Samples that were not measured in an assay simply have no rows for it.
- `annotation`: values starting with `unknown` or `name from mass` count as not identified ("identified only" on the
  Metabolomics page).

For Expobiome, `metabolomics_export.py` writes these tables from the final sheets of the lab's Excel workbooks.

## Options

| option | default | meaning |
|---|---|---|
| `--data` | – | release directory |
| `--metab` | – | folder with the tables |
| `--assays` | all | only these assays (comma-separated) |
| `--id-map`, `--id-regex`, `--id-case` | – | map the sample ids of the tables to the release ids, as for [export](export.md) |
| `--format` | auto | Parquet (with pyarrow) or TSV.gz |

The log reports how many samples are in the release, and lists the ones that aren't.

## In the release

- **Tables:** `metab_assays` (with `n_samples`, `n_features`, `n_in_release`), `metab_features`, `metab_values`
  (`sample_id`, `assay`, `feature_id`, `value`, `source_id`), `metab_samples`.
- **Manifest:** a `metabolomics` section with the assays, units and counts.

Running the command again replaces the tables.

## In the cohort analysis

`imp-explorer cohort` analyses every assay as a layer `metab.<assay>`. They're included by default when the release
has metabolomics; with `--layers`, add `metab` for all assays or `metab.<assay>` for single ones. The statistics are
the same as for the other layers, adapted to continuous values:

| step | metabolomics layers |
|---|---|
| features | detected (> 0) in ≥ `--prevalence` of the samples and not constant; features missing (no value) in > 20 % of the samples are left out, then samples with a missing value |
| group comparisons | van Elteren / Mann-Whitney / Kruskal-Wallis on the values (rank-based, so scale-free) |
| adjusted models | OLS with HC3 on log2(value + half the smallest positive value of the feature): estimates in log2 units, e.g. 1 = twofold |
| beta diversity, PERMANOVA, dispersion | Euclidean distances of the autoscaled (z-scored) log2 values |

The results appear on the app's Cohort page as "Metabolomics: …" layers, and in the report.

## In the app

- **Metabolomics page:**
  - profiles of each assay as a clustered heatmap (all features or identified only);
  - one feature by group, in the units of the assay;
  - the feature's Spearman correlation with the taxa or functions of the IMP3 profiles of the same samples
    (exploratory);
  - a table of all features with annotation, detection and median per group.
- **Search:** finds metabolites by name or feature id and plots the values by group.
- **Overview:** a row of key numbers, and one column per assay in the sample table.
