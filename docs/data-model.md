# Release tables (data model)

A release is a directory of long-format tables that join on the same keys, plus `manifest.json`. Every layer of a
study (IMP3 here; sRNA, proteomics, metabolomics, immune data in a study portal) uses the same keys, so any layer can
be combined with any other and with the metadata.

| key | meaning | example |
|---|---|---|
| `sample_id` | one sample (one IMP3 run) | `S01` |
| `omic` | data layer of the sample | `mg` metagenome, `mt` metatranscriptome, `mgmt` assembly |
| feature | the measured thing, with its namespace | `taxon` + `rank`, `bin_id`, `feature` + `db` |

## samples

| column | type | content |
|---|---|---|
| `sample_id` | str | sample id |
| `imp_dir` | str | IMP3 output folder, relative to `--imp-root` / the sample sheet by default (`--paths`); not shown in the app |
| `has_mg`, `has_mt` | bool | metagenome / metatranscriptome present |
| *metadata* | any | every column of `--metadata` |

## qc

`sample_id`, `omic`, `metric`, `value`. Metrics:

| metric | omic | source |
|---|---|---|
| `raw:pairs`, `trimmed:pairs`, `*:avg_len`, `*:q30_pct`, `*:gc_pct` | mg, mt | seqkit stats of R1 (`Stats/<omic>/<omic>.read_stats.tsv`) |
| `nonpareil_diversity` | mg, mt | mean of R1 and R2 Nonpareil diversity |
| `flagstat:total`, `flagstat:primary`, `flagstat:mapped`, `flagstat:primary_mapped`, `mapped_fraction` | mg, mt | samtools flagstat of the reads on the assembly |
| `contigs`, `assembly_bp`, `contigs_ge_1kb`, `n50`, `longest_contig` | mgmt | contig lengths |
| `imp:<step>` | mg, mt, mgmt | every row of IMP3's `Stats/all_stats.tsv` (reads per filtering step, annotation counts, bins) |

## taxonomy

| column | content |
|---|---|
| `sample_id`, `omic` | |
| `source` | `kraken2_bracken` or `motus` |
| `rank` | domain, phylum, class, order, family, genus, species. mOTUs: each mOTU (a species-level cluster) is a `species` row with the mOTU id as `taxid`; mOTUs with the same species name are summed in the app. Releases exported before 0.1.1 have all mOTUs rows as rank `motu` with the full lineage as name; the cohort analysis still reads them, a new export gives the ranks. |
| `taxon`, `taxid` | name and id (of the Kraken2 database set in IMP3, NCBI or GTDB; mOTU id for mOTUs) |
| `reads` | Bracken estimated reads / mOTUs read counts |
| `rel_abundance` | reads / all reads of the profile (sample, omic, source, rank); unclassified reads excluded |

## bins

One row per MAG: `bin_id` (`<sample>:<bin>`), `sample_id`, `bin`, `completeness`, `contamination` (CheckM2),
`quality` (MIMAG: high ≥ 90 / ≤ 5, medium ≥ 50 / ≤ 10, else low), `contig_n50`, `genome_size`, `gc`, the GTDB
ranks `domain` … `species`, `classification` (full GTDB string), `closest_genome_reference`, `closest_genome_ani`,
`grid` and `coverage` (GRiD growth rate, when IMP3 ran it).

## bin_abundance

| column | content |
|---|---|
| `sample_id`, `bin_id`, `omic` | |
| `reads` | reads of the omic on the bin's contigs |
| `length_bp` | summed contig length |
| `mean_depth` | length-weighted mean contig depth |
| `share_mapped` | reads / all reads of the omic on the assembly |
| `rpkm` | reads per kb of bin per million mapped reads |

## functions

`sample_id`, `omic`, `db`, `feature`, `reads` (featureCounts with `-M --fraction`, paired + single reads). The app
normalises within (sample, omic, db).

## Cohort tables

Written into the release by [`imp-explorer cohort`](cohort.md) (optional; the app's Cohort page needs them):
`cohort_samples` (metrics, group, block per sample), `cohort_qc`, `cohort_covariates`, `cohort_metric_tests`,
`cohort_tests` (layer, feature, comparison, kind, n_a, n_b, median_a, median_b, effect, p, q), `cohort_models` (model,
layer, feature, n, term, focus, estimate, se, t, p, q), `cohort_models_overview`, `cohort_permanova`,
`cohort_dispersion`, `cohort_ordination` and `cohort_mag_prevalence`. The settings are stored in the manifest under
`cohort`.

## sRNA tables

Written by [`imp-explorer srna`](srna.md) (optional): `srna_samples`, `srna_qc`, `srna_features` (`sample_id`,
`omic` = sRNA fraction, `set`, `feature`, `reads`, `rel`), `srna_tests`, `srna_models`, `srna_models_overview`,
`srna_permanova`, `srna_dispersion`, `srna_ordination`, `srna_covariates`; settings in the manifest under `srna`.
The app combines the cohort results of both as views `all_tests`, `all_models`, ... (`UNION ALL BY NAME`).

## manifest.json

```json
{"release": "release_v0.1", "study": "My study", "created": "2026-10-31T14:54:45", "impexplorer": "0.1.0",
 "tables": {"samples": {"rows": 180, "file": "samples.parquet", "doc": "..."}, "...": {}},
 "sources": {"pipeline": "IMP3", "function_dbs": ["kegg_ko", "..."]},
 "failed_samples": {}, "metadata_without_output": []}
```

## Metabolomics tables

Written by [`imp-explorer metabolomics`](metabolomics.md) (optional): `metab_assays` (`assay`, `label`, `platform`,
`matrix`, `acquisition_date`, `value_type`, `units`, `n_samples`, `n_features`, `n_in_release`, ...), `metab_features`
(`assay`, `feature_id`, `name`, `annotation`, ...), `metab_values` (`sample_id`, `assay`, `feature_id`, `value`,
`source_id`), `metab_samples`; summary in the manifest under `metabolomics`. The cohort analysis treats each assay as a
layer `metab.<assay>`.
