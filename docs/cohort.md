# Cohort analysis

`imp-explorer cohort` runs the cohort statistics of a study on a release. It covers QC and covariates, taxa and
diversity, group comparisons within blocks (e.g. site), covariate-adjusted models, MAGs, functions and
metatranscriptome / metagenome ratios.

It writes:

- TSV tables;
- a PDF and HTML report;
- `cohort_*` tables in the release, which the app shows on its **Cohort** page.

The statistics and settings are the same as in the [sRNA-IMP](https://github.com/genomeanalysislcsb/sRNA-IMP) cohort
analysis. The design of a study (groups, blocks, models) can live in one [settings file](#settings-file) that serves
every data type of the study; for studies with sRNA data the `COHORT_*` block of the sRNA-IMP config can be that file.

```bash
imp-explorer cohort --data /path/to/release_v0.1 --metadata /path/to/metadata.tsv \
    --group-col group --group-order control,case_a,case_b --block site \
    --model "group+site, group+site+sex+age, dose+sex+age | group=case_a" \
    --model-notes "M3=dose response in case_a"

# the design in a settings file of the study
imp-explorer cohort --write-config study.cohort.yaml --block site    # template, then edit it
imp-explorer cohort --data /path/to/release_v0.1 --config study.cohort.yaml
```

It needs numpy, pandas and matplotlib (for the report); scipy only makes it faster. For about 200 samples and the
default layers it takes 2–4 minutes with 999 permutations.

## Design: groups and blocks

`--group-col` holds the groups and `--reference` names the reference group (default: the first group of
`--group-order`). `--block` is a column such as site, centre or batch.

With a block, every group is compared with the reference group **within the blocks they share**, and the reference
group split by block shows how large block differences are without disease.

!!! example "Nested design"
    Case group A was sampled only at site 1, case group B only at site 2, controls at both. The comparisons are:

    - A vs control (site 1)
    - B vs control (site 2)
    - control: site 2 vs site 1

    A difference between A and its controls that is about as large as the site contrast in the controls is likely a
    site effect.

PERMANOVA permutes within blocks, and is also run as group after block and as block after group (partial R²). In the
adjusted models the block is an ordinary term (`group+site`).

## Settings file

`--config FILE` reads the settings of a study from a file; options on the command line override it. Two formats are
understood, also mixed:

- `key: value` lines, as written by `--write-config FILE` (a commented template with every key and the values given
  on the command line):

    ```yaml
    name: My study
    metadata: /path/to/metadata.tsv
    group_col: group
    group_order: control,case_a,case_b
    block: site
    model: group+site, group+site+sex+age, dose+sex+age | group=case_a
    model_notes: M3=dose response in case_a
    permutations: 999
    ```

- bash assignments `COHORT_GROUP_COL="group"` or `COHORT_GROUP_COL="${COHORT_GROUP_COL-group}"`, e.g. the cohort
  block of an sRNA-IMP config (`srna_pipeline.config.sh`). Values set in the environment are used as in bash; the file
  is read, not executed.

Keys are the option names below, in any case, with or without the prefix `COHORT_` (`block`, `BLOCK`,
`COHORT_BLOCK`, also the old `COHORT_MATCH_COL`). Keys that only matter for sRNA (`COHORT_MIN_MIRNA_READS`,
`COHORT_NMDS_STARTS`, `COHORT_DIR`) and other variables of a pipeline config are ignored; unknown cohort keys are
reported. `--srna-config` is an alias of `--config`.

A study with several data types keeps one settings file (or the sRNA-IMP config) and passes it to each cohort
analysis, so groups, reference, blocks and models are the same for all of them.

## Settings

| option | key in the settings file | default | meaning |
|---|---|---|---|
| `--metadata` | `metadata` | – | extra metadata (TSV/CSV) merged on the sample id. Ids match exactly, case-insensitively, or before the first `.` (`S01.t0` → `S01`). |
| `--id-col` | `id_col` | first column | id column of `--metadata` |
| `--group-col` | `group_col` | `group` | group column |
| `--group-order` | `group_order` | sorted | order of the groups |
| `--reference` | `reference` | first group | reference group |
| `--block` | `block` | – | block column (site, centre, batch) |
| `--covariates` | `covariates` | all columns | covariates checked for balance |
| `--model` | `model` | `group(+block)` | adjusted models, see below |
| `--model-notes` | `model_notes` | – | question per model: `M3=dose response; M4=...` |
| `--permutations` | `permutations` | 999 | PERMANOVA permutations |
| `--sig-q` | `sig_q` | 0.1 | FDR threshold for highlights and key findings |
| `--outlier-z` | `outlier_z` | 4 | robust z of a QC metric that gives WARN (0 = off) |
| `--qc-drop` | `qc_drop` | `FAIL` | samples left out: `FAIL`, `WARN` (= WARN + FAIL) or `none` |
| `--term-names` | `term_names` | – | `term<TAB>name` TSVs or KEGG flat files for readable KEGG labels |
| `--name` | `name` | study name | report title |
| `--min-mg-pairs` | `min_mg_pairs` | 1e6 | FAIL below this many trimmed metagenome read pairs |
| `--min-mt-pairs` | `min_mt_pairs` | 1e6 | below this, the sample is left out of the metatranscriptome layers (WARN) |
| `--prevalence` | `prevalence` | 0.1 | features with ≥ 10 reads in at least this share of the samples are tested |
| `--layers` | `layers` | see below | feature sets |
| `--nmds` | `nmds` | `kraken.genus.mg,kegg_pathway.mg` | layers with an NMDS |
| `--out` | – | `<release>/cohort` | output directory |
| `--no-release-tables`, `--no-report` | – | | skip the `cohort_*` release tables / the report |

Also: `--config FILE` (settings file), `--write-config FILE` (template), `--data` (release), `--out`.

## Layers

A layer is one feature set: samples × features.

| layer | content |
|---|---|
| `kraken.<rank>.<omic>` | Kraken2/Bracken reads per taxon (rank: phylum … species) |
| `motus.<rank>.<omic>` | mOTUs reads per taxon |
| `mag.<rank>.<omic>` | reads on the MAGs, summed per GTDB taxon of the MAGs |
| `kegg_ko`, `kegg_module`, `kegg_pathway`, `cazy`, `enzyme_ec` `.<omic>` | featureCounts reads per functional term. KEGG global maps (ko011xx, ko012xx) are left out; Mantis `mapNNNNN` duplicates of `koNNNNN` are dropped. |
| `ratio.<layer>` | expression relative to gene copies: CLR(mt) − CLR(mg) in log2 units, per feature tested in both omics, e.g. `ratio.kegg_module`, `ratio.mag.genus` |
| `prot.<rank or db>.<fraction>` | metaproteomics imported with [`imp-explorer proteomics`](proteomics.md): AUC per taxon or KEGG term of one fraction (tests on shares, models on CLR, Bray-Curtis); `prot` = genus and KEGG KO / module / pathway of every fraction, included by default when the release has them |
| `metab.<assay>` | metabolomics assays imported with [`imp-explorer metabolomics`](metabolomics.md): values per feature (tests on the values, models on log2, Euclidean distances of autoscaled log2 values); `metab` = all assays, included by default when the release has them |

`<omic>` is `mg` (metagenome: present, potential) or `mt` (metatranscriptome: active, expressed). Default:
Kraken2 genus and species, mOTUs species, MAG genus, KEGG module, pathway and KO, and CAZy, each for mg and mt, plus
the ratios of MAG genera, KEGG modules and KEGG pathways.

## Analyses

| report section | what |
|---|---|
| Samples and design | samples per group and block; QC by group, with χ² of exclusion vs group |
| Key findings | results below the thresholds, by area (a guide to the report, not a result list) |
| Covariates by group | Kruskal-Wallis / χ² overall and for every comparison |
| QC | QC rules (see below); depth, host and rRNA share, mapping, Nonpareil, assembly, MAGs by group and block |
| Alpha diversity | Shannon and richness of Kraken2 genus and species and of mOTUs species |
| Metrics by group | every per-sample metric: van Elteren within blocks (effect = P(group > reference)), the block contrast, Kruskal-Wallis |
| Taxonomic composition | mean relative abundance of the top taxa per group and block (Kraken2, MAGs; mg and mt) |
| Beta diversity, NMDS | Bray-Curtis PCoA / NMDS; PERMANOVA (permutations within blocks, partial R²) for all groups and every comparison |
| Group comparisons per feature | the same tests per feature of every layer, on relative abundances; BH-FDR per layer and comparison |
| Adjusted models | OLS of CLR abundances (ratio layers: log2 ratio; metrics: own scale) with HC3 robust SE; volcano per model, hits of the term of interest |
| MAGs | MAGs per sample; share of samples with a MAG of a genus by group, van Elteren on recovery |
| Functions and MAGs: mt / mg | median expression ratios of the most variable features per group; group comparisons of the ratios |
| Heterogeneity | dispersion: distance to the own centroid (PCoA, all axes), ratio and Mann-Whitney p per comparison |
| Settings and notes | all settings, layers and notes of the run |

### QC rules

| level | rule | consequence |
|---|---|---|
| FAIL | no metagenome reads, or fewer than `--min-mg-pairs` trimmed pairs | left out of all statistics (`--qc-drop FAIL`) |
| WARN | no metatranscriptome reads, or fewer than `--min-mt-pairs` | left out of the mt and ratio layers only |
| WARN | robust outlier (\|robust z\| > `--outlier-z` and a relevant absolute difference) in mapping, host share, rRNA share, Nonpareil or assembly N50 | kept, listed (left out with `--qc-drop WARN`) |

### Adjusted models

`--model` takes one or more models, separated by commas, with terms joined by `+`. Terms are metadata columns or
per-sample metrics (e.g. `mg_pairs_M` as a depth covariate).

| syntax | example | meaning |
|---|---|---|
| `a+b+c` | `group+site+sex+age` | first term = term of interest; the others adjust |
| `… \| col=v1/v2` | `dose+sex+age \| group=case_a` | fit on the samples with these values only |

Numeric terms are scaled to SD units; categorical terms are treatment dummies against the reference level (group:
the reference group; yes/no variables: `no`). Collinear columns, e.g. a covariate fully confounded with the group,
are dropped and reported.

Features: ≥ 10 reads in ≥ `--prevalence` of the model's samples, CLR with pseudocount 0.5 over these features.

A fit is skipped if it has fewer than 20 samples, fewer than 3 samples in a level, or a sample with leverage > 0.95.

BH-FDR is per model, layer and term. The models are the `M1 … Mn` of the sRNA-IMP cohort analysis when its config is
used.

## Output

```text
<release>/cohort/
├── cohort.report.pdf, cohort.report.html
├── cohort.settings.json            settings, layers, models, notes (also in the release manifest under "cohort")
├── tables/
│   ├── cohort.metadata.tsv         merged metadata
│   ├── cohort.sample_metrics.tsv   per-sample metrics + alpha diversity, group, block
│   ├── cohort.qc_flags.tsv         QC level and reasons per sample
│   └── cohort.<layer>.counts.tsv   feature x sample matrices (ratio layers: .log2ratio.tsv), for other tools
└── stats/
    ├── cohort.design.tsv, cohort.qc_summary.tsv, cohort.covariates.tsv
    ├── cohort.metrics_tests.tsv    metric tests
    ├── cohort.feature_tests.tsv    per layer, feature and comparison: n, medians, effect, p, q
    ├── cohort.permanova.tsv, cohort.dispersion.tsv, cohort.ordination.tsv
    ├── cohort.models.tsv           model, layer, feature, n, term, focus, estimate, se, t, p, q
    ├── cohort.models_overview.tsv  per model and layer: n, features fitted, hits of the term of interest, design notes
    ├── cohort.model_designs.tsv, cohort.mag_prevalence.tsv
    └── cohort.findings.tsv
```

In the release, the tables `cohort_models`, `cohort_models_overview`, `cohort_tests`, `cohort_metric_tests`,
`cohort_permanova`, `cohort_dispersion`, `cohort_ordination`, `cohort_samples`, `cohort_covariates`, `cohort_qc` and
`cohort_mag_prevalence` feed the app's [Cohort page](app.md#cohort). Re-running the command replaces them.

!!! note "What the statistics are and are not"
    The group comparisons are non-parametric and unadjusted. The models adjust for the covariates given, but not for
    everything that differs between groups (diet, medication, stool consistency) unless it is in the model.
    Relative abundances are compositional. Results are candidates to check in the report pages and the app, not
    final effects.
