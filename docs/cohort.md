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
`COHORT_DIR`) and other variables of a pipeline config are ignored; unknown cohort keys are
reported. `--srna-config` is an alias of `--config`.

A study with several data types keeps one settings file (or the sRNA-IMP config) and passes it to each cohort
analysis, so groups, reference, blocks and models are the same for all of them.

## Settings

| option | key in the settings file | default | meaning |
|---|---|---|---|
| `--metadata` | `metadata` | – | extra metadata (TSV/CSV) merged on the sample id. Ids match exactly, case-insensitively, or before the first `.` (`S01.t0` → `S01`). |
| `--dictionary` | `dictionary` | – | metadata dictionary ([Sample metadata](metadata.md)): its types decide numeric / categorical columns, its missing codes become NA |
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
| `--term-names` | `term_names` | – | `term<TAB>name` TSVs, KEGG flat files, `Pfam-A.clans.tsv(.gz)` or NCBI `cog-24.def.tab` for readable labels |
| `--name` | `name` | study name | report title |
| `--min-mg-pairs` | `min_mg_pairs` | 1e6 | FAIL below this many trimmed metagenome read pairs |
| `--min-mt-pairs` | `min_mt_pairs` | 1e6 | below this, the sample is left out of the metatranscriptome layers (WARN) |
| `--prevalence` | `prevalence` | 0.1 | features with ≥ 10 reads in at least this share of the samples are tested |
| `--layers` | `layers` | see below | feature sets |
| `--nmds` | `nmds` | `kraken.genus.mg,kegg_pathway.mg` | layers with an NMDS besides the PCoA; `all` = every layer (about 30 s per layer for 180 samples with 10 starts). sRNA-IMP names (`COHORT_NMDS=genus,fn_kegg_pathway,...`) mean the metagenome layers (`kraken.genus.mg`, `kegg_pathway.mg`); sRNA-only sets are left out |
| `--nmds-starts` | `nmds_starts` | 5 | random starts per NMDS (+ 1 PCoA start); `COHORT_NMDS_STARTS` of the sRNA-IMP config |
| `--threads` | – | 1 | processes for the beta-diversity step (PCoA, NMDS, PERMANOVA, dispersion per layer); the results are the same as with 1. With many layers and `--nmds all`, use the cores of the node, e.g. `--threads 16` |
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
| `kegg_ko`, `kegg_module`, `kegg_pathway`, `cazy`, `enzyme_ec`, `pfam`, `cog` `.<omic>` | featureCounts reads per functional term. KEGG global maps (ko011xx, ko012xx) are left out; Mantis `mapNNNNN` duplicates of `koNNNNN` are dropped. |
| `ratio.<layer>` | expression relative to gene copies: CLR(mt) − CLR(mg) in log2 units, per feature tested in both omics, e.g. `ratio.kegg_module`, `ratio.mag.genus`. Missing where the feature has < 10 metagenome reads in a sample (no gene copies, no ratio); a feature is compared only on the samples with a ratio, and only with ≥ `--ratio-min-n` (10) of them on each side |
| `prot.<rank or db>.<fraction>` | metaproteomics imported with [`imp-explorer proteomics`](proteomics.md): AUC per taxon or KEGG term of one fraction (tests on shares, models on CLR, Bray-Curtis); `prot` = genus and KEGG KO / module / pathway of every fraction, included by default when the release has them |
| `met.reactions.<omic>`, `ratio.met.reactions` | reaction groups of [IMP-metabolic](https://github.com/genomeanalysislcsb/IMP-metabolic) imported with `imp-explorer networks`: reads per reaction group (metaG, metaT) and their ratio; included by default when the release has them. The properties of the community networks are added as metrics `net_<set>_<property>` (reactions, largest component, degree, path length, scope, share in MAGs) |
| `met.taxon.<omic>`, `ratio.met.taxon` | taxon × KEGG pathway reads of IMP-metabolic's `activity` (genes on the MAGs of a genus summed per pathway), imported with `imp-explorer networks --activity`; the ratio layer is the taxon-resolved relative expression; included by default when the release has it |
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

### Paper analysis (`--paper-model`)

A second analysis that follows a paper's Methods: one multivariable linear model per layer, with the same formula
for every layer. It runs next to the default analysis and does not change it.

```bash
imp-explorer cohort --data release_v0.4 \
    --paper-model "group+stool_consistency+laxative+age+BMI+sex+diet+bacterial_count+spikes+depth"
```

| part | what | response / distance |
|---|---|---|
| alpha | Shannon and Pielou of the layer ~ model, ANOVA per term | count and protein layers |
| features | each tested feature ~ model, ANOVA per term, BH per layer and term | count layers: relative abundance (`--paper-transform tss`, default) or CLR (`clr`); proteins and metabolomics: CLR, pseudocount = smallest positive value |
| beta | PERMANOVA with all terms, marginal sums of squares per term and permutations under the reduced model, as `adonis2(by = "margin")` (vegan's default `model = "reduced"`: per term the residuals of the model without it are permuted, i.e. the distance-based Freedman–Lane test) | count layers: Bray–Curtis; proteins and metabolomics: Euclidean distance of CLR values (Aitchison); in addition, for count and protein layers, Jaccard on presence / absence (who is there, not how much) |

- **ANOVA**: sums of squares of each term given all the others (type II). With main effects only, as here, type II
  and type III give the same values. F test on the residual degrees of freedom; `r2` = SS(term) / SS(total), the
  share of the variance explained by the term. Terms that are correlated share variance, so the r2 of the terms
  need not add up to 1 − r2 of the residual.
- **`depth`**: log10 of the total reads (proteins: AUC) of the sample in that layer, as the paper's sequencing-depth
  covariate. Other terms are metadata columns or sample metrics; samples with a missing value in a term are left
  out of that layer (complete cases), at least 20 samples per layer.
- **Samples**: per layer, those flagged FAIL in [`sample_qc`](sample-qc.md) are left out (`--paper-qc WARN`: also
  WARN; `none`: no layer QC).
- Settings file keys: `paper_model`, `paper_transform`, `paper_qc`.

Release tables: `cohort_paper_anova` (layer, analysis = shannon / pielou / feature, feature, term, n, df, ss, r2, F,
p, q), `cohort_paper_permanova` (layer, term, n, distance, df, ss, r2, F, p, q; BH per distance) and
`cohort_paper_variance` (per layer and term: r2 and p of Shannon and of the PERMANOVA, Bray–Curtis / Euclidean and
Jaccard, median r2 over the features, share of features with q <
`--sig-q`); the same as `stats/cohort.paper_*.tsv`.

### Agreement between layers (`--concordance`)

How the omics of the same samples agree, in sets of layers with the same features:

```bash
imp-explorer cohort --data release_v0.4 --concordance auto
imp-explorer cohort --data release_v0.4 \
    --concordance "genus=kraken.genus.mg,kraken.genus.mt,prot.genus.intra,prot.genus.extra,srna.genus.srna_int;
                   pfam=pfam.mg,pfam.mt,prot.pfam.intra,prot.pfam.extra"
```

`auto` = these two sets with the layers the release has; `srna.<set>.<fraction>` is an imported sRNA-IMP profile
(e.g. its rRNA-based genera). Per set, on the samples two layers share (FAIL in a layer's `sample_qc` left out):

| table | what |
|---|---|
| `cohort_concord_presence` | per layer and feature: samples with the feature, prevalence, median relative abundance and median CLR (overlap, sparsity, top features) |
| `cohort_concord_alpha` | Spearman ρ of Shannon diversity between every two layers |
| `cohort_concord_mantel` | Mantel test between every two layers: Pearson r of the Bray–Curtis distances, samples permuted (`--permutations`), one-sided p |
| `cohort_concord_features` | Spearman ρ of each feature between two layers (detected in ≥ `--prevalence` of the shared samples in both) |
| `cohort_concord_within` | per sample and pair of layers: Spearman ρ of its two profiles over the shared features, and the median ρ of its profile with the other samples' profiles in the other layer (individuality) |

BH-FDR over the pairs of a set (alpha, Mantel) and per pair (features). The app shows them on the Concordance page.
Feature names must match between layers (GTDB genus names, Pfam accessions). `auto` adds the extracellular
metagenome (`mg_ext`) and the contig taxonomy when the release has them; the extracellular metatranscriptome is left out.

### Intracellular vs extracellular (`--fraction-pairs`)

Which taxa and functions are enriched in one fraction of the same samples:

```bash
imp-explorer cohort --data release_v0.4 --fraction-pairs auto
imp-explorer cohort --data release_v0.4 \
    --fraction-pairs "kraken.genus.mg,kraken.genus.mg_ext; prot.genus.intra,prot.genus.extra; srna.genus.srna_int,srna.genus.srna_ext"
```

`auto` = mg / mg_ext (Kraken and contig genera, KEGG KO and modules, CAZy, Pfam), proteins intra / extra (genera,
KEGG KO, Pfam) and sRNA-Int / sRNA-Ext, as far as the release has them. Per pair, on the samples with both (FAIL in a
layer's `sample_qc` left out), the features detected in ≥ `--prevalence` of the samples in either fraction; both
profiles are closed to these features and CLR-transformed. Per feature a paired Wilcoxon signed-rank test of
CLR(second) − CLR(first), the matched-pairs rank-biserial correlation r as effect size (−1 every sample higher in the
first fraction, +1 every sample higher in the second; the paired counterpart of Cliff's delta) and the median log2
ratio; BH-FDR per pair. A feature counts as enriched in a fraction with q < `--sig-q` and |r| ≥ `--fraction-min-r`
(default 0.5, i.e. at least three in four samples go the same way): with ~180 paired samples almost any shift is
significant, so q alone marks most features. 0 counts by q only.

| table | what |
|---|---|
| `cohort_fraction_diff` | per pair and feature: n, prevalence and median relative abundance per fraction, median log2 ratio, r, z, p, q |
| `cohort_fraction_summary` | per pair: samples, features tested, features enriched in each fraction (q < `--sig-q` and |r| ≥ `--fraction-min-r`; `_any`: q only) |

The app shows them on the Intra vs extra page.

### Contig taxonomy layers

`contig.<rank>.<omic>` layers (MMseqs2 LCA of the contigs, `--contig-taxonomy` of export / fraction) have a long tail
of rare taxa assigned to single short contigs. `--contig-min-share` (default 1e-5) keeps a taxon in a sample only from
that share of the sample's reads on; the log reports how many taxa each layer keeps. 0 keeps all.

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
