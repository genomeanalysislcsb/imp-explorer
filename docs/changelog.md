# Changelog

## Unreleased

- **Sample QC per layer** (`imp-explorer qc`, table `sample_qc`): every sample keeps its place in the release and
  gets a status (PASS, WARN, FAIL, missing) in each layer, from automatic rules, the sRNA-IMP verdicts and a manual
  flags file (e.g. a paper's exclusions with their reasons). The app has a QC switch in the top bar (all samples /
  without FAIL / PASS only) that every page applies per layer, and the Overview shows the sample x layer grid and
  the flagged samples.
- **Paper analysis** (`imp-explorer cohort --paper-model`): one multivariable linear model per layer as in a paper's
  Methods, next to the default analysis: type II ANOVA per term for Shannon, Pielou and every feature, variance
  explained per term, PERMANOVA with all terms (marginal R², Bray–Curtis or Euclidean on CLR), the layer's flagged
  samples left out (`sample_qc`). Tables `cohort_paper_anova`, `cohort_paper_permanova`, `cohort_paper_variance`.
- **Fractions** (`imp-explorer fraction`): a second fraction of the same samples (e.g. extracellular) as omics of its
  own: its IMP3 run (`mg_ext`: taxa, MAGs, functions, QC) and its reads remapped to the main assembly
  (`mt_ext_remap`: functions); QC layers, default cohort layers and concordance sets include them.
- **Agreement between layers** (`imp-explorer cohort --concordance`): overlap and sparsity of the features, top
  features by median CLR, Spearman of Shannon diversity, Mantel tests and per-feature Spearman between every two
  layers of a set (e.g. genera in mg, mt, proteins and sRNA; Pfam in mg, mt and proteins); Concordance page.
- **Variance page** (app): variance explained per layer and term (beta diversity, Shannon, features), terms ranked
  across layers, unexplained variance per layer, significant features per term.

## 0.3.0 (2026-11)

Optional data layers next to the IMP3 metagenome (mg) and metatranscriptome (mt), each imported into a release by
its own command and shown in the app and the cohort analysis:

- **sRNA** (`imp-explorer srna`): results of the sRNA-IMP cohort analysis; sRNA page, sRNA in Cohort, Search and
  Overview; sRNA compared with the metaproteome and with metabolites.
- **Metaproteomics** (`imp-explorer proteomics`): taxa and functions per protein fraction (e.g. extracellular /
  intracellular) as `mp_<fraction>` next to mg and mt on the Taxonomy, Functions and Search pages; protein
  catalogue; Proteomics page; cohort layers `prot.*` (CLR, Bray-Curtis).
- **Metabolomics** (`imp-explorer metabolomics`): generic assay / feature / value tables for any platform;
  Metabolomics page, Search and Overview; cohort layers `metab.*` (autoscaled log2, Euclidean).
- **Pfam and COG** for mg / mt and proteins (links, names, cohort layers), besides KEGG.

Cohort analysis (`imp-explorer cohort`):

- study settings file (`--config`, `--write-config`), also the `COHORT_*` block of an sRNA-IMP config;
- NMDS for chosen or all layers (`--nmds all`, `--nmds-starts`), settings saved with the results;
- `--threads`: beta diversity (PCoA, NMDS, PERMANOVA) per layer in parallel processes; PERMANOVA reproducible
  (permutations independent of string hash order);
- `--dictionary`: column types and missing codes from a metadata dictionary.

Sample metadata (`imp-explorer metadata`): dictionary and template, check of a metadata table (ids, types, units,
categories, missing values, group sizes), clean copy.

App: clustered heatmaps on the Taxonomy and Functions pages (colour bar below the heatmap), Cohort page (model list,
feature plot of the top hit), Overview with key numbers of every layer.

Export: worker processes are spawned (Python 3.12), plotly < 6 with dash 2, generic sample ids and id translation
(`--id-map`, `--id-regex`, `--id-case`).

## 0.1.1 (2026-11)

Generic sample ids, paths and examples; id translation tested; mOTUs ranks; docs.

## 0.1.0 (2026-10)

Export of IMP3 results into a release (Parquet / TSV tables + manifest) and the interactive Dash app.
