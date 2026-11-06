# Changelog

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
