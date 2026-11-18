# Changelog

## Unreleased

- **Intra vs extra: effect-size threshold** (`cohort --fraction-min-r`, default 0.5; selector on the page): a
  feature counts as enriched with q < `--sig-q` and |rank-biserial r| ≥ the minimum, since with ~180 paired samples
  almost any shift is significant. The page's summary, volcano (dotted lines), strongest differences and table
  (enriched features by default, or all with the column *enriched*) follow it.
- Intra vs extra page: summary plot drawn; volcano point size = median abundance; p and q shown with 3 significant
  digits. A test checks every `%`-formatted string of the package.
- Labels: the metagenome and metatranscriptome of the main export are called *intracellular* once the release has an
  extracellular fraction (as the protein fractions already were).

## 0.4.0 (2026-11)

The analyses of a multi-omics cohort paper (Expobiome, cross-sectional part) in the explorer: QC flags per sample and
layer, extracellular fractions, a paper's linear models, agreement between the layers and intracellular vs
extracellular differences, each with its page in the app.

Data

- **Sample QC per layer** (`imp-explorer qc`, table `sample_qc`): every sample keeps its place in the release and
  gets a status (PASS, WARN, FAIL, missing) in each layer, from automatic rules, the sRNA-IMP verdicts and a manual
  flags file (e.g. a paper's exclusions with their reasons; `sample_id` `*` = every sample of a layer; tabs or
  spaces). The app has a QC switch in the top bar (all samples / without FAIL / PASS only) that every page applies
  per layer, and the Overview shows the sample x layer grid and the flagged samples.
- **Fractions** (`imp-explorer fraction`): a second fraction of the same samples (e.g. extracellular) as omics of its
  own: its IMP3 run (`mg_ext`: taxa, MAGs, functions, QC) and its reads remapped to the main assembly
  (`mt_ext_remap`: functions); QC layers and default cohort layers include them.
- **Kraken2 database reruns** (`--kraken-db <db>` for `export` and `fraction`): read the Bracken profiles of a rerun
  (`Analysis/taxonomy/bracken_<db>/<omic>/<omic>.combined.<rank>.bracken`, domain from `R1`) instead of IMP3's own,
  so all fractions use the same database; only that database is read (no mixing). `fraction` takes the release's
  database by default and also reads the rerun of a fraction without an IMP3 run (`--remap-only`, e.g. the single-end
  extracellular mt -> `mt_ext`). The manifest records the database (`sources.kraken_db`, per fraction) and which
  tables a fraction has, and the cohort adds only those layers.
- **Contig taxonomy** (`--contig-taxonomy PATTERN` for `export` and `fraction`): the MMseqs2 `easy-taxonomy` LCA
  table of each sample's contigs against GTDB (`<prefix>_lca.tsv`), weighted by the reads on the contigs
  (`Stats/<omic>/*.assembly.contig_read_count.txt`), as taxonomy source `mmseqs2_gtdb`; cohort layers
  `contig.<rank>.<omic>`. Names in GTDB spelling (`Blautia_A`), so they match the Kraken2 GTDB profiles; the log
  reports per sample how many contigs matched the read counts.

Analyses (`imp-explorer cohort`)

- **Paper analysis** (`--paper-model`): one multivariable linear model per layer as in a paper's Methods, next to
  the default analysis: type II ANOVA per term for Shannon, Pielou and every feature, variance explained per term,
  PERMANOVA with all terms (marginal R², Bray–Curtis or Euclidean on CLR), the layer's flagged samples left out
  (`sample_qc`). Tables `cohort_paper_anova`, `cohort_paper_permanova`, `cohort_paper_variance`.
- **Agreement between layers** (`--concordance`): overlap and sparsity of the features, top features, Spearman of
  Shannon diversity, Mantel tests and per-feature Spearman between every two layers of a set (e.g. genera in mg, mt,
  mg_ext, contigs, proteins and sRNA; Pfam in mg, mt, mg_ext and proteins; `auto` leaves out the extracellular
  metatranscriptome), and within each sample (profiles of the same sample vs other samples: individuality).
- **Intracellular vs extracellular** (`--fraction-pairs auto|pairs`, tables `cohort_fraction_diff` /
  `cohort_fraction_summary`): per feature a paired Wilcoxon on CLR, matched-pairs rank-biserial r and median log2
  ratio, BH-FDR per pair.
- **Contig taxonomy filter** (`--contig-min-share`, default 1e-5): rare taxa of single short contigs left out of the
  `contig.*` layers.

App

- **Variance page**: variance explained per layer and term (beta diversity, Shannon, features), terms ranked across
  layers, unexplained variance per layer, significant features per term.
- **Concordance page**: UpSet of the layers, top features (percentile within each layer), sparsity, Shannon and
  Mantel heatmaps, same feature in two layers, within each sample, features of one pair.
- **Intra vs extra page**: features enriched per fraction, volcano, strongest differences, table.
- Dropdowns with long names no longer overlap.

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
