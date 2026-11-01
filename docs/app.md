# The app

The top bar holds three controls that every page follows:

- **group by**: any metadata column with 2–12 values (default `group` if present). Colours, box plots, group tests
  and PERMANOVA use it.
- **filter** / **values**: keep only samples with the chosen values of one column, e.g. `site = site_1` to compare
  groups within one site.

Tables can be sorted and filtered in their column headers; clicking a row of a group-test table plots that feature.
Plots can be zoomed, and saved as SVG with the camera icon. Feature ids link to their reference page (KEGG, CAZy,
GTDB).

## Overview

Samples per group, MAGs per sample (MIMAG tiers) and a table of all samples with their metadata and key numbers
(trimmed read pairs, mapping to the assembly, assembly N50, Nonpareil diversity, high- and medium-quality MAGs).

![Overview](img/overview.png)

## QC

Any QC metric by group, and any two metrics against each other (Spearman ρ). The curated metrics come first
(read pairs, Q30, mapping to the assembly, Nonpareil, assembly size and N50, number of bins); every statistic of
IMP3's `Stats/all_stats.tsv` is listed after them (prefix `imp:`).

![QC](img/qc.png)

## Taxonomy

Choose a profile: Kraken2/Bracken per rank (with the Kraken2 database set in the IMP3 config, NCBI or GTDB) or mOTUs, for the metagenome or the metatranscriptome.

- **Composition**: stacked bars of the top taxa, samples sorted by group and then by the most abundant taxon.
- **Beta diversity**: Bray-Curtis PCoA of the relative abundances; PERMANOVA (499 permutations) by the grouping.
- **Differences between groups**: Kruskal-Wallis (two groups: Mann-Whitney) per taxon present in ≥ 10 % of the
  samples, Benjamini-Hochberg q, group medians in %.
- **One taxon**: by group (log scale; samples without the taxon at the bottom line) and metagenome vs
  metatranscriptome abundance per sample: points above the diagonal are more active than abundant.

## MAGs

The MAGs of all samples (IMP3: several binners refined by MetaWrap, CheckM2 quality, GTDB-Tk taxonomy).

- **Quality**: completeness vs contamination, coloured by MIMAG tier (high ≥ 90 % / ≤ 5 %, medium ≥ 50 % / ≤ 10 %).
- **Taxonomy of the MAGs**: MAGs per phylum and chosen rank.
- **Abundance and activity**: for one taxon, the summed share of mapped metagenome and metatranscriptome reads on its
  MAGs in each sample, by group; samples without such a MAG count as 0.
- **Activity vs abundance per MAG**: log2(mt share / mg share) of each MAG; > 0 means more transcription than the
  genome's copy number suggests.
- **MAG catalogue**: all bins with quality, size, N50, taxonomy and closest reference genome.

![MAGs](img/mags.png)

## Functions

Reads per functional term from featureCounts on IMP3's Mantis annotation, relative within database and omic.

- **Most variable terms**: z-scores of log10 relative abundance of the 40 most variable terms; term ids link to KEGG.
- **Differences between groups**: as for taxonomy, at the chosen minimum prevalence; a link column opens KEGG / CAZy.
- **One term**: metagenome and metatranscriptome share by group, and the expression ratio log2(mt / mg) per sample.

![Functions](img/functions.png)

## Search

Type part of a name: a taxon at any rank (Kraken2, mOTUs, GTDB names of the MAGs) or a functional id (KO, module,
pathway, CAZyme, EC). Picking a hit plots it in every layer where it occurs, side by side by group, e.g. a genus in
Kraken2 mg and mt, mOTUs mg and mt, and the MAGs mg and mt.

## Cohort

Shown when the release contains the results of [`imp-explorer cohort`](cohort.md). The results are fixed: the global
'group by' and filter controls do not apply.

- **Group comparisons**: per layer and comparison (each group vs the reference within shared blocks, the block
  contrast in the reference group, Kruskal-Wallis across groups), a volcano and the feature table. Clicking a row plots
  the feature by group and block.
- **Adjusted models**: per model (with its question), layer and term (★ = term of interest), a volcano of the
  estimates and the feature table, with design notes. The plot shows the raw values; the model itself uses CLR and the
  covariates.
- **Beta diversity**: PCoA or NMDS of a layer by group and block, PERMANOVA and dispersion.
- **QC and covariates**: covariate balance, metric tests, any metric by group and block, QC flags.

![Cohort](img/cohort.png)

## Statistics: what they are and are not

The tests of the other pages are exploratory: one non-parametric test per feature with BH correction, no covariates,
on the samples and grouping chosen in the top bar. Covariate-adjusted models and comparisons within blocks come from
the [cohort analysis](cohort.md). Relative abundances are compositional; group differences in one taxon can come from
changes in others.
