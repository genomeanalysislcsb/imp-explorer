# IMP explorer

[![tests](https://github.com/genomeanalysislcsb/imp-explorer/actions/workflows/tests.yml/badge.svg)](https://github.com/genomeanalysislcsb/imp-explorer/actions/workflows/tests.yml)
[![docs](https://img.shields.io/badge/docs-genomeanalysislcsb.github.io%2Fimp--explorer-blue)](https://genomeanalysislcsb.github.io/imp-explorer/)

**Documentation: <https://genomeanalysislcsb.github.io/imp-explorer/>**
Changes per version: [CHANGELOG.md](CHANGELOG.md).

Interactive exploration of [IMP3](https://github.com/genomeanalysislcsb/IMP3) results
(metagenomics + metatranscriptomics), with sRNA-IMP, metaproteomics and metabolomics results as optional layers. Stand-alone tool, and the core of
study portals that add their own layers and pages.

Two steps, plus an optional cohort analysis:

1. **export** (on the HPC, needs only pandas): IMP3 sample directories → one *release* directory of tables
   (Parquet with pyarrow, else TSV.gz) + `manifest.json`.
2. **run** (laptop / VM / container): the Dash app on a release directory.
3. **cohort** (optional, HPC or laptop): cohort statistics of the release — QC and covariates, diversity, group
   comparisons within blocks (site), covariate-adjusted models (OLS on CLR, HC3), MAGs, mt/mg ratios — as TSVs, a
   PDF/HTML report and tables for the app's Cohort page. Same statistics and settings as the sRNA-IMP cohort analysis.

```bash
micromamba create -y -p ./env -f environment.yml && micromamba activate ./env
pip install --no-deps -e .  # or: pip install -e ".[app]"

# all sample folders below an IMP3 output root; folder name = sample id, metadata joined by sample_id
imp-explorer export --imp-root /path/to/IMP3_outputs --metadata /path/to/metadata.tsv \
    --study "My study" --out /path/to/release_v0.1 --threads 8
# folder names that differ from the metadata ids: translate them (see the docs)
imp-explorer export --imp-root /path/to/IMP3_outputs --id-regex 'run_(\w+)' --id-case upper ...
# or an explicit sheet: sample_id, imp_dir, (metadata columns ...)
imp-explorer export --samples /path/to/samples.tsv --out /path/to/release_v0.1

# sample metadata: dictionary + template, check before export / cohort (types, units, categories, ids)
imp-explorer metadata template --out metadata/ --data /path/to/release_v0.1
imp-explorer metadata check --metadata clinical.tsv --dictionary metadata/metadata.dictionary.tsv --data /path/to/release_v0.1

# optional: cohort statistics (group comparisons within blocks, adjusted models, report, app page)
imp-explorer cohort --data /path/to/release_v0.1 --block site --model "group+site, group+site+sex+age"
# ... or with the design in a settings file of the study (template: --write-config; an sRNA-IMP config works too)
imp-explorer cohort --data /path/to/release_v0.1 --config study.cohort.yaml

# optional: small-RNA results of sRNA-IMP (its Cohort folder) as an extra layer
imp-explorer srna --data /path/to/release_v0.1 --srna-cohort /path/to/srna_outdir/Cohort
# optional: metaproteomics (taxa / KEGG per fraction, protein catalogue) next to mg and mt
imp-explorer proteomics --data /path/to/release_v0.1 --prot /path/to/proteomics_export
# optional: metabolomics (generic assay / feature / value tables) as an extra layer
imp-explorer metabolomics --data /path/to/release_v0.1 --metab /path/to/metabolomics_export

imp-explorer run --data /path/to/release_v0.1
```

## The app

`imp-explorer run --data <release> [--port 8050] [--host 0.0.0.0]`, then open http://127.0.0.1:8050. On the HPC:
start it on a node and open it through an SSH tunnel (`ssh -L 8050:<node>:8050 <login>`).

| page | what it shows |
|---|---|
| Overview | samples per group, MAGs per sample (MIMAG tiers), sample table with metadata and key numbers (plus key numbers of the sRNA, proteomics and metabolomics layers when present) |
| QC | any QC metric by group (curated metrics first, all IMP3 summary statistics), metric vs metric |
| Taxonomy | Kraken2/Bracken or mOTUs per omic and rank (and the metaproteome, with an across-omics view of one taxon): composition (stacked bars or clustered heatmap), Bray-Curtis PCoA + PERMANOVA, group tests (Kruskal-Wallis / Mann-Whitney, BH-FDR), one taxon by group and mg vs mt |
| MAGs | CheckM2 quality, GTDB taxonomy of the MAGs, abundance (mg) and activity (mt) of a taxon's MAGs by group, log2(mt/mg) per MAG, MAG catalogue |
| Functions | KEGG KO / module / pathway, CAZy, EC per omic (mg, mt and the metaproteome fractions): clustered heatmap of the most variable terms (rows / samples clustered, CLR or log, site colour bar), group tests, one term by group, expression ratio mt/mg |
| Search | one name across all layers (Kraken2, mOTUs, MAGs, functions, proteins, sRNA, metabolites) by group |
| sRNA | sRNA-IMP results (optional): composition per fraction, QC, RNA classes, one feature next to mg and mt |
| Proteomics | metaproteomics (optional, `imp-explorer proteomics`): library QC, which genera make the proteins of a KEGG term, a genus's main KEGG terms, the proteins of a term |
| Metabolomics | metabolomics assays (optional, `imp-explorer metabolomics`): clustered heatmap per assay, one metabolite by group in the assay's units, its correlation with taxa / functions of the same samples, feature table |
| Cohort | results of `imp-explorer cohort` (and of the sRNA-IMP cohort analysis; proteins and metabolomics as layers `prot.*` / `metab.*`): group comparisons within blocks, adjusted models (volcano + table per model, layer, term), beta diversity, QC and covariates |

The top bar sets the grouping variable (any metadata column with 2–12 values) and a sample filter for all pages, e.g.
`site = site_1` to compare groups within one site.

**Plug-ins.** A study portal adds pages without touching the core: a module with `PAGES = [Page(...)]`
(`impexplorer.app.ui.Page`: path, name, order, layout(data), register(app, data), required tables), loaded with
`--plugin my_portal.pages` (or the entry-point group `impexplorer.pages`); extra release tables with `--table name`.
A page only appears when its tables are in the release, so e.g. the sRNA, Proteomics and Metabolomics pages show up only for studies with these data.

## Release tables

| table | rows | content |
|---|---|---|
| `samples` | sample | id, IMP3 directory, `has_mg` / `has_mt`, all metadata columns |
| `qc` | sample × omic × metric | reads per preprocessing step, R1 quality, mapping to the assembly, Nonpareil, assembly stats (`omic = mgmt`) |
| `taxonomy` | sample × omic × source × rank × taxon | Kraken2/Bracken (domain … species) and mOTUs, reads + relative abundance |
| `bins` | bin | CheckM2 completeness / contamination, MIMAG tier, N50, size, GC, GTDB-Tk ranks, closest reference, GRiD |
| `bin_abundance` | sample × bin × omic | reads on the bin's contigs, mean depth, share of mapped reads, RPKM |
| `functions` | sample × omic × db × feature | featureCounts on the Mantis annotation (default kegg_ko, kegg_module, kegg_pathway, cazy, enzyme_ec) |
| `cohort_*` | (optional) | results of `imp-explorer cohort`: tests, models, PERMANOVA, ordination, QC, covariates (see the docs) |

Keys: `sample_id`, `omic`, feature (`taxon` + `rank`, `bin_id`, `feature` + `db`). Study layers (sRNA, proteomics,
metabolomics, immune) add tables with the same keys.

## IMP3 files read

| table | files (per sample directory) |
|---|---|
| qc | `Stats/all_stats.tsv`, `Stats/{mg,mt}/{omic}.read_stats.tsv`, `*_nonpareil.diversity.txt`, `Stats/{omic}/mgmt.assembly.contig_flagstat.txt`, `Stats/mgmt/mgmt.assembly.length.txt` |
| taxonomy | `Analysis/taxonomy/bracken/{omic}/{rank}.bracken` (IMP3 ≥ 3.0) or `Analysis/taxonomy/kraken/{omic}.kraken_bracken_{rank}.report`; `Analysis/taxonomy/mOTUs/{omic}.mOTU.counts.tsv` |
| bins | `Stats/all_bin_stats.tsv` (or `Stats/mgmt/mgmt.bins.tsv`) |
| bin_abundance | `Binning/MetaWrap/bins.contigs`, `Stats/{omic}/mgmt.assembly.contig_read_count.txt`, `…contig_depth.txt` |
| functions | `Analysis/annotation/{omic}.mantis.{db}.counts.tsv` |

All files may be gzipped; missing files give empty tables (partly finished samples still load).

## Tests

```bash
python3 tests/run_tests.py      # synthetic IMP3 tree in the real file formats (tests/make_fixture.py)
```

## Licence

MIT, LCSB BioCore, University of Luxembourg. Related: [IMP3](https://github.com/genomeanalysislcsb/IMP3),
[sRNA-IMP](https://github.com/genomeanalysislcsb/sRNA-IMP).
