# IMP explorer

IMP explorer turns the results of [IMP3](https://github.com/genomeanalysislcsb/IMP3) (integrated metagenomics and
metatranscriptomics) for a whole study into an interactive web app: QC, read-based taxonomy, MAGs, functions and a
search across all of them, by group and site.

![Taxonomy page](img/taxonomy.png)

It works in two steps, plus an optional third:

```mermaid
flowchart LR
  imp[(IMP3 sample folders)] -->|imp-explorer export<br/>HPC, pandas only| rel[(release directory<br/>Parquet tables + manifest)]
  meta[(sample metadata)] --> rel
  rel -->|imp-explorer cohort| coh[(cohort tables<br/>+ PDF/HTML report)]
  coh --> rel
  rel -->|imp-explorer run| app[Dash app<br/>laptop / VM / container]
  srna[(sRNA-IMP, proteomics, ...)] -.->|study layers| rel
  plug[study plug-in pages] -.-> app
```

1. **export** reads the IMP3 output folder of every sample and writes one *release*: six long-format tables with shared
   keys (sample, omic, feature) plus a manifest. It needs only pandas and runs where the data are.
2. **cohort** (optional) runs the cohort statistics: group comparisons within blocks, covariate-adjusted models,
   beta diversity, MAGs, mt/mg ratios; TSVs, a PDF/HTML report and tables for the app ([Cohort analysis](cohort.md)).
3. **run** starts the app on a release. The release is small (tens to a few hundred MB for ~200 samples), so the app
   runs on a laptop, a VM or in a container.

The app is the core of study portals: a study adds its own layers (sRNA from [sRNA-IMP](https://github.com/genomeanalysislcsb/sRNA-IMP),
proteomics, metabolomics, immune data) as extra tables and pages, without changing the IMP core
([plug-ins](plugins.md)).

| page | question it answers |
|---|---|
| [Overview](app.md#overview) | what is in the study, how many samples, MAGs, features |
| [QC](app.md#qc) | are depth, quality, mapping and assembly comparable between groups? |
| [Taxonomy](app.md#taxonomy) | who is there (metagenome) and who is active (metatranscriptome)? |
| [MAGs](app.md#mags) | which genomes were recovered, how abundant and how active are they? |
| [Functions](app.md#functions) | which KEGG / CAZy / EC functions differ, and are they expressed? |
| [Search](app.md#search) | one taxon or function across every layer |
| [sRNA](srna.md) | small-RNA results of sRNA-IMP (optional): QC, RNA classes, miRNAs, taxa, next to mg and mt |
| [Cohort](app.md#cohort) | results of the [cohort analysis](cohort.md): group comparisons within blocks, adjusted models, beta diversity, QC and covariates |

## Quick start

```bash
micromamba create -y -p ./env -f environment.yml && micromamba activate ./env
pip install --no-deps -e .

imp-explorer export --imp-root /path/to/IMP3_outputs --metadata /path/to/metadata.tsv \
    --study "My study" --out /path/to/release_v0.1
imp-explorer cohort --data /path/to/release_v0.1 --group-col group --block site \
    --model "group+site, group+site+sex+age"         # optional: cohort statistics + report
imp-explorer run --data /path/to/release_v0.1          # http://127.0.0.1:8050
```

## Citation and licence

MIT licence. LCSB BioCore, University of Luxembourg.
