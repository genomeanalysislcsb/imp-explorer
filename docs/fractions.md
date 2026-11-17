# Fractions (e.g. extracellular)

A study can sequence a second fraction of the same samples, e.g. the extracellular DNA and RNA next to the
intracellular material of the main IMP3 runs. `imp-explorer fraction` adds it to a release as omics of their own,
so every page, the QC table and the cohort analysis treat it like mg and mt.

```bash
# extracellular DNA: its own IMP3 run (mg only) + its reads remapped to the intracellular assembly
imp-explorer fraction --data release_v0.3 --imp-root /path/IMP/imp_ext_mg --suffix ext
# extracellular RNA: remapped reads only (e.g. single-end reads that IMP3 could not assemble)
imp-explorer fraction --data release_v0.3 --imp-root /path/IMP/imp_ext_mt --suffix ext --remap-only
imp-explorer qc --data release_v0.3 [--flags flags.tsv]     # QC layers for the new omics
imp-explorer cohort --data release_v0.3 ...                 # their layers are added to the default layers
```

The folders below `--imp-root` are named like the release's sample ids (or translated with `--id-map`,
`--id-regex`, `--id-case`); other folders (e.g. `contigs/`, `ids/`) are skipped.

| what the folder holds | omic in the release | tables |
|---|---|---|
| an IMP3 run of the fraction (`Stats/`, `Analysis/`, `Binning/`) | `mg_ext` (and `mt_ext`) | qc, taxonomy (Kraken2 / Bracken, mOTUs), bins (ids `<sample>:ext:<bin>`), bin_abundance, functions |
| `remapping/annotation/<omic>.mantis.<db>.counts.tsv` (featureCounts on the main assembly's annotation) | `mg_ext_remap`, `mt_ext_remap` | functions |

- Rows of these omics already in the release are replaced, so the command can be run again.
- The samples table gets `has_mg_ext`, `has_mt_ext_remap`, ... (not used as metadata); the manifest records the
  fractions under `fractions` with their labels (*extracellular metagenome*, *extracellular metatranscriptome (reads
  on the main assembly)*), which the app and the cohort report use.
- Cohort layers: `kraken.genus.mg_ext`, `kraken.species.mg_ext`, `mag.genus.mg_ext` and the function layers for an
  IMP3 run; function layers only for remapped reads (`kegg_ko.mt_ext_remap`, `pfam.mt_ext_remap`, ...).
- A fraction that is not usable (e.g. degraded RNA) stays in the release and can be flagged for every sample in the
  [QC flags file](sample-qc.md) (`layer` = `mt_ext_remap`), so it is shown but left out with the QC switch.
