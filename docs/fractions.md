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

## Kraken2 database reruns

When the taxonomy was rerun with another Kraken2 database (so that all fractions use the same one), export and add
the fractions with that database:

```bash
imp-explorer export   --imp-root /path/IMP/intra_MGMT --kraken-db k2_gtdb_genome_reps_20250609 ...
imp-explorer fraction --data release --imp-root /path/IMP/imp_ext_mg --suffix ext
imp-explorer fraction --data release --imp-root /path/IMP/imp_ext_mt --suffix ext --remap-only
```

Per sample and omic the profiles are read from
`Analysis/taxonomy/bracken_<db>/<omic>/<omic>.combined.<rank>.bracken` (paired and single reads together; ranks
R1 = domain, P, C, O, F, G, S). With `--kraken-db` only this database is read: a sample without the rerun has no
Kraken taxa rather than those of the other database. `fraction` uses the release's database unless `--kraken-db` is
given, and reads it also for a fraction without an IMP3 run (`--remap-only`): the single-end extracellular mt then
gets taxa as `mt_ext` next to its remapped functions (`mt_ext_remap`).

## Contig taxonomy (MMseqs2, GTDB)

The contigs of an assembly classified with `mmseqs easy-taxonomy` against GTDB give a taxonomic profile of their own:
per rank, the reads mapped to the contigs assigned at that rank or below. Pass the LCA tables as a file pattern with
`{folder}` (the sample folder), `{sample}` and `{omic}`:

```bash
imp-explorer fraction --data release --imp-root /path/IMP/imp_ext_mg --suffix ext \
  --contig-taxonomy '/path/IMP/imp_ext_mg/contigs/{folder}ext/{folder}ext.{omic}.GTDB_lca.tsv'
```

Contig names are matched to the read counts exactly, or by their last `_` field (`SBHC01_merged_C0000100779` vs
`C0000100779`). The profile is taxonomy source `mmseqs2_gtdb`, cohort layers `contig.genus.mg_ext` etc.; `export`
takes the same option for the main assemblies. Without `{omic}` in the pattern one table serves both omics: the mg and
the mt reads on the (co-)assembly give the profiles `contig.genus.mg` and `contig.genus.mt`. Per sample the log shows how
many contigs were matched to read counts; check it on the first run.
