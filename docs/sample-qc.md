# Sample QC per layer

A release keeps every sample. `imp-explorer qc` gives each sample a status in each layer, the app marks the flagged
samples, and a switch in the top bar leaves them out of the plots of that layer. A paper and its explorer can then
show the same sample sets, while readers can still look at all samples.

```bash
imp-explorer qc --data release_v0.4 --list-layers                 # the layers of the release
imp-explorer qc --data release_v0.4 --flags paper_exclusions.tsv  # write the table sample_qc
```

## Layers

One layer per omic and fraction, named as in the release tables:

| layer | from |
| --- | --- |
| `mg`, `mt` | IMP3 export (metagenome, metatranscriptome) |
| `srna_int`, `srna_ext`, ... | [sRNA layer](srna.md), one per sRNA-IMP fraction |
| `mp_extra`, `mp_intra`, ... | [proteomics layer](proteomics.md), one per protein fraction |
| `metab.<assay>` | [metabolomics layer](metabolomics.md), one per assay |

## Status and where it comes from

| status | meaning |
| --- | --- |
| PASS | data present, no flag |
| WARN | kept by default; left out with the switch at "PASS only" |
| FAIL | left out with the switch at "without FAIL" or "PASS only" |
| missing | no data in this layer for the sample |

The flags are collected in this order:

1. **auto**: no data in a layer is `missing`; mg / mt read pairs after trimming below `--min-mg-pairs` /
   `--min-mt-pairs` (default 1 M) are FAIL; outliers of the mg / mt metrics found by
   [`imp-explorer cohort`](cohort.md) (table `cohort_qc`) are WARN.
2. **srna**: the QC verdict of sRNA-IMP for each library (PASS, WARN, FAIL with its rules).
3. **manual**: `--flags FILE`, a TSV or CSV with the columns `sample_id`, `layer`, `status`, `reason`. A manual row
   replaces the automatic status of that sample and layer, so it can also clear a flag with PASS. `layer` `*` means
   every layer of that sample, `sample_id` `*` every sample in that layer (a row for one sample overrides it). A
   sample without data in a layer stays `missing`.

```text
sample_id	layer	status	reason
KPL16	metab.bile_acids	FAIL	measured twice, labels to confirm
BRC27	*	WARN	outlier in the paper's QC
BHC22	mt	FAIL	RNA quality too low (paper)
*	mt_ext_remap	FAIL	extracellular mt: single-end, too shallow to use
*	mt_ext	FAIL	extracellular mt: single-end, too shallow to use
```

Use the manual file for the exclusions of a paper, with the paper's reasons, so the explorer and the paper agree.
Status and reason per sample and layer are written to the table `sample_qc` (`sample_id`, `layer`, `status`,
`reasons`, `source`); the counts per layer and the flags file go to the manifest under `sample_qc`. Run
`imp-explorer qc` again after adding a layer or changing the flags.

## In the app

- **QC switch** (top bar, shown when the release has `sample_qc`): *all samples* (default), *without FAIL*,
  *PASS only*. Every page follows it, per layer: a sample that fails in `mt` is left out of the mt plots but stays
  in the mg plots. Views with several layers (Taxonomy *Across omics*, Functions by omic, Search) filter each layer
  on its own.
- **Overview**: *Data and QC per layer*, a grid of every sample in every layer coloured by status (reasons on
  hover), and the table of flagged samples.

The precomputed [cohort statistics](cohort.md) do not follow the switch; they use the sample set of the cohort
run.
