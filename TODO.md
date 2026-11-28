# imp-explorer: to do

Open points that are known but not scheduled. Remove an item when it is done (the commit says so).

## Cohort

- **Read processing, left over**: record the host filter method of each run (IMP3 config `filtering.method`;
  `bwa_any` before 3.1.0) next to the host share and mark host shares of `bwa_any` runs in the report; a rule for an
  extreme metatranscriptome / metagenome library ratio; Deacon's per-step summaries once IMP3 keeps them in the logs.
- **smORFs across samples** (moved here from the IMP3 smORF module): read `imp/<sample>/Analysis/smORF/` of the
  cohort (smORF faa / table / read counts), cluster all smORFs with MMseqs2 at 100 % and 95 % identity and coverage,
  presence per sample and group, core / singleton numbers, length distributions per group, read support (MG only /
  MT only / MG+MT) and metaP status per cluster; add to `cohort.html` and as tables; test cluster presence and
  expression between groups like the other features.
