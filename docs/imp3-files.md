# IMP3 files read

The export reads these files from each IMP3 (3.0) sample folder. Every file may be gzipped. `mgmt` is the assembly
prefix of integrated runs; single-omic runs (`mg` or `mt` assembly) are detected automatically.

| table | files |
|---|---|
| qc | `Stats/all_stats.tsv` · `Stats/{mg,mt}/{omic}.read_stats.tsv` · `Stats/{omic}/{omic}.r{1,2}_nonpareil.diversity.txt` · `Stats/{omic}/mgmt.assembly.contig_flagstat.txt` · `Stats/mgmt/mgmt.assembly.length.txt` |
| taxonomy | `Analysis/taxonomy/bracken/{omic}/{rank}.bracken` (current IMP3) or `Analysis/taxonomy/kraken/{omic}.kraken_bracken_{rank}.report` (earlier runs) · `Analysis/taxonomy/mOTUs/{omic}.mOTU.counts.tsv` (plain or hierarchical `k__…\|p__…` profiles) |
| bins | `Stats/all_bin_stats.tsv` (space-padded cells are fine) or `Stats/mgmt/mgmt.bins.tsv` |
| bin_abundance | `Binning/MetaWrap/bins.contigs` (contig → bin) · `Stats/{omic}/mgmt.assembly.contig_read_count.txt` (hts-nim-tools count-reads: contig, start, end, reads) · `Stats/{omic}/mgmt.assembly.contig_depth.txt` |
| functions | `Analysis/annotation/{omic}.mantis.{db}.counts.tsv` (featureCounts; the paired and single BAM columns are summed) |

Not exported (yet): per-gene counts, mosdepth coverage, SNPs, Bakta annotation, single-binner results, the
Kraken2 contig classification. They are large or per-gene and better queried from the IMP3 folders directly.

## Sample folder layout

For orientation, the main folders of an integrated IMP3 run (sizes of one gut sample with ~35 M metagenome and
~47 M metatranscriptome read pairs):

| folder | content | size |
|---|---|---|
| `Preprocessing/` | trimmed, filtered (PhiX, human), rRNA-depleted reads | ~12 GB |
| `Assembly/` | hybrid assembly, indices, gene calls | ~9 GB |
| `Analysis/annotation/` | Bakta / Mantis annotation, featureCounts tables | ~1.7 GB |
| `Analysis/taxonomy/` | Kraken2, Bracken, mOTUs | ~1.9 GB |
| `Analysis/snps/` | variant calls | ~60 MB |
| `Binning/` | binners (CONCOCT, MetaBAT, MetaDecoder, SemiBin, COMEBin, binny), DAS Tool, MetaWrap | ~0.8 GB |
| `Stats/` | read, mapping, coverage and bin statistics | ~1.3 GB |
| `status/` | `*.done` markers of the IMP3 steps | – |
