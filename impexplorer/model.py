"""Release data model of the IMP explorer.

Every table is long format and joins on the same keys, so a study portal can add its own layers
(proteomics, metabolomics, immune, sRNA) without changing the IMP core:

  sample_id   one stool / environmental sample (the IMP3 run)
  omic        mg (metagenome), mt (metatranscriptome); optional layers add their own (srna_int, srna_ext, ...)
  feature     taxon, bin, KO, module, ... (with its namespace: rank or db)

Tables (one file each in the release directory, Parquet or TSV.gz):
"""
TABLES = {
    "samples": {
        "doc": "one row per sample: id, IMP3 directory, which omics exist, plus every metadata column given",
        "keys": ["sample_id"],
        "columns": ["sample_id", "imp_dir", "has_mg", "has_mt"],
    },
    "qc": {
        "doc": "per sample and omic: read numbers per preprocessing step, read quality, mapping to the assembly, "
               "Nonpareil diversity, assembly statistics (omic = mgmt)",
        "keys": ["sample_id", "omic", "metric"],
        "columns": ["sample_id", "omic", "metric", "value"],
    },
    "taxonomy": {
        "doc": "read-based taxonomic profiles: Kraken2/Bracken per rank and mOTUs, reads and relative abundance",
        "keys": ["sample_id", "omic", "source", "rank", "taxon"],
        "columns": ["sample_id", "omic", "source", "rank", "taxon", "taxid", "reads", "rel_abundance"],
    },
    "bins": {
        "doc": "MAGs / bins of each sample (MetaWrap refined): CheckM2 quality, MIMAG tier, GTDB-Tk taxonomy, GRiD",
        "keys": ["bin_id"],
        "columns": ["bin_id", "sample_id", "bin", "completeness", "contamination", "quality", "contig_n50",
                    "genome_size", "gc", "domain", "phylum", "class", "order", "family", "genus", "species",
                    "classification", "closest_genome_reference", "closest_genome_ani", "grid", "coverage"],
    },
    "bin_abundance": {
        "doc": "reads of each omic on the contigs of each bin; share of the reads mapped to the assembly and "
               "reads per kb per million mapped (RPKM); mt/mg ratio is computed in the app",
        "keys": ["sample_id", "bin_id", "omic"],
        "columns": ["sample_id", "bin_id", "omic", "reads", "length_bp", "mean_depth", "share_mapped", "rpkm"],
    },
    "functions": {
        "doc": "reads per functional term (featureCounts on the Mantis annotation; paired + single reads), "
               "per omic and database (kegg_ko, kegg_module, kegg_pathway, cazy, enzyme_ec, ...)",
        "keys": ["sample_id", "omic", "db", "feature"],
        "columns": ["sample_id", "omic", "db", "feature", "reads"],
    },
}

MIMAG = (("high", 90, 5), ("medium", 50, 10))      # (tier, min completeness, max contamination)
GTDB_RANKS = ("domain", "phylum", "class", "order", "family", "genus", "species")
KRAKEN_RANKS = {"D": "domain", "P": "phylum", "C": "class", "O": "order", "F": "family", "G": "genus",
                "S": "species"}


def mimag_tier(completeness, contamination):
    try:
        c, x = float(completeness), float(contamination)
    except (TypeError, ValueError):
        return "unknown"
    for tier, cmin, xmax in MIMAG:
        if c >= cmin and x <= xmax:
            return tier
    return "low"
