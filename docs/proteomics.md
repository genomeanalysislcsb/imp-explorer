# Proteomics layer

`imp-explorer proteomics` adds metaproteomics results of the same samples to a release. Proteins identified against
the protein catalogue of the IMP3 assemblies are quantified by their peak area (AUC), per fraction, e.g.
extracellular (fecal water) and intracellular (pellet). Each fraction becomes an omic `mp_<fraction>` next to `mg`
and `mt`.

```bash
imp-explorer proteomics --data /path/to/release_v0.1 --prot /path/to/proteomics_export
# group comparisons, models and PERMANOVA of the protein profiles next to the IMP3 layers:
imp-explorer cohort --data /path/to/release_v0.1 --config study.cohort.yaml
```

## Input

Each table can be `.parquet`, `.tsv` or `.tsv.gz` in the `--prot` folder.

| table | columns | |
|---|---|---|
| `proteomics_taxonomy` | `sample_id`, `fraction`, `rank`, `taxon`, `auc`, `rel`, `n_proteins` | required |
| `proteomics_functions` | `sample_id`, `fraction`, `db` (`kegg_ko`, `kegg_module`, `kegg_pathway`, `pfam`, `cog`), `feature`, `auc`, `rel`, `n_proteins` | required |
| `proteomics_samples` | `sample_id`, `fraction`, `n_proteins`, `n_peptides`, `total_auc` | optional |
| `proteomics_proteins` | `protein_id`, `taxon`, `rank`, `lineage` (GTDB, `d_…;p_…;…;g_…;s_…`), `product`, `kegg_ko`, `kegg_module`, `kegg_pathway`, … (several terms `,`) | optional: catalogue |
| `proteomics_abundance` | `sample_id`, `fraction`, `protein_id`, `auc`, `n_peptides` | optional |

Conventions:

- **`rel`:** the share of the total AUC of the sample and fraction. It is computed from `total_auc`, or from the AUC
  of the table, when missing.
- **Taxonomy:** a protein counts for every rank of its lineage, down to the rank it is assigned to. Taxon names
  follow the GTDB style of IMP3 (`Blautia_A`).
- **Functions:** a protein with several terms of a database counts fully for each of them, as featureCounts does for
  mg and mt.

For Expobiome, `proteomics_export.py` writes these tables from the files of the proteomics analysis.

## Options

| option | default | meaning |
|---|---|---|
| `--data`, `--prot` | – | release; folder with the tables |
| `--no-proteins` | off | leave out the catalogue and the per-protein abundances (smaller release, no taxon × function views) |
| `--id-map`, `--id-regex`, `--id-case` | – | map the sample ids to the release, as for [export](export.md) |
| `--format` | auto | Parquet (with pyarrow) or TSV.gz |

**In the release:** `prot_taxonomy`, `prot_functions`, `prot_samples`, `prot_proteins`, `prot_abundance` (`omic` =
`mp_<fraction>`), and a `proteomics` section in the manifest.

## In the cohort analysis

`imp-explorer cohort` analyses the protein profiles as layers `prot.<rank or database>.<fraction>`, e.g.
`prot.genus.extra` or `prot.kegg_module.intra`. When the release has proteomics, genus, KEGG KO, KEGG module and KEGG
pathway (and Pfam / COG when imported) of every fraction are included by default; with `--layers`, add `prot` for the same set, or single layers.

| step | protein layers |
|---|---|
| features | detected (AUC > 0) in ≥ `--prevalence` of the samples of the fraction |
| group comparisons | on the shares of the layer, the same tests as for the read-based layers |
| adjusted models | CLR of the shares, with half the smallest positive share as pseudocount |
| beta diversity, PERMANOVA | Bray-Curtis of the shares |

Samples without a fraction are not part of that fraction's layers.

## In the app

- **Taxonomy and Functions:** the omics `mp_extra` / `mp_intra` (source "Metaproteomics") are profiles like mg and mt.
  - Taxonomy has an **Across omics** panel: one taxon in Kraken2 mg / mt and in the protein fractions, by group.
  - Functions shows a term's share in mg, mt and the protein fractions side by side.
- **Proteomics page:**
  - library QC (proteins, peptides, total AUC) by group;
  - **who makes it**: the genera behind a KEGG term's protein abundance;
  - **what the genus makes**: a genus's main KEGG terms;
  - the proteins of a term, with taxon and product.
- **Search:** taxa and KEGG terms are found in the protein profiles too.
- **Overview:** a row of key numbers, and one column per fraction in the sample table.
- **Metabolomics page:** a metabolite can be correlated with the protein genera or KEGG modules.
