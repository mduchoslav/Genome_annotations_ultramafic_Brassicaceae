Annotation of ultramafic Brassicaceae genomes
================
Miloš Duchoslav
2026-02

This repository contains scripts and data for annotation of the newly assembled genomes of *Aethionema saxatile*, *Cardamine glauca*, *Erysimum linariifolium*, *Noccaea praecox*, and *Odontarrhena muralis* for the following publication:

> Mahnaz Nezamivand-Chegini, Miloš Duchoslav, Raúl Wijfjes, Gabriela Šrámková, Jana Nosková, Vít Bureš, Tereza Koberová, Terezie Mandáková, Tomica Mišljenović, Ksenija Jakovljević, Panayiotis G. Dimitrakopoulos, Stanislav Španiel, Bruno Huettel, Martin A. Lysak, Levi Yant, Korbinian Schneeberger, Filip Kolář: Chromosome-level reference genome assemblies of five Brassicaceae species inhabiting challenging ultramafic substrates, 2026 (under review)

## Scripts and detailed description of methods

Each species has it's own folder with with some intermediate data (intermediate data that are small enough and valuable enough to share), final results and scripts used.

I use for documentation RMarkdown combining snippets of R code and Bash code. These Rmd files are then knitted to GitHub markdown files for better viewing.

The main files with scripts and other information about annotation (GitHub md file first, original Rmd file second):

- [*Aethionema saxatile*](Aethionema_saxatile_2025_06/annotation_Aethionema_saxatile.md) | [*Aethionema saxatile* (Rmd file)](Aethionema_saxatile_2025_06/annotation_Aethionema_saxatile.rmd)
- [*Cardamine glauca*](Cardamine_glauca_2025_06/annotation_Cardamine_glauca.md) | [*Cardamine glauca* (Rmd file)](Cardamine_glauca_2025_06/annotation_Cardamine_glauca.rmd)
- [*Erysimum linariifolium*](Erysimum_linariifolium_2025_09/annotation_Erysimum_linariifolium.md) | [*Erysimum linariifolium* (Rmd file)](Erysimum_linariifolium_2025_09/annotation_Erysimum_linariifolium.rmd)
- [*Noccaea praecox*](Noccaea_praecox_2024_12/annotation_Noccaea_praecox.md) | [*Noccaea praecox* (Rmd file)](Noccaea_praecox_2024_12/annotation_Noccaea_praecox.rmd)
- [*Odontarrhena muralis*](Odontarrhena_muralis_2025_10/annotation_Odontarrhena_muralis.md) | [*Odontarrhena muralis* (Rmd file)](Odontarrhena_muralis_2025_10/annotation_Odontarrhena_muralis.rmd)

Installation and versions of tools used:  
[GitHub md file](Installation_of_SW.md) | [original RMarkdown file](Installation_of_SW.Rmd)

If you want to get inspired and do such annotation for other species, I recommend using my [script for annotation of *Odontarrhena muralis*](Odontarrhena_muralis_2025_10/annotation_Odontarrhena_muralis.rmd), which I run recently and which is the most developed.

## Final results

### Results for each species

The `final_files` folder for each species contains:

- **Soft-masked reference genome** (`*_masked.fa.gz`)
	- This is version used for annotation, note that the version in GenBank database has different scaffold names.
- **Genome annotation** (`*.gff.gz`)
	- The genome annotation contains protein-coding genes (*gene* feature with *transcript* as a child feature) and non-coding RNA genes (*ncrna_gene* feature with *ncrna* as a child feature). The non-coding RNA genes are derived from assembled transcripts that don't have any overlap with predicted protein-coding genes.
- **Coding sequences (CDS)** for protein-coding genes (`*_cds.fasta.gz`)
- **Protein sequences** for protein-coding genes (`*_proteins.fasta.gz`)
- **Table with external evidence** for protein-coding genes (`*_protein_coding_genes_support.tsv`)
	- Table shows for each gene if there is overlap with transcripts assembled from RNA-seq data used for annotation or with aligned protein sequences from *Arabidopsis thaliana*, *Arabidopsis lyrata* and *Brassica rapa*).
	
Link to final results for each species:
- [*Aethionema saxatile*](Aethionema_saxatile_2025_06/final_files)
- [*Cardamine glauca*](Cardamine_glauca_2025_06/final_files)
- [*Erysimum linariifolium*](Erysimum_linariifolium_2025_09/final_files)
- [*Noccaea praecox*](Noccaea_praecox_2024_12/final_files)
- [*Odontarrhena muralis*](Odontarrhena_muralis_2025_10/final_files)

### Statistics for all species

Multi-species comparison of number of predicted potein-coding genes and their support by external evidence is in table [prot_coding_genes_support.tsv](01_multispecies_stats/prot_coding_genes_support.tsv).

## Orthology and functional annotations

Tables with orthologues in *Arabidopsis thaliana* and other species and functional annotations obtained both through InterProScan and from *Arabidopsis thaliana* orthologues can be found in my repository [Brassicaceae_orthology](https://github.com/mduchoslav/Brassicaceae_orthology).