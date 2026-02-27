Annotation - *Cardamine glauca*
================
Miloš Duchoslav
2025-06-05

- [Introduction](#introduction)
- [Instructions for use for different
  species](#instructions-for-use-for-different-species)
- [Making new folders](#making-new-folders)
- [Preparation for annotation using
  Braker](#preparation-for-annotation-using-braker)
  - [Genome assembly](#genome-assembly)
  - [Preparation of RNAseq data](#preparation-of-rnaseq-data)
  - [Alignment of RNAseq reads with
    HISAT2](#alignment-of-rnaseq-reads-with-hisat2)
  - [Preparation of protein
    sequences](#preparation-of-protein-sequences)
- [Running Braker](#running-braker)
- [Quality check of the Braker
  annotation](#quality-check-of-the-braker-annotation)
  - [Running IGV on Metacentrum](#running-igv-on-metacentrum)
  - [Number of genes](#number-of-genes)
  - [Checking proteins predicted by Braker for stop
    codons](#checking-proteins-predicted-by-braker-for-stop-codons)
  - [Alignment of protein sequences to the genome
    assembly](#alignment-of-protein-sequences-to-the-genome-assembly)
  - [Assembly of transcripts from RNAseq
    alignment](#assembly-of-transcripts-from-rnaseq-alignment)
  - [Extraction of transcript
    sequences](#extraction-of-transcript-sequences)
  - [Busco](#busco)
  - [Intersect between gene models and aligned
    proteins](#intersect-between-gene-models-and-aligned-proteins)
  - [Intersect between gene models and assembled
    transcripts](#intersect-between-gene-models-and-assembled-transcripts)
  - [Getting the list of “reliable”
    genes](#getting-the-list-of-reliable-genes)
  - [Move this old version of intersect
    files](#move-this-old-version-of-intersect-files)
  - [Comparing predicted genes with previous
    run](#comparing-predicted-genes-with-previous-run)
  - [Conversion of Braker GTF to GFF](#conversion-of-braker-gtf-to-gff)
- [Adding unknown expressed
  features](#adding-unknown-expressed-features)
  - [Merging annotations](#merging-annotations)
- [Changing IDs in R](#changing-ids-in-r)
  - [Statistics of annotation using
    AGAT](#statistics-of-annotation-using-agat)
  - [Extraction of protein and coding sequences using
    AGAT](#extraction-of-protein-and-coding-sequences-using-agat)
- [Final files](#final-files)
- [Reliability of the predicted
  genes](#reliability-of-the-predicted-genes)
  - [Intersect between gene models and aligned
    proteins](#intersect-between-gene-models-and-aligned-proteins-1)
  - [Intersect between gene models and assembled
    transcripts](#intersect-between-gene-models-and-assembled-transcripts-1)
  - [Generating table of reliability of
    genes](#generating-table-of-reliability-of-genes)

# Introduction

This RMarkdown file (or its markdown version for GitHub) documents the
annotation of the newly assembled genome of *Cardamine glauca*.

This file includes:

1.  BASH code that I ran at MetaCentrum (Czech national grid
    infrastructure) running PBS scheduling system for batch jobs.
2.  R code that I ran locally.

### SW installation and versions

The SW installation instructions and versions of SW used is described in
[Installation_of_SW.md](../Installation_of_SW.md).

# Instructions for use for different species

If you want to run this for new species, here is a (partial) list of
what you should change: - Change folder for the species
(`Cardamine_glauca_2025_06`) throughout the document. - Change
genome_assembly variable
(`genome_assembly="Cardamine_glauca_CUNI_V1_2024_10_masked.fa"`)
throughout the document. - Change pattern for finding of RNAseq data
(`-iname XC*.gz`). - Consider changing of species with good genomes
selected for validation by alignment of protein sequences. - Change
species in Braker run (`--species=Cardamine_glauca`). - Change prefix
for gene IDs. - Adjust GFF header.

# Making new folders

``` sh
cd /storage/brno12-cerit/home/duchmil/annotations/
# make new folder and subfolders
mkdir Cardamine_glauca_2025_06
cd Cardamine_glauca_2025_06
mkdir genome_assembly
mkdir rnaseq
mkdir metacentrum_scripts
```

# Preparation for annotation using Braker

## Genome assembly

In the previous version of this script, I used assembly named
`$genome_assembly`. This was now renamed to
`Cardamine_glauca_CUNI_V1_2024_10_masked.fa`, but it is the same
assembly.

``` sh
# Copying assembly from Mahnaz
cd /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/genome_assembly
# Assembly copied from google drive through my PC.
```

### Checking assembly statistics

``` sh
cd /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/genome_assembly
genome_assembly="Cardamine_glauca_CUNI_V1_2024_10_masked.fa"

# Current assembly
grep -c '>' $genome_assembly # 53


module load quast
quast.py $genome_assembly
less quast_results/latest/report.txt

module load seqtk
# info about columns
seqtk comp
# info about sequences
seqtk comp $genome_assembly | column -t
```

    Assembly                    Cardamine_glauca_CUNI_V1_2024_10_masked
    # contigs (>= 0 bp)         53
    # contigs (>= 1000 bp)      53
    # contigs (>= 5000 bp)      53
    # contigs (>= 10000 bp)     53
    # contigs (>= 25000 bp)     53
    # contigs (>= 50000 bp)     45
    Total length (>= 0 bp)      301015024
    Total length (>= 1000 bp)   301015024
    Total length (>= 5000 bp)   301015024
    Total length (>= 10000 bp)  301015024
    Total length (>= 25000 bp)  301015024
    Total length (>= 50000 bp)  300693161
    # contigs                   53
    Largest contig              39236745
    Total length                301015024
    GC (%)                      37.85
    N50                         34778859
    N75                         30582462
    L50                         4
    L75                         7
    # N's per 100 kbp           0.80

## Preparation of RNAseq data

### Copying our RNAseq data

One sample was sequenced both by Macrogen and Novogene to compare them.

``` sh
# Copy RNAseq data
cd /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/rnaseq
mkdir 1_raw_reads
cd 1_raw_reads

# find the RNAseq data for Cardamine ("XC" code)
find /storage/brno12-cerit/home/filip_kolar/00_shared_space_ecolgen/raw_sequencing_data/RNA_shortread/ -iname XC*.gz
# Copy these data
find /storage/brno12-cerit/home/filip_kolar/00_shared_space_ecolgen/raw_sequencing_data/RNA_shortread/ -iname XC*.gz -exec cp -v -t . {} \;

# Rename the Macrogen version of sample that I do not have two samples with the same name
mv XC174YL_1.fastq.gz XC174YL_Macrogen_1.fastq.gz
mv XC174YL_2.fastq.gz XC174YL_Macrogen_2.fastq.gz
```

### Script to run both FastQC and MultiQC on RNAseq data

``` sh
#!/bin/bash
#PBS -N FastQC_RNAseq
#PBS -l select=1:ncpus=6:mem=96gb:scratch_local=10gb:brno=True
#PBS -l walltime=2:00:00
#PBS -j oe

trap 'clean_scratch' TERM EXIT

# define a DATADIR variable: directory where the input files are taken from
DATADIR=/storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/rnaseq/1_raw_reads
# directory for output
OUTDIR=/storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/rnaseq/1_raw_reads
# report name
report_name="Cardamine_RNAseq_MultiQC_report"

# create the fastqc dir if it does not exists
if [ ! -d $OUTDIR/fastqc ]; then 
mkdir $OUTDIR/fastqc
fi

# create the multiqc dir if it does not exists
if [ ! -d $OUTDIR/multiqc ]; then
mkdir $OUTDIR/multiqc
fi

# test if scratch directory is set
# if scratch directory is not set, issue error message and exit
#test -n "$SCRATCHDIR" || { echo >&2 "Variable SCRATCHDIR is not set!"; exit 1; }

# load Java
module add openjdk/

# FastQC run

# Version 1 for files in multiple folders
time find $DATADIR -type f \( -iname "*.fastq.gz" -o -iname "*.fq.gz" \) -print0 | xargs -0 /storage/brno12-cerit/home/duchmil/SW/fastqc/FastQC/fastqc -t 6 -o $OUTDIR/fastqc -f fastq

# Version 2 for files in single folder (not going to subfolders)
# time /storage/brno12-cerit/home/duchmil/SW/fastqc/FastQC/fastqc -t 10 -o $OUTDIR/fastqc -f fastq $DATADIR/*.fastq.gz


# Running MultiQC
# activation
source /storage/brno2/home/duchmil/SW/mambaforge/bin/activate Multiqc

# run MultiQC
time multiqc --filename $report_name --outdir $OUTDIR/multiqc $OUTDIR/fastqc

# Resources: 7 min, 24 GB memory, 61 % CPU.
```

``` sh
cd /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/metacentrum_scripts
qsub run_fastqc_multiqc_Milos_RNAseq.sh
```

#### Result of FastQC - comparison of Macrogen and Novogene

Macrogen has higher Quality Scores, but the newer Novogene data are
similar. On the other hand, they have ~ 15 % of reads with adapters,
some overrepresented sequences and some small bumps in “Per Sequence GC
Content”. As Novogene has no reads with adapters, I suspect that they
filtered the results and thus it is not well comparable.

### Trimming

Trimming based on quality, removal of the adaptors. Adaptors are there
only in case that the insert is too short and it is read through. It
means that they are usually at the 3’ terminus.

**Warning: Trim Galore uses the first 1 million sequences of the first
file to detect adapter type. If there are files from different sources
(with different adapters), they should be run separately!**

``` sh
cd /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/rnaseq/
mkdir 2_trimmed_reads
```

Trimming and FastQC and MultiQC after trimming

``` sh
### Script for Metacentrum

#PBS -N trim_galore_for_RNAseq
#PBS -l select=1:ncpus=4:mem=32gb:scratch_local=10gb:brno=True
#PBS -l walltime=24:00:00 
#PBS -m ae

## trim_galore

# define a DATADIR variable: directory where the input files are taken from
DATADIR=/storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/rnaseq/1_raw_reads
# directory for output
OUTDIR=/storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/rnaseq/2_trimmed_reads

# append a line to a file "jobs_info.txt" containing the ID of the job, the hostname of node it is run on and the path to a scratch directory
# this information helps to find a scratch directory in case the job fails and you need to remove the scratch directory manually 
echo "$PBS_JOBID is running on node `hostname -f` in a scratch directory $SCRATCHDIR" >> $PBS_O_WORKDIR/jobs_info.txt

# move into data directory
cd $DATADIR

# running trim_galore
module load trim_galore/0.6.2_py3
trim_galore --version # version 0.6.2

# There is needed "sort" because the manual says:
# "Trim Galore! expects paired-end files to be supplied in a pairwise fashion, e.g. file1_1.fq file1_2.fq SRR2_1.fq.gz SRR2_2.fq.gz ... ."
time find $DATADIR -type f \( -iname "*.fastq.gz" -o -iname "*.fq.gz" \) | sort | xargs trim_galore --paired --cores 4 -o $OUTDIR

# Remove trim_galore module
# The MultiQC otherwise does not work.
module remove trim_galore/0.6.2_py3

## FastQC and MultiQC after trim_galore

#### Files after trim_galore trimming end with *.fq.gz and not *.fastq.gz!

# report name
report_name="Cardamine_RNAseq_MultiQC_report_after_trimming"

# create the fastqc dir if it does not exists
if [ ! -d $OUTDIR/fastqc ]; then 
mkdir $OUTDIR/fastqc
fi

# create the multiqc dir if it does not exists
if [ ! -d $OUTDIR/multiqc ]; then
mkdir $OUTDIR/multiqc
fi

# test if scratch directory is set
# if scratch directory is not set, issue error message and exit
#test -n "$SCRATCHDIR" || { echo >&2 "Variable SCRATCHDIR is not set!"; exit 1; }

# load Java
module add openjdk/

# FastQC run

# Version 1 for files in multiple folders
# time find $DATADIR -type f \( -iname "*.fastq.gz" -o -iname "*.fq.gz" \) -print0 | xargs -0 /storage/brno12-cerit/home/duchmil/SW/fastqc/FastQC/fastqc -t 10 -o $OUTDIR/fastqc -f fastq

# Version 2 for files in single folder (not going to subfolders)
time /storage/brno12-cerit/home/duchmil/SW/fastqc/FastQC/fastqc -t 4 -o $OUTDIR/fastqc -f fastq $OUTDIR/*.fq.gz


# Running MultiQC
# activation
source /storage/brno2/home/duchmil/SW/mambaforge/bin/activate Multiqc

# run MultiQC
time multiqc --filename $report_name --outdir $OUTDIR/multiqc $OUTDIR/fastqc

# clean the SCRATCH directory
clean_scratch

# Resources: 1 h, 79 % CPU, 26 GB memory.
```

Checking trimming reports

``` sh
cd /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/rnaseq/2_trimmed_reads
grep 'Reads with adapters:' *trimming_report.txt # 48-49 % for Macrogen, 35 % for Novogene
grep 'Total written (filtered):' *trimming_report.txt # 96 % for Macrogen, over 99 % for Novogene
```

MultiQC report: Adapters from Macrogen samples were succesfully removed.

## Alignment of RNAseq reads with HISAT2

Braker instructions
(<https://github.com/Gaius-Augustus/BRAKER#braker-with-rna-seq-and-protein-data>):

“GeneMark-ETP utilizes Stringtie2 to assemble RNA-Seq data, which
requires that the aligned reads (BAM files) contain the XS (strand) tag
for spliced reads. Therefore, if you align your reads with HISAT2, you
must enable the –dta option.”

HISAT2 manual

<http://daehwankimlab.github.io/hisat2/manual/>

Folder for results

``` sh
cd /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/rnaseq
mkdir 3_aligned_reads
```

### Running mapping with Hisat2

Conversion to bam, sorting, indexing and merging of the bam and indexing
of the merged bam using samtools is included in the script.

Test

``` sh
cd /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/rnaseq/2_trimmed_reads
find . -type f -name "*_2_val_2.fq.gz" | sort | sed 's,./,,' | sed 's,_2_val_2.fq.gz,,' | while read my_sample
do
echo $my_sample
done
```

``` sh
### Script for Metacentrum

#PBS -N hisat2_alignment_RNAseq
#PBS -l select=1:ncpus=16:mem=128gb:scratch_local=1000gb
#PBS -l walltime=24:00:00 
#PBS -m ae


# define a DATADIR variable: directory where the input files are taken from and where output will be copied to
DATADIR=/storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06
# assembly name
genome_assembly="Cardamine_glauca_CUNI_V1_2024_10_masked.fa"

# append a line to a file "jobs_info.txt" containing the ID of the job, the hostname of node it is run on and the path to a scratch directory
# this information helps to find a scratch directory in case the job fails and you need to remove the scratch directory manually 
echo "$PBS_JOBID is running on node `hostname -f` in a scratch directory $SCRATCHDIR" | ts '[%Y-%m-%d %H:%M:%S]' >> $PBS_O_WORKDIR/jobs_info.txt

# loading modules needed
module load hisat2
hisat2 --version
module load samtools
samtools --version

echo "Modules loaded" | ts '[%Y-%m-%d %H:%M:%S]'


# test if scratch directory is set
# if scratch directory is not set, issue error message and exit
test -n "$SCRATCHDIR" || { echo >&2 "Variable SCRATCHDIR is not set!"; exit 1; }

# copy reference genome
cp -v -r $DATADIR/genome_assembly/$genome_assembly $SCRATCHDIR || { echo >&2 "Error while copying index file(s)!"; exit 2; }

# copy RNAseq files
cp -v -r $DATADIR/rnaseq/2_trimmed_reads/*.fq.gz $SCRATCHDIR || { echo >&2 "Error while copying fastq file(s)!"; exit 2; }

echo "Input files copied." | ts '[%Y-%m-%d %H:%M:%S]'

# move into scratch directory
cd $SCRATCHDIR

# Build the index of reference sequence
hisat2-build -f -p 16 $genome_assembly hisat2_index

echo "Hisat2 index done." | ts '[%Y-%m-%d %H:%M:%S]'

# running HISAT2
find . -type f -name "*_2_val_2.fq.gz" | sort | sed 's,./,,' | sed 's,_2_val_2.fq.gz,,' | while read my_sample
do
  hisat2 -q -t --dta -p 16 -x hisat2_index \
   -1  $my_sample"_1_val_1.fq.gz"\
   -2  $my_sample"_2_val_2.fq.gz"\
   -S $my_sample"_trimmed.sam" | ts '[%Y-%m-%d %H:%M:%S]'
   
  echo "Alignment for $my_sample done." | ts '[%Y-%m-%d %H:%M:%S]'
  
  # conversion to bam 
  samtools view -bS --threads 16 $my_sample"_trimmed.sam" > $my_sample"_trimmed.bam"
  # sorting of bam
  samtools sort --threads 16 $my_sample"_trimmed.bam" -o $my_sample"_trimmed_sorted.bam"
  # indexing of bam
  samtools index -@ 16 $my_sample"_trimmed_sorted.bam"
  
  echo "Samtools conversion to bam, sorting and indexing for $my_sample done." | ts '[%Y-%m-%d %H:%M:%S]'
  
  # move the output to user's DATADIR or exit in case of failure
  cp -v $my_sample"_trimmed_sorted.bam" $DATADIR/rnaseq/3_aligned_reads/ || { echo >&2 "Result file(s) copying failed (with a code $?) !!"; exit 4; }
  cp -v $my_sample"_trimmed_sorted.bam.bai" $DATADIR/rnaseq/3_aligned_reads/ || { echo >&2 "Result file(s) copying failed (with a code $?) !!"; exit 4; }
  
  echo "Copying output files for $my_sample done." | ts '[%Y-%m-%d %H:%M:%S]'

done

# Merging of bams
samtools merge -@ 16 -o RNAseq_trimmed_merged.bam *_trimmed_sorted.bam
# indexing of merged bam
samtools index -@ 16 RNAseq_trimmed_merged.bam
cp -v RNAseq_trimmed_merged.bam* $DATADIR/rnaseq/3_aligned_reads/ || { echo >&2 "Result file(s) copying failed (with a code $?) !!"; exit 4; }

echo "Merging of bams done." | ts '[%Y-%m-%d %H:%M:%S]'

# clean the SCRATCH directory
clean_scratch

# Resources: 1 h 15 min, 43 % CPU, 100 % memory
# It seems it was somehow stucked this time, according to log the copying of files was very slow. The previous run was:
# Resources: 40 min, 70 % CPU, 100 % memory
```

### Checking logs

Output statistics: hisat2_alignment_RNAseq.e11116576

``` sh
cd /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/metacentrum_scripts
# Most important statistics
o_file=hisat2_alignment_RNAseq.o11116576
e_file=hisat2_alignment_RNAseq.e11116576
paste <(grep 'Alignment for ' $o_file | sed 's/^.*Alignment for //' | sed 's/ done.//') <(grep 'overall alignment rate' $e_file)  <(grep 'aligned concordantly exactly 1 time' $e_file) <(grep 'aligned concordantly >1 times' $e_file) | column -t | less -S
```

    Cardamine

    XC174YL_Macrogen                        88.61%  overall  alignment  rate  15073594  (80.53%)  aligned  concordantly  exactly  1  time  821422   (4.39%)  aligned  concordantly  >1  times
    XC174_OL                                89.65%  overall  alignment  rate  16788949  (78.12%)  aligned  concordantly  exactly  1  time  644914   (3.00%)  aligned  concordantly  >1  times
    XC174_R                                 88.84%  overall  alignment  rate  13775354  (79.71%)  aligned  concordantly  exactly  1  time  460806   (2.67%)  aligned  concordantly  >1  times
    XC174_S                                 88.40%  overall  alignment  rate  13977235  (76.99%)  aligned  concordantly  exactly  1  time  543087   (2.99%)  aligned  concordantly  >1  times
    XC174_YL                                89.91%  overall  alignment  rate  14933580  (77.70%)  aligned  concordantly  exactly  1  time  581419   (3.03%)  aligned  concordantly  >1  times
    XC618_FB_MKRN250017681-1A_22VTVNLT4_L6  87.68%  overall  alignment  rate  22868036  (77.23%)  aligned  concordantly  exactly  1  time  827162   (2.79%)  aligned  concordantly  >1  times
    XC618_F_MKRN250017680-1A_22VTVNLT4_L4   85.99%  overall  alignment  rate  29896456  (75.43%)  aligned  concordantly  exactly  1  time  1160215  (2.93%)  aligned  concordantly  >1  times


    Noccaea

    SRR24947461  62.94%  overall  alignment  rate  52688337  (52.47%)  aligned  concordantly  exactly  1  time  3947151   (3.93%)   aligned  concordantly  >1  times
    SRR24947462  62.06%  overall  alignment  rate  65494606  (52.63%)  aligned  concordantly  exactly  1  time  3078652   (2.47%)   aligned  concordantly  >1  times
    SRR24947463  57.88%  overall  alignment  rate  61563211  (50.41%)  aligned  concordantly  exactly  1  time  2530332   (2.07%)   aligned  concordantly  >1  times
    SRR24947464  74.88%  overall  alignment  rate  75143196  (65.64%)  aligned  concordantly  exactly  1  time  3670500   (3.21%)   aligned  concordantly  >1  times
    SRR24947465  78.91%  overall  alignment  rate  73881131  (69.84%)  aligned  concordantly  exactly  1  time  3821305   (3.61%)   aligned  concordantly  >1  times
    SRR24947466  70.50%  overall  alignment  rate  59916569  (49.11%)  aligned  concordantly  exactly  1  time  20972725  (17.19%)  aligned  concordantly  >1  times
    XN580_FB     80.79%  overall  alignment  rate  13836893  (68.56%)  aligned  concordantly  exactly  1  time  470979    (2.33%)   aligned  concordantly  >1  times
    XN580_F      80.49%  overall  alignment  rate  11639449  (68.42%)  aligned  concordantly  exactly  1  time  416401    (2.45%)   aligned  concordantly  >1  times
    XN580_L      81.19%  overall  alignment  rate  12380457  (68.45%)  aligned  concordantly  exactly  1  time  547587    (3.03%)   aligned  concordantly  >1  times
    XN580_R      76.58%  overall  alignment  rate  11855116  (64.63%)  aligned  concordantly  exactly  1  time  390891    (2.13%)   aligned  concordantly  >1  times
    XN580_S      80.39%  overall  alignment  rate  11335150  (68.29%)  aligned  concordantly  exactly  1  time  370158    (2.23%)   aligned  concordantly  >1  times



    Alyssum
    AM08M_CF        82.26% overall alignment rate       15395977 (72.53%) aligned concordantly exactly 1 time           1103653 (5.20%) aligned concordantly >1 times
    AM08M_RO        80.78% overall alignment rate       14909449 (72.43%) aligned concordantly exactly 1 time           710062 (3.45%) aligned concordantly >1 times
    AM08N_CF        75.13% overall alignment rate       13619900 (66.20%) aligned concordantly exactly 1 time           987032 (4.80%) aligned concordantly >1 times
    AM08N_OF        68.77% overall alignment rate       12482837 (59.82%) aligned concordantly exactly 1 time           1164331 (5.58%) aligned concordantly >1 times
    AM08N_OL        71.28% overall alignment rate       12701049 (59.93%) aligned concordantly exactly 1 time           1622373 (7.66%) aligned concordantly >1 times
    AM08N_RO        75.34% overall alignment rate       13749781 (67.05%) aligned concordantly exactly 1 time           790056 (3.85%) aligned concordantly >1 times
    AM08N_YL        69.47% overall alignment rate       12755339 (60.27%) aligned concordantly exactly 1 time           1153072 (5.45%) aligned concordantly >1 times

## Preparation of protein sequences

Protein sequences downloaded from
<https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/>
(<https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/Viridiplantae.fa.gz>).

Protein sequences should be decompressed (it should be normal fasta).

``` sh
# OrthoDB Viridiplantae proteins are already downloaded (they are shared for all species)

# cd /storage/brno12-cerit/home/duchmil/annotations/OrthoDB_proteins/
# wget https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/Viridiplantae.fa.gz
# gunzip Viridiplantae.fa.gz
```

# Running Braker

``` sh
### Script for Metacentrum

#PBS -N braker_Cardamine_01
#PBS -l select=1:ncpus=12:mem=96gb:scratch_local=1000gb
#PBS -l walltime=24:00:00 
#PBS -m ae

# define a DATADIR variable: directory where the input files are taken from and where output will be copied to
DATADIR=/storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06
# Name of genome assembly
genome_assembly="Cardamine_glauca_CUNI_V1_2024_10_masked.fa"

# append a line to a file "jobs_info.txt" containing the ID of the job, the hostname of node it is run on and the path to a scratch directory
# this information helps to find a scratch directory in case the job fails and you need to remove the scratch directory manually 
echo "$PBS_JOBID is running on node `hostname -f` in a scratch directory $SCRATCHDIR" | ts '[%Y-%m-%d %H:%M:%S]' >> $PBS_O_WORKDIR/jobs_info.txt

# test if scratch directory is set
# if scratch directory is not set, issue error message and exit
test -n "$SCRATCHDIR" || { echo >&2 "Variable SCRATCHDIR is not set!"; exit 1; }

# copy files
cp -r $DATADIR/rnaseq/3_aligned_reads/RNAseq_trimmed_merged.bam $SCRATCHDIR || { echo >&2 "Error while copying input file(s)!"; exit 2; }
cp -r /storage/brno12-cerit/home/duchmil/annotations/OrthoDB_proteins/Viridiplantae.fa $SCRATCHDIR || { echo >&2 "Error while copying input file(s)!"; exit 2; }
cp -r $DATADIR/genome_assembly/$genome_assembly $SCRATCHDIR || { echo >&2 "Error while copying input file(s)!"; exit 2; }


# move into scratch directory
cd $SCRATCHDIR 
mkdir results_braker_01

# running BRAKER
export BRAKER_SIF=/storage/brno12-cerit/home/duchmil/SW/braker_sw/braker3.sif

singularity exec -B ${PWD}:${PWD} ${BRAKER_SIF} braker.pl --bam=RNAseq_trimmed_merged.bam --genome=$genome_assembly --prot_seq=Viridiplantae.fa --threads=12 --species=Cardamine_glauca --workingdir=$SCRATCHDIR/results_braker_01
# Note: Protein sequences shouldn't be compressed. It should be plain fasta.


# move the output to user's DATADIR or exit in case of failure
cp -r results_braker_01 $DATADIR/ || { echo >&2 "Result file(s) copying failed (with a code $?) !!"; exit 4; }

# clean the SCRATCH directory
clean_scratch

# Resources: The job was running 8 h, using 93 GB memory and 30% of CPU time.
```

# Quality check of the Braker annotation

## Running IGV on Metacentrum

``` sh
# Download IGV
wget https://data.broadinstitute.org/igv/projects/downloads/2.16/IGV_2.16.2.zip
# Unzip IGV
unzip IGV_2.16.2.zip
# Start interactive job
qsub -I -l select=1:ncpus=2:mem=32gb:scratch_local=40gb -l walltime=5:00:00
# load and start GUI
module add gui
gui start
# Open the URL in browser.
# Open terminal in GUI.
# Run these commands:
# load Java
module load openjdk
# run IGV
/storage/brno12-cerit/home/duchmil/SW/igv/IGV_2.16.2/igv.sh

# run with direct opening of the right files
/storage/brno12-cerit/home/duchmil/SW/igv/IGV_2.16.2/igv.sh -g /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/genome_assembly/Cardamine_glauca_CUNI_V1_2024_10_masked.fa /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/results_braker_01/braker.gtf /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/rnaseq/3_aligned_reads/*_sorted.bam /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/rnaseq/4_assembled_transcripts/Assembled_transcripts.gtf /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/annot_processing/repeat_annotation/uppercase_cardamine_polished_filppedSc_4.fa.out.gff /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2024_12/annot_processing/Cardamine_glauca_v1_annotation.gff
```

## Number of genes

``` sh
cd /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/results_braker_01/
grep -c -P "\tgene\t" braker.gtf # 24828 (24376 previous run)
```

## Checking proteins predicted by Braker for stop codons

``` sh
# converting from folded fasta to unfolded fasta for better counting and counting internal stop codons
awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' < braker.aa | grep \*[[:alpha:]] | wc -l # 0

awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' < braker.aa | grep -B1 '\*[[:alpha:]]'
```

## Alignment of protein sequences to the genome assembly

Miniprot

<https://github.com/lh3/miniprot>

The alignment of the proteins will be needed to check the annotation by
Braker.

Generally proteins of closely related species with good annotation
should be used. We used protein sequences of B. rapa, A. thaliana and A.
lyrata.

Data:

B. rapa web page: <http://brassicadb.cn>

<http://39.100.233.196:82/download_genome/Brassica_Genome_data/Brara_Chiifu_V3.0/Brapa_genome_v3.0_pep.fasta.gz>

<http://39.100.233.196:82/download_genome/Brassica_Genome_data/Brara_Chiifu_V4.1/Brapa_chiifu_v41_gene20230413.gff3.pep.fa.gz>

A. thaliana

<https://www.arabidopsis.org/download/file?path=Proteins/Araport11_protein_lists/Araport11_pep_20220914.gz>

A. lyrata NCBI 101 annotation

<https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/255/GCF_000004255.2_v.1.0/GCF_000004255.2_v.1.0_protein.faa.gz>

### Getting protein sequences

``` sh
cd /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/
mkdir protein_seqs_input
cd protein_seqs_input
mkdir 1_downloaded
cd 1_downloaded



# download protein sequences
cd /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/protein_seqs_input/1_downloaded
wget http://39.100.233.196:82/download_genome/Brassica_Genome_data/Brara_Chiifu_V4.1/Brapa_chiifu_v41_gene20230413.gff3.pep.fa.gz

# this command does not work, I had to download it to my computer and then copy to Metacentrum
# wget "https://www.arabidopsis.org/download/file?path=Proteins/Araport11_protein_lists/Araport11_pep_20220914.gz"
# copy from Alyssum annotation
cp -v /storage/brno12-cerit/home/duchmil/annotations/alyssum_2024_Mahnaz_assembly/protein_seqs_input/1_downloaded/Araport11_pep_20220914.gz .

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/255/GCF_000004255.2_v.1.0/GCF_000004255.2_v.1.0_protein.faa.gz


# checking number of sequences
zgrep -c '>' Araport11_pep_20220914.gz # 48266
zgrep -c '>' Brapa_chiifu_v41_gene20230413.gff3.pep.fa.gz # 83470
zgrep -c '>' GCF_000004255.2_v.1.0_protein.faa.gz # 39161


# check stop codons
# converting from folded fasta to unfolded fasta for better counting and counting internal stop codons
zcat GCF_000004255.2_v.1.0_protein.faa.gz | awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' | grep \*[[:alpha:]] | wc -l # 0

# make folder for results
cd /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/protein_seqs_input
mkdir 2_aligned
```

### Miniprot alignment

``` sh
### Script for Metacentrum

#PBS -N miniprot_alignment_01
#PBS -l select=1:ncpus=4:mem=8gb:scratch_local=1000gb
#PBS -l walltime=2:00:00 
#PBS -m ae

# define a DATADIR variable: directory where the input files are taken from and where output will be copied to
DATADIR=/storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06
# Name of genome assembly
genome_assembly="Cardamine_glauca_CUNI_V1_2024_10_masked.fa"

MINIPROTDIR=/storage/brno12-cerit/home/duchmil/SW/miniprot

# append a line to a file "jobs_info.txt" containing the ID of the job, the hostname of node it is run on and the path to a scratch directory
# this information helps to find a scratch directory in case the job fails and you need to remove the scratch directory manually 
echo "$PBS_JOBID is running on node `hostname -f` in a scratch directory $SCRATCHDIR" | ts '[%Y-%m-%d %H:%M:%S]' >> $PBS_O_WORKDIR/jobs_info.txt

# test if scratch directory is set
# if scratch directory is not set, issue error message and exit
test -n "$SCRATCHDIR" || { echo >&2 "Variable SCRATCHDIR is not set!"; exit 1; }

# copy files
cp -r $DATADIR/protein_seqs_input/1_downloaded/Araport11_pep_20220914.gz \
$DATADIR/protein_seqs_input/1_downloaded/Brapa_chiifu_v41_gene20230413.gff3.pep.fa.gz \
$DATADIR/protein_seqs_input/1_downloaded/GCF_000004255.2_v.1.0_protein.faa.gz \
$SCRATCHDIR || { echo >&2 "Error while copying input file(s)!"; exit 2; }

cp -r $DATADIR/genome_assembly/$genome_assembly $SCRATCHDIR || { echo >&2 "Error while copying input file(s)!"; exit 2; }

echo "Input files copied." | ts '[%Y-%m-%d %H:%M:%S]'


# move into scratch directory
cd $SCRATCHDIR 

# make index for miniprot
$MINIPROTDIR/miniprot -t4 -d miniprot_index.mpi $genome_assembly

echo "Miniprot index prepared." | ts '[%Y-%m-%d %H:%M:%S]'

# running miniprot
$MINIPROTDIR/miniprot -Iut4 --gff miniprot_index.mpi Araport11_pep_20220914.gz > A.thaliana.pep.gff
echo "A. thaliana done." | ts '[%Y-%m-%d %H:%M:%S]'
$MINIPROTDIR/miniprot -Iut4 --gff miniprot_index.mpi Brapa_chiifu_v41_gene20230413.gff3.pep.fa.gz > B.rapa.pep.gff
echo "B. rapa done." | ts '[%Y-%m-%d %H:%M:%S]'
$MINIPROTDIR/miniprot -Iut4 --gff miniprot_index.mpi GCF_000004255.2_v.1.0_protein.faa.gz > A.lyrata.pep.gff
echo "A. lyrata done." | ts '[%Y-%m-%d %H:%M:%S]'




# move the output to user's DATADIR or exit in case of failure
cp -v *.gff $DATADIR/protein_seqs_input/2_aligned | ts '[%Y-%m-%d %H:%M:%S]' || { echo >&2 "Result file(s) copying failed (with a code $?) !!"; exit 4; }

# clean the SCRATCH directory
clean_scratch

# Resources: 11 min, 97 % CPU, only 3 GB memory
```

## Assembly of transcripts from RNAseq alignment

The assembled transcipts will be used to:

- check Braker annotation visually using IGV,
- extract sequences from them using AGAT and run Busco on them,
- check the overlaps of them with annotated genes.

stringtie to produce transcripts from aligned RNAseq reads

- default parameters should be fine
- set the number of threads

``` sh
# folder for results
cd /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/rnaseq
mkdir 4_assembled_transcripts
```

``` sh
### Script for Metacentrum

#PBS -N stringtie_transcript_assembly
#PBS -l select=1:ncpus=1:mem=40gb:scratch_local=1000gb
#PBS -l walltime=2:00:00 
#PBS -m ae

# define a DATADIR variable: directory where the input files are taken from and where output will be copied to
DATADIR=/storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/rnaseq

STRINGTIEDIR=/storage/brno12-cerit/home/duchmil/SW/stringtie

# append a line to a file "jobs_info.txt" containing the ID of the job, the hostname of node it is run on and the path to a scratch directory
# this information helps to find a scratch directory in case the job fails and you need to remove the scratch directory manually 
echo "$PBS_JOBID is running on node `hostname -f` in a scratch directory $SCRATCHDIR" | ts '[%Y-%m-%d %H:%M:%S]' >> $PBS_O_WORKDIR/jobs_info.txt

# test if scratch directory is set
# if scratch directory is not set, issue error message and exit
test -n "$SCRATCHDIR" || { echo >&2 "Variable SCRATCHDIR is not set!"; exit 1; }

# copy files
cp -rv $DATADIR/3_aligned_reads/RNAseq_trimmed_merged.bam* $SCRATCHDIR || { echo >&2 "Error while copying index file(s)!"; exit 2; }

echo "Input files copied." | ts '[%Y-%m-%d %H:%M:%S]'


# move into scratch directory
cd $SCRATCHDIR 


# running Stringtie
$STRINGTIEDIR/stringtie -p 1 -o Assembled_transcripts.gtf RNAseq_trimmed_merged.bam || { echo >&2 "Running of main command failed (with a code $?) !!"; exit 3; }

echo "Stringtie assembly done." | ts '[%Y-%m-%d %H:%M:%S]'

# move the output to user's DATADIR or exit in case of failure
mkdir -p $DATADIR/4_assembled_transcripts
cp -v Assembled_transcripts.gtf $DATADIR/4_assembled_transcripts || { echo >&2 "Result file(s) copying failed (with a code $?) !!"; exit 4; }

echo "Output files copied." | ts '[%Y-%m-%d %H:%M:%S]'

# clean the SCRATCH directory
clean_scratch


# Note to computational resources:
# The script ran 13 min, with 86 % CPU time and 26 GB memory used.
```

## Extraction of transcript sequences

Extracted transcript sequences will be used to run Busco on them.

### Extraction of transcript (mRNA) sequences

``` sh
# interactive job
qsub -I -l select=1:ncpus=1:mem=8gb:scratch_local=10gb -l walltime=2:00:00

# Name of genome assembly
genome_assembly="Cardamine_glauca_CUNI_V1_2024_10_masked.fa"

cd /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/rnaseq/4_assembled_transcripts

# run the container
singularity run /storage/brno12-cerit/home/duchmil/SW/agat/agat_1.4.0--pl5321hdfd78af_0.sif

# extraction of mRNA (UTRs + CDS)
agat_sp_extract_sequences.pl -g Assembled_transcripts.gtf -f /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/genome_assembly/$genome_assembly -t exon --merge -o Assembled_transcripts.fasta | tee log_transcripts_AGAT.txt

exit

# Number of transcripts
grep -P -c '\ttranscript\t' Assembled_transcripts.gtf # 45030
grep -c '>' Assembled_transcripts.fasta # 45030


# Counting of genes
grep -P '\ttranscript\t' Assembled_transcripts.gtf | cut -f 9 | sed 's/; transcript_id.*$//' | head -n 50
grep -P '\ttranscript\t' Assembled_transcripts.gtf | cut -f 9 | sed 's/; transcript_id.*$//' | sort | uniq | wc -l # 31004 (previous run: 27054)
grep -P '\ttranscript\t' Assembled_transcripts.gtf | cut -f 9 | sed 's/; transcript_id.*$//' | tail # the highest number is 31004 (previous run: 27054)
```

For this run, we have 45030 transcript. The number of transcripts in
previous run (without flowers and flower buds RNAseq) was 39064, so the
increase is significant. The flowers and flower buds were sequenced
deeper compared to other samples (60M reads vs. 30M reads).

Similarly, the number of genes from StringTie (alternative transcripts
in one locus are assigned to one gene) increased from 27054 to 31004.

## Busco

``` sh
cd /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/
mkdir busco_results
cd busco_results
```

### Running Busco

``` sh
### Script for Metacentrum

#PBS -N Busco_proteins_CDS_assembly_transcripts
#PBS -l select=1:ncpus=4:mem=8gb:scratch_local=1000gb
#PBS -l walltime=4:00:00 
#PBS -m ae

# define a DATADIR variable: directory where the input files are taken from and where output will be copied to
DATADIR=/storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06
# Name of genome assembly
genome_assembly="Cardamine_glauca_CUNI_V1_2024_10_masked.fa"

# append a line to a file "jobs_info.txt" containing the ID of the job, the hostname of node it is run on and the path to a scratch directory
# this information helps to find a scratch directory in case the job fails and you need to remove the scratch directory manually 
echo "$PBS_JOBID is running on node `hostname -f` in a scratch directory $SCRATCHDIR" | ts '[%Y-%m-%d %H:%M:%S]' >> $PBS_O_WORKDIR/jobs_info.txt

# test if scratch directory is set
# if scratch directory is not set, issue error message and exit
test -n "$SCRATCHDIR" || { echo >&2 "Variable SCRATCHDIR is not set!"; exit 1; }

# copy files
cp -v $DATADIR/results_braker_01/braker.aa $SCRATCHDIR || { echo >&2 "Error while copying index file(s)!"; exit 2; }
cp -v $DATADIR/results_braker_01/braker.codingseq $SCRATCHDIR || { echo >&2 "Error while copying index file(s)!"; exit 2; }
cp -v $DATADIR/genome_assembly/$genome_assembly $SCRATCHDIR || { echo >&2 "Error while copying index file(s)!"; exit 2; }
cp -v $DATADIR/rnaseq/4_assembled_transcripts/Assembled_transcripts.fasta $SCRATCHDIR || { echo >&2 "Error while copying index file(s)!"; exit 2; }



echo "Input files copied." | ts '[%Y-%m-%d %H:%M:%S]'


# move into scratch directory
cd $SCRATCHDIR 


# activate Busco
source /storage/brno2/home/duchmil/SW/mambaforge/bin/activate busco_5_7_1

# run Busco
busco --in braker.aa --mode proteins --lineage_dataset brassicales_odb10 --cpu 4
echo "Busco for proteins done." | ts '[%Y-%m-%d %H:%M:%S]'

busco --in braker.codingseq --mode transcriptome --lineage_dataset brassicales_odb10 --cpu 4
echo "Busco for CDS done." | ts '[%Y-%m-%d %H:%M:%S]'

busco --in $genome_assembly --mode genome --lineage_dataset brassicales_odb10 --cpu 4
echo "Busco for assembly done." | ts '[%Y-%m-%d %H:%M:%S]'

busco --in Assembled_transcripts.fasta --mode transcriptome --lineage_dataset brassicales_odb10 --cpu 4
echo "Busco for transcripts done." | ts '[%Y-%m-%d %H:%M:%S]'


echo "Main computation done." | ts '[%Y-%m-%d %H:%M:%S]'

# make folder for results (if not existing)
mkdir -p $DATADIR/busco_results
# move the output to user's DATADIR or exit in case of failure
cp -vR BUSCO_* $DATADIR/busco_results || { echo >&2 "Result file(s) copying failed (with a code $?) !!"; exit 4; }

echo "Output files copied." | ts '[%Y-%m-%d %H:%M:%S]'

# clean the SCRATCH directory
clean_scratch

# Resources: The script ran 2,5 h, with 73 % CPU time and 100% memory used.
```

### Plot Busco results

``` sh
cd /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/busco_results/

mkdir summaries_BUSCO

# copy summaries from all Busco results folders
find BUSCO* -name "short_summary.specific.*.txt" -exec cp {} summaries_BUSCO/ \;

# Rename some summaries (avoid dot in specific part of the name) so that Busco script will distinguish them
cd summaries_BUSCO/
mv short_summary.specific.brassicales_odb10.BUSCO_braker.aa.txt short_summary.specific.brassicales_odb10.BUSCO_braker_aa.txt
mv short_summary.specific.brassicales_odb10.BUSCO_braker.codingseq.txt short_summary.specific.brassicales_odb10.BUSCO_braker_codingseq.txt


# activate Busco
source /storage/brno2/home/duchmil/SW/mambaforge/bin/activate busco_5_7_1

# Use script to generate plot
generate_plot.py -wd .
```

### Busco results

#### Assembly

    C:99.0%[S:94.2%,D:4.8%],F:0.2%,M:0.8%,n:4596,E:1.8%    
    4549    Complete BUSCOs (C) (of which 82 contain internal stop codons)         
    4329    Complete and single-copy BUSCOs (S)    
    220 Complete and duplicated BUSCOs (D)     
    9   Fragmented BUSCOs (F)              
    38  Missing BUSCOs (M)             
    4596    Total BUSCO groups searched        

#### Transcripts (RNAseq reads aligned to genome and assembled to transcripts)

    C:91.9%[S:58.3%,D:33.6%],F:1.3%,M:6.8%,n:4596      
    4222    Complete BUSCOs (C)            
    2679    Complete and single-copy BUSCOs (S)    
    1543    Complete and duplicated BUSCOs (D)     
    61  Fragmented BUSCOs (F)              
    313 Missing BUSCOs (M)             
    4596    Total BUSCO groups searched 

#### Annotation (proteins predicted by Braker)

    C:98.1%[S:83.9%,D:14.2%],F:0.1%,M:1.8%,n:4596      
    4506    Complete BUSCOs (C)            
    3854    Complete and single-copy BUSCOs (S)    
    652 Complete and duplicated BUSCOs (D)     
    6   Fragmented BUSCOs (F)              
    84  Missing BUSCOs (M)             
    4596    Total BUSCO groups searched   

#### Annotation (CDS predicted by Braker)

    C:98.0%[S:83.9%,D:14.1%],F:0.2%,M:1.8%,n:4596      
    4505    Complete BUSCOs (C)            
    3857    Complete and single-copy BUSCOs (S)    
    648 Complete and duplicated BUSCOs (D)     
    7   Fragmented BUSCOs (F)              
    84  Missing BUSCOs (M)             
    4596    Total BUSCO groups searched            

### Busco results from previous run without RNA from flowers and flower buds

#### Transcripts (RNAseq reads aligned to genome and assembled to transcripts)

    C:87.4%[S:56.7%,D:30.7%],F:2.2%,M:10.4%,n:4596     
    4019    Complete BUSCOs (C)            
    2607    Complete and single-copy BUSCOs (S)    
    1412    Complete and duplicated BUSCOs (D)     
    100 Fragmented BUSCOs (F)              
    477 Missing BUSCOs (M)             
    4596    Total BUSCO groups searched            

#### Annotation (proteins predicted by Braker)

    C:97.5%[S:82.9%,D:14.6%],F:0.2%,M:2.3%,n:4596      
    4482    Complete BUSCOs (C)            
    3811    Complete and single-copy BUSCOs (S)    
    671 Complete and duplicated BUSCOs (D)     
    7   Fragmented BUSCOs (F)              
    107 Missing BUSCOs (M)             
    4596    Total BUSCO groups searched    

#### Annotation (CDS predicted by Braker)

    C:97.5%[S:83.0%,D:14.5%],F:0.2%,M:2.3%,n:4596      
    4481    Complete BUSCOs (C)            
    3813    Complete and single-copy BUSCOs (S)    
    668 Complete and duplicated BUSCOs (D)     
    8   Fragmented BUSCOs (F)              
    107 Missing BUSCOs (M)             
    4596    Total BUSCO groups searched        

#### Comments to results

- All Busco results are very similar to Noccaea, even that Cardamine has
  much less annotated genes.
- Including RNAseq from flowers and flower buds increased the Busco
  scores from 97.5% to 98.1%.

## Intersect between gene models and aligned proteins

bedtools intersect

<https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html>

``` sh
cd /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06
mkdir intersects
cd intersects

module load bedtools2/2.30.0-gcc-10.2.1-5acjqve

ls ../protein_seqs_input/2_aligned/

# Intersects with proteins from single species
bedtools intersect -a ../results_braker_01/braker.gtf -b ../protein_seqs_input/2_aligned/A.thaliana.pep.gff -u -f 0.9 -r > braker_annotation_x_A.thaliana.pep.gff_intersect.tab

bedtools intersect -a ../results_braker_01/braker.gtf -b ../protein_seqs_input/2_aligned/A.lyrata.pep.gff -u -f 0.9 -r > braker_annotation_x_A.lyrata.pep.gff_intersect.tab

bedtools intersect -a ../results_braker_01/braker.gtf -b ../protein_seqs_input/2_aligned/B.rapa.pep.gff -u -f 0.9 -r > braker_annotation_x_B.rapa.pep.gff_intersect.tab

# Intersect with all proteins together
bedtools intersect -a ../results_braker_01/braker.gtf -b ../protein_seqs_input/2_aligned/*  -u -f 0.9 -r > braker_annotation_x_all.pep.gff_intersect.tab
```

This command finds the intersects between Braker output (`-a`) and at
least one of the aligned proteins (files in `-b`). It will report the
genes in `-a` only once (option `-u`). The overlap should be at least 90
% (`-f 0.9`) for both `-a` and `-b` (option `-r` as reciprocal).

Checking the outputs

``` sh
grep -P -c "\tgene\t" ../results_braker_01/braker.gtf # 24828
grep -P -c "\tmRNA\t" ../protein_seqs_input/2_aligned/A.thaliana.pep.gff # 51011
grep -P -c "\tgene\t" braker_annotation_x_A.thaliana.pep.gff_intersect.tab # 19889

for FILE in braker_annotation_x_*
do
COUNT=$(grep -P -c "\tgene\t" $FILE)
echo "$FILE $COUNT"
done
```

**Cardamine:** braker_annotation_x_A.lyrata.pep.gff_intersect.tab 19974
braker_annotation_x_A.thaliana.pep.gff_intersect.tab 19889
braker_annotation_x_B.rapa.pep.gff_intersect.tab 18005
braker_annotation_x_all.pep.gff_intersect.tab 21160

**Cardamine (older run without RNA from flowers and flower buds):**
braker_annotation_x_A.lyrata.pep.gff_intersect.tab 19833
braker_annotation_x_A.thaliana.pep.gff_intersect.tab 19751
braker_annotation_x_B.rapa.pep.gff_intersect.tab 17925
braker_annotation_x_all.pep.gff_intersect.tab 20990

All of the numbers increased in the current run compared to previous
run.

**Noccaea:** braker_annotation_x_A.lyrata.pep.gff_intersect.tab 19105
braker_annotation_x_A.thaliana.pep.gff_intersect.tab 18949
braker_annotation_x_B.rapa.pep.gff_intersect.tab 17914
braker_annotation_x_all.pep.gff_intersect.tab 20757

**Alyssum:** braker_annotation_x_A.alpina.pep.gff_intersect.tab 16637
braker_annotation_x_A.lyrata.pep.gff_intersect.tab 18772
braker_annotation_x_A.saxatilis.pep.gff_intersect.tab 17095
braker_annotation_x_A.thaliana.pep.gff_intersect.tab 18676
braker_annotation_x_B.rapa.pep.gff_intersect.tab 17642
braker_annotation_x_all.pep.gff_intersect.tab 21660

Comment: The intersects for Noccaea are higher than for Alyssum.

## Intersect between gene models and assembled transcripts

Again bedtools intersect. This time we do not use the `-r` option - the
overlap should be at least 90 % of Braker output, but not necesarilly of
the Stringtie output. It is because the Stringtie tries to predict whole
transcripts including UTRs.

``` sh
cd /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/intersects

module load bedtools2/2.30.0-gcc-10.2.1-5acjqve

ls ../rnaseq/4_assembled_transcripts/

# Intersects with Stringtie assembled transcripts
bedtools intersect -a ../results_braker_01/braker.gtf -b ../rnaseq/4_assembled_transcripts/Assembled_transcripts.gtf -u -f 0.9 > braker_annotation_x_Assembled_transcripts.gtf_intersect.tab

# Number of genes predicted by Braker
grep -P -c "\tgene\t" ../results_braker_01/braker.gtf # 24828
# Number of transcripts assembled from RNAseq
grep -P -c "\ttranscript\t" ../rnaseq/4_assembled_transcripts/Assembled_transcripts.gtf # 45030
# Number of genes assembled from RNAseq
# (Gene features are not annotated separately by StringTie, but there is a field gene_id for transcripts which can be used.)
grep -P "\ttranscript\t" ../rnaseq/4_assembled_transcripts/Assembled_transcripts.gtf | sed -E 's/.*\tgene_id//' | sed -E 's/; transcript_id.*//' | sort -u | wc -l # 31004
# Intersect of genes predicted by Braker and transcripts assembled from RNAseq
grep -P -c "\tgene\t" braker_annotation_x_Assembled_transcripts.gtf_intersect.tab # 20404 (previous run: 18050)


# intersect between two intersects: (1) intersect between annotation and protein alignments and (2) intersect between annotation and transcripts
bedtools intersect -a braker_annotation_x_Assembled_transcripts.gtf_intersect.tab -b braker_annotation_x_all.pep.gff_intersect.tab -u -f 1.0 -r > braker_transcript_and_braker_protein_meta_intersect.tab
grep -P -c "\tgene\t" braker_transcript_and_braker_protein_meta_intersect.tab # 18298 (previous run: 16457)
# There should be '-f 1.0', because we are comparing genes, mRNAs and exons from the same Braker prediction. If there is 0.9, the number of genes is slightly higher than it should be, because it probably takes also some intersects of genes with transcripts or something like that.

# Alternative way to calculate number of genes in intersect of intersects (metaintersect):
comm -12 <(grep -P "\tgene\t" braker_annotation_x_all.pep.gff_intersect.tab | cut -f 9 | sort -u) <(grep -P "\tgene\t" braker_annotation_x_Assembled_transcripts.gtf_intersect.tab | cut -f 9 | sort -u) | wc -l # 18298 (previous run: 16457)
```

We will rather take union of the two intersects than the intersect of
them as the “confident” set of gene models. These will be gene models
that are either supported by RNAseq reads (more precisely assembled
transcripts) or by proteins from related species.

## Getting the list of “reliable” genes

Genes that are either supported by RNAseq reads (more precisely
assembled transcripts) or by proteins from related species.

``` sh
# checks
grep -P "\tgene\t" braker_annotation_x_all.pep.gff_intersect.tab | cut -f 9 | sort | head -n 20
grep -P "\tgene\t" braker_annotation_x_Assembled_transcripts.gtf_intersect.tab | cut -f 9 | sort | head -n 20

cat <(grep -P "\tgene\t" braker_annotation_x_all.pep.gff_intersect.tab | cut -f 9) <(grep -P "\tgene\t" braker_annotation_x_Assembled_transcripts.gtf_intersect.tab | cut -f 9) | sort -u | head -n 20

cat <(grep -P "\tgene\t" braker_annotation_x_all.pep.gff_intersect.tab | cut -f 9) <(grep -P "\tgene\t" braker_annotation_x_Assembled_transcripts.gtf_intersect.tab | cut -f 9) | sort -u | wc -l

# List of reliable genes (union)
cat <(grep -P "\tgene\t" braker_annotation_x_all.pep.gff_intersect.tab | cut -f 9) <(grep -P "\tgene\t" braker_annotation_x_Assembled_transcripts.gtf_intersect.tab | cut -f 9) | sort -u > reliable_genes.txt

wc -l reliable_genes.txt # 23266 (previous run: 22583)
```

## Move this old version of intersect files

In the end of this script, I use new version of intersect files, from
which I make “reliability table”.

``` sh
cd /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/intersects

mkdir -p intersects_old_version
mv -t ./intersects_old_version *
```

### Results of intersects

    Cardamine:
    Genes predicted by Braker                                                                              24828
    Genes predicted by Braker supported by aligned proteins                                                21160
    Genes predicted by Braker supported by transcripts assembled from RNAseq                               20404
    Genes predicted by Braker supported by both aligned proteins and transcripts assembled from RNAseq     18298
    Genes predicted by Braker supported by either aligned proteins or transcripts assembled from RNAseq    23266

    Cardamine (older run without RNA from flowers and flower buds):
    Genes predicted by Braker                                                                              24376
    Genes predicted by Braker supported by aligned proteins                                                20990
    Genes predicted by Braker supported by transcripts assembled from RNAseq                               18050
    Genes predicted by Braker supported by both aligned proteins and transcripts assembled from RNAseq     16457
    Genes predicted by Braker supported by either aligned proteins or transcripts assembled from RNAseq    22583

    Noccaea:
    Genes predicted by Braker                                                                              30373
    Genes predicted by Braker supported by aligned proteins                                                20757
    Genes predicted by Braker supported by transcripts assembled from RNAseq                               21454
    Genes predicted by Braker supported by both aligned proteins and transcripts assembled from RNAseq     17177
    Genes predicted by Braker supported by either aligned proteins or transcripts assembled from RNAseq    25034

    Alyssum:
    Genes predicted by Braker                                                                              32073
    Genes predicted by Braker supported by aligned proteins                                                21660
    Genes predicted by Braker supported by transcripts assembled from RNAseq                               21841
    Genes predicted by Braker supported by both aligned proteins and transcripts assembled from RNAseq     16978
    Genes predicted by Braker supported by either aligned proteins or transcripts assembled from RNAseq    26523

## Comparing predicted genes with previous run

I will compare genes predicted by Braker for the current run and the
previous run without RNA from flowers and flower buds.

``` sh
cd /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/intersects

module load bedtools2/2.30.0-gcc-10.2.1-5acjqve

bedtools intersect -a ../results_braker_01/braker.gtf -b ../../Cardamine_glauca_2024_12/results_braker_01/braker.gtf -u -f 1.0 -r > braker_current_and_old_run_intersect.tab

grep -P -c "\tgene\t" braker_current_and_old_run_intersect.tab # 23480
```

The current run of Braker predicted 24828 genes, the previous one
without RNA from flowers and flower buds predicted 24376 genes. Of
these, 23480 have exactly the same coordinates.

## Conversion of Braker GTF to GFF

There should be script for conversion in Augustus, but it seems that it
does not work well
(<https://github.com/Gaius-Augustus/BRAKER/issues/275>). Thus, I will
rather use AGAT.

``` sh
# interactive job
qsub -I -l select=1:ncpus=1:mem=8gb:scratch_local=10gb -l walltime=2:00:00

cd /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/
mkdir -p annot_processing
cd annot_processing

# run the container
singularity run /storage/brno12-cerit/home/duchmil/SW/agat/agat_1.4.0--pl5321hdfd78af_0.sif

agat_convert_sp_gxf2gxf.pl --gff ../results_braker_01/braker.gtf -o braker.gff
```

# Adding unknown expressed features

Bedtools subtract
<https://bedtools.readthedocs.io/en/latest/content/tools/subtract.html>
I will use -A option to get only transcripts with no overlap with Braker
annotation.

### Note for future

Some tools used to annotate noncoding genes in A. thaliana genomes are
descibed in this publication:

Lian, Qichao, Bruno Huettel, Birgit Walkemeier, Baptiste Mayjonade,
Céline Lopez-Roques, Lisa Gil, Fabrice Roux, Korbinian Schneeberger, and
Raphael Mercier. “A Pan-Genome of 69 Arabidopsis Thaliana Accessions
Reveals a Conserved Genome Structure throughout the Global Species
Range.” Nature Genetics 56, no. 5 (May 2024): 982–91.
<https://doi.org/10.1038/s41588-024-01715-9>.

It might be good to try them.

``` sh
cd /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/annot_processing


module load bedtools2/2.30.0-gcc-10.2.1-5acjqve

ls ../rnaseq/4_assembled_transcripts/

# Get assembled transcripts with no overlap with Braker annotation. I will use them as unknown expressed features.
bedtools subtract -a ../rnaseq/4_assembled_transcripts/Assembled_transcripts.gtf -b  braker.gff -A > Assembled_transcripts_with_no_annotation_overlap.gtf

grep -P -c "\ttranscript\t" Assembled_transcripts_with_no_annotation_overlap.gtf # 10991 (previous run: 7903)

# There are also exons without parental transcripts (when only part of the transcript overlaps with annotation). I will have to remove them.

# make a keep list with IDs of transcripts
grep -P "\ttranscript\t" Assembled_transcripts_with_no_annotation_overlap.gtf | sed -E 's/(^.*; transcript_id ")|("; cov.*$)//g' > Assembled_transcripts_with_no_annotation_overlap.txt

wc -l Assembled_transcripts_with_no_annotation_overlap.txt # 10991

### Filtering assembled transcripts based on keep list
# This should ensure that there will be no orphan exons.

# interactive job
qsub -I -l select=1:ncpus=1:mem=8gb:scratch_local=10gb -l walltime=2:00:00

cd /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/annot_processing


# run the container
singularity run /storage/brno12-cerit/home/duchmil/SW/agat/agat_1.4.0--pl5321hdfd78af_0.sif

# Filter based on keep list and convert to GFF
agat_sp_filter_feature_from_keep_list.pl --type transcript,RNA --gff Assembled_transcripts_with_no_annotation_overlap.gtf --keep_list Assembled_transcripts_with_no_annotation_overlap.txt -o Assembled_transcripts_with_no_annotation_overlap_almost_clean.gff

grep -P -c "\ttranscript\t" Assembled_transcripts_with_no_annotation_overlap_almost_clean.gff # 10991
grep -P -c "\tRNA\t" Assembled_transcripts_with_no_annotation_overlap_almost_clean.gff # 188
# For genes where some transcripts overlap with Braker annotation and some not, sometimes the overlapping exons are kept by AGAT and are asigned a RNA feature. I will need to remove them in following steps.

# Make a kill list of IDs of RNA features that I want to remove.
grep -P "\tRNA\t" Assembled_transcripts_with_no_annotation_overlap_almost_clean.gff | sed -E 's/^.*;transcript_id=//g' > transcript_kill_list.txt

agat_sp_filter_feature_from_kill_list.pl --gff Assembled_transcripts_with_no_annotation_overlap_almost_clean.gff --kill_list transcript_kill_list.txt -o Assembled_transcripts_with_no_annotation_overlap_clean.gff

grep -P -c "\ttranscript\t" Assembled_transcripts_with_no_annotation_overlap_clean.gff # 10991
grep -P -c "\tRNA\t" Assembled_transcripts_with_no_annotation_overlap_clean.gff # 0
# Now the RNA features and their exons should be removed, there shouldn't be any transcripts overlapping with Braker annotation.
```

## Merging annotations

``` sh
# interactive job
qsub -I -l select=1:ncpus=1:mem=8gb:scratch_local=10gb -l walltime=2:00:00

cd /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/annot_processing


# run the container
singularity run /storage/brno12-cerit/home/duchmil/SW/agat/agat_1.4.0--pl5321hdfd78af_0.sif

# Complement annotations
# The file with unknown expressed features will be used as a reference, so the header of this file will be kept.
agat_sp_complement_annotations.pl --ref Assembled_transcripts_with_no_annotation_overlap_clean.gff --add braker.gff -o merged_annotation_01.gff

grep -P -c "\tgene\t" merged_annotation_01.gff # 33729
grep -P -c "\ttranscript\t" merged_annotation_01.gff # 38864
grep -P -c "\tRNA\t" merged_annotation_01.gff # 0
grep -P -c "\tmRNA\t" merged_annotation_01.gff # 1136


# Check numbers of lines (without header)
grep -P -v -c "^#" braker.gff # 542601
grep -P -v -c "^#" Assembled_transcripts_with_no_annotation_overlap_clean.gff # 46404
grep -P -v -c "^#" merged_annotation_01.gff # 589005
```

``` r
# Check whether the number of lines of merged annotation is sum of the files that were merged
542601 + 46404 == 589005 # TRUE
```

# Changing IDs in R

``` r
setwd("D:/!ecolgen/annotations/Cardamine_glauca_2025_06")

# read the gff file
gff.1 <- read.table(file = "annot_processing/merged_annotation_01.gff", header = F, sep = "\t", comment.char = "#" 
                     #, nrows = 1500
                    )
summary(gff.1)
head(gff.1)
gff.1[1:50, ]
```

### What feature types are there?

``` r
# Types of features
levels(as.factor(gff.1$V3))
table(as.factor(gff.1$V3))

# scaffold names
levels(as.factor(gff.1$V1))
```

### Adjusting column 2 and 3

``` r
## Changes in column 3

# mRNA is there only for genes predicted by GeneMark.hmm3
gff.1[grepl(pattern = "mRNA|transcript", x = gff.1$V3), ][1:40, ]
# It seems that mRNA is redundant with transcript
gff.1[grepl(pattern = "mRNA|transcript", x = gff.1$V3) & grepl(pattern = "GeneMark.hmm3", x = gff.1$V2), ][1:40, ]
# It doesn't have any childs features.
gff.1[grepl(pattern = "GeneMark.hmm3", x = gff.1$V2), ][1:40, ]

# I will remove mRNA features.
gff.2 <- gff.1[gff.1$V3 != "mRNA", ]
dim(gff.1)
dim(gff.2)


## Changes in column 2

table(as.factor(gff.1$V2))
table(as.factor(gff.1$V2), as.factor(gff.1$V3))


# Genes added by AGAT to assembled transcript - change col 2 to "StringTie"
gff.2[gff.2[, 2] == "AGAT", 2] <- "StringTie"

# Add "BRAKER" before other sources
gff.2$V2[gff.2$V2 %in% c("AUGUSTUS", "GeneMark.hmm3", "gmst")] <- paste0("BRAKER-", gff.2$V2[gff.2$V2 %in% c("AUGUSTUS", "GeneMark.hmm3", "gmst")])

## other changes

# For features "gene", remove cov, fPKM, tPM and transcript_id, because AGAT copied it just from the first transcript
gff.2[grepl(pattern = "gene", x = gff.2$V3), 9][1:40]
gff.2[grepl(pattern = "gene", x = gff.2$V3) & gff.2$V2 == "StringTie", 9][1:40]
# For some genes, there is also "exon_number", which I will also remove.


gene.col9.1 <- gff.2[grepl(pattern = "gene", x = gff.2$V3) & gff.2$V2 == "StringTie", 9]
gene.col9.2 <- gsub(pattern = ";*cov=[0-9.]+;*", replacement = ";", x = gene.col9.1)
gene.col9.3 <- gsub(pattern = ";*fPKM=[0-9.]+;*", replacement = ";", x = gene.col9.2)
gene.col9.4 <- gsub(pattern = ";*tPM=[0-9.]+;*", replacement = ";", x = gene.col9.3)
gene.col9.5 <- gsub(pattern = ";*exon_number=[0-9.]+;*", replacement = ";", x = gene.col9.4)
gene.col9.6 <- gsub(pattern = ";*transcript_id=STRG[0-9.]+;*", replacement = "", x = gene.col9.5)
gene.col9.7 <- gsub(pattern = ";*gene_id=STRG[0-9.]+;*", replacement = "", x = gene.col9.6)
gff.2[grepl(pattern = "gene", x = gff.2$V3) & gff.2$V2 == "StringTie", 9] <- gene.col9.7
```

### Changing IDs

``` r
# old gene names
old.genes <- gsub(pattern = "ID=|;gene_id.*$", replacement = "", x = gff.2[grepl(pattern = "gene", x = gff.2$V3), "V9"])
scaff.genes <- as.integer(gsub(pattern = "scaffold_", replacement = "", x = gff.2[grepl(pattern = "gene", x = gff.2$V3), "V1"]))
levels(as.factor(scaff.genes))

# New gene names with species code and readable ID
# COmpared to the very first version, the IDs are one digit shorter (I am adding just one 0 in the end).
new.genes <- paste0("Cg", formatC(x = scaff.genes, width = 3, flag = "0"), "G", formatC(x = 1:length(old.genes), width = 5, flag = "0"), "0")
tail(new.genes)

# Checking if all gene names have the same length
table(nchar(new.genes))

genes <- cbind.data.frame(old.genes, new.genes, scaff.genes)
head(genes)
tail(genes)


# Row numbers where data for one gene begin and end are needed (see the loop below)
# GFF needs to be ordered (rows for one gene and its subfeatures should't be interspersed with data for other gene).
genes$rows <- grep(pattern = "gene", x = gff.2$V3)
genes$end.rows <- NA
genes$end.rows[1:(nrow(genes)-1)] <- genes$rows[2:nrow(genes)]-1
genes$end.rows[nrow(genes)] <- nrow(gff.2)
tail(genes)

# New gff
gff.3 <- gff.2

## Loop to make new feature ids in GFF

# Gene names will be changed to new ones.
# Transcript names will be in the form gene_name.t1.
# Other features ("exon", "intron", "cds", "start_codon", "stop_codon") will have IDs in form "exon-AM_transcript_id-1", where the last number will count exons in that particular transcript.

# Warning:
# This loops expect the GFF to be ordered (rows for one gene and its subfeatures should't be interspersed with data for other gene, the same should be valid for transcripts of one gene). First, I tried to make it more universal, but it was always searching for the rows for one gene using grep, which was slow like a hell.

i=1
i=5

for(i in 1:nrow(genes)) {
  # Extract rows for one gene and its subfeatures.
  descr1 <- gff.2[genes$rows[i]:genes$end.rows[i], "V9"]
  
  # Change the old gene ID for new gene ID
  descr2 <- gsub(pattern = genes$old.genes[i], 
                 replacement = genes$new.genes[i], 
                 x = descr1)
  
  # For stringtie assembled transcripts, add "t" to make the transcript ID in form "Np001G000010.t1"
  descr2 <- gsub(pattern = paste0("(", genes$new.genes[i], "\\.)([[:digit:]]+)"), 
                 replacement = "\\1t\\2", 
                 x = descr2)
  
  # # Patterns defining how to find old gene name within different context.
  # patterns <- c(paste0("=", genes$old.genes[i], ";"),
  #               paste0("=", genes$old.genes[i], "\\.t"),
  #               paste0("=", genes$old.genes[i], "$")#,
  #               #"ID=agat-exon-[[:digit:]]+;"
  #               )
  # 
  # # Strings defining what to exchange the previous patterns for
  # replacements <- c(paste0("=", genes$new.genes[i], ";"),
  #                       paste0("=", genes$new.genes[i], ".t"),
  #                       paste0("=", genes$new.genes[i])#,
  #                   #paste0("ID=exon-", genes$new.genes[i], "-#;")
  #                   )
  # 
  # descr2 <- descr1
  # 
  # # Loop exchanging patterns for replacements
  # for(j in 1:length(patterns)) {
  #   descr2 <- gsub(pattern = patterns[j], 
  #                replacement = replacements[j], 
  #                x = descr2)
  # }
  # 
  # Loop to change IDs of features other than genes or transcripts
  for(feature in c("exon", "intron", "cds", "start_codon", "stop_codon")) {
    # Extract rows for one feature
    extr.feature <- descr2[grep(pattern = paste0("ID=(agat|IDmodified)-", feature, "-[[:digit:]]+;"), x = descr2)]
    # For some genes, there are no introns and it would make mess if we try to exchange something
    if(length(extr.feature) > 0) {
      # Transcript IDs
      transcripts <- gsub(pattern = ".*transcript_id=", replacement = "", x = extr.feature)
      # Sequence of numbers for features within transcript and gene (1 to n for every transcript)
      # feat.numbers <- unlist(mapply(seq, from = 1, to = table(transcripts))) # this doesn't work if the exons or other features are not in alphabetical order
      feat.numbers <- ave(x = seq_along(transcripts), transcripts, FUN = seq_along)
      # Modify the feature IDs
      mod.feature <- mapply(sub, 
                            pattern = paste0("ID=(agat|IDmodified)-", feature, "-[[:digit:]]+;"), 
                            replacement = paste0("ID=", feature, "-", transcripts, "-", feat.numbers, ";"), 
                            x = extr.feature, USE.NAMES = F)
      # Replace rows with modified ones
      descr2[grep(pattern = paste0("ID=(agat|IDmodified)-", feature, "-[[:digit:]]+;"), x = descr2)] <- mod.feature
    }
  }
  
  # Replace info in column 9 of GFF for the modified one
  gff.3[genes$rows[i]:genes$end.rows[i], "V9"] <- descr2
  # grep version (slow)
  # gff.3[grep(pattern = paste0("=", genes$old.genes[i], "(;|\\.t|$)"), x = gff.2$V9), "V9"] <- descr2
  
  # For every 100 genes, print the counting
  if(i %% 100 == 0) print(paste(i, "genes out of", nrow(genes), "done"))
  # i=i+1
}


tail(gff.3, n = 40)
# tail(gff.3, n = 100)
# 
# gff.3[grep(pattern = "=g365(;|\\.t|$)", x = gff.2$V9), ] 
#   
# i=1
# i=2
```

### Further adjustments in column 3

``` r
gff.4 <- gff.3

# Changing feature names
gff.4$V3[gff.4$V3 == "gene" & gff.4$V2 == "StringTie"] <- "ncRNA_gene"
gff.4$V3[gff.4$V3 == "transcript" & gff.4$V2 == "StringTie"] <- "ncRNA"



levels(as.factor(gff.4$V3))
table(gff.4$V3)
```

### Preparing header

It is needed to produce statistics (that are computed later) and then to
get back here and update the statistics in the header.

``` r
# read the header of gff file
gff.raw <- readLines(con = "annot_processing/merged_annotation_01.gff" 
                     , n = 1500
                    )
gff.header <- gff.raw[grep(pattern = "^#", x = gff.raw)]

# new header with new date
gff.header.3 <- c("##gff-version 3",
                  "# Cardamine glauca genome annotation",
                  "# Version 1.1",
                  "# 2025-06-10",
                  "#",
                  "# Version history",
                  "#",
                  "# Version 1.1 (2025-06-10)",
                  "# - Included RNAseq data from flowers and flower buds.",
                  "# - Removed one zero from the end of the gene IDs.",
                  "# - The genome assembly was renamed, but it is the same as for version 1.",
                  "# - Beware that the gene IDs might not correspond to the gene IDs of previous version.",
                  "#",
                  "# Version 1 (2024-12-30)",
                  "# - Initial version",
                  "#",
                  "# Genome assembly",
                  "#",
                  "# File Cardamine_glauca_CUNI_V1_2024_10_masked.fa",
                  "# # contigs (>= 0 bp)         53",
                  "# # contigs (>= 1000 bp)      53",
                  "# # contigs (>= 5000 bp)      53",
                  "# # contigs (>= 10000 bp)     53",
                  "# # contigs (>= 25000 bp)     53",
                  "# # contigs (>= 50000 bp)     45",
                  "# Total length (>= 0 bp)      301015024",
                  "# Total length (>= 1000 bp)   301015024",
                  "# Total length (>= 5000 bp)   301015024",
                  "# Total length (>= 10000 bp)  301015024",
                  "# Total length (>= 25000 bp)  301015024",
                  "# Total length (>= 50000 bp)  300693161",
                  "# # contigs                   53",
                  "# Largest contig              39236745",
                  "# Total length                301015024",
                  "# GC (%)                      37.85",
                  "# N50                         34778859",
                  "# N75                         30582462",
                  "# L50                         4",
                  "# L75                         7",
                  "# # N's per 100 kbp           0.80",
                  "#",
                  "# Annotation",
                  "#",
                  "# Annotation generated by Milos Duchoslav (Group of Filip Kolar, ",
                  "# Department of Botany, Faculty of Science, Charles University, Prague, Czechia).",
                  "# Email: duchmil[at]gmail.com",
                  "# ",
                  "# Protein coding genes were predicted by",
                  "# BRAKER version 3.0.8 based on RNA-seq data from Cardamine glauca ",
                  "# and OrthoDB proteins from the whole Viridiplantae.",
                  "# ",
                  "# To include non-coding RNA genes and similar features,",
                  "# transcripts assembled from RNA-seq data using StringTie v2.2.2",
                  "# were filtered using bedtools v2.30.0 to keep only those that",
                  "# do not have any overlap with protein-coding gene annotations.",
                  "# They were included in annotation as 'ncRNA_gene' and 'ncRNA'.",
                  "# However, they are simply regions that are transcribed and ",
                  "# probably are not protein-coding genes.",
                  "# ",
                  "# GTF converted to GFF using AGAT v1.4.0.",
                  "# Feature IDs changed using custom R script.",
                  "# ",
                  "# Message for Vitek from Milos: If you find this message, I will buy you a beer.",
                  "# ",
                  "# Annotation statistics",
                  "# ",
                  "# Type (3rd column)       Number  Size total (kb)  Size mean (bp)  % of the genome",
                  "# cds                     153637         36170.64          235.43            12.02",
                  "# exon                    180149         46155.78          256.21            15.33",
                  "# gene                     24828         53008.90         2135.05            17.61",
                  "# intron                  125764         26300.04          209.12             8.74",
                  "# ncrna                    10991         29818.57         2713.00             9.91",
                  "# ncrna_gene                8901         20901.72         2348.24             6.94",
                  "# start_codon              27865            83.57            3.00             0.03",
                  "# stop_codon               27861            83.57            3.00             0.03",
                  "# transcript               27873         62470.68         2241.26            20.75",
                  "# Total                   587869        274993.47          467.78            91.36",
                  "# "
)

# gff.header.2 <- gff.header
```

### Writing GFF file

``` r
# write header to file
write.table(x = gff.header.3, 
            file = "annot_processing/Cardamine_glauca_CUNI_V1_annotation_v1.1.gff", 
            sep = "\t",
            row.names = F, col.names = F, quote = F)
# append the gff itself
write.table(x = gff.4, 
            file = "annot_processing/Cardamine_glauca_CUNI_V1_annotation_v1.1.gff", 
            sep = "\t",
            row.names = F, col.names = F, quote = F, append = T)
```

## Statistics of annotation using AGAT

``` sh
# interactive job
qsub -I -l select=1:ncpus=1:mem=8gb:scratch_local=10gb -l walltime=2:00:00

mkdir -p /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/stat_annotation

cd /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/stat_annotation

# run the container
singularity run /storage/brno12-cerit/home/duchmil/SW/agat/agat_1.4.0--pl5321hdfd78af_0.sif

# Name of genome assembly
genome_assembly="Cardamine_glauca_CUNI_V1_2024_10_masked.fa"

# basic statistics
agat_sq_stat_basic.pl -i ../annot_processing/Cardamine_glauca_CUNI_V1_annotation_v1.1.gff -g ../genome_assembly/$genome_assembly > stat_basic_Cardamine_glauca_CUNI_V1_annotation_v1.1.gff.txt

# detailed statistics
agat_sp_statistics.pl -i ../annot_processing/Cardamine_glauca_CUNI_V1_annotation_v1.1.gff -g ../genome_assembly/$genome_assembly > stat_detailed_Cardamine_glauca_CUNI_V1_annotation_v1.1.gff.txt

# view the basic statistics
column -t stat_basic_Cardamine_glauca_CUNI_V1_annotation_v1.1.gff.txt
tail -n +11 stat_basic_Cardamine_glauca_CUNI_V1_annotation_v1.1.gff.txt | column -s $'\t' -t -R 2,3,4,5
```

    Type (3rd column)       Number  Size total (kb)  Size mean (bp)  % of the genome
    cds                     153637         36170.64          235.43            12.02
    exon                    180149         46155.78          256.21            15.33
    gene                     24828         53008.90         2135.05            17.61
    intron                  125764         26300.04          209.12             8.74
    ncrna                    10991         29818.57         2713.00             9.91
    ncrna_gene                8901         20901.72         2348.24             6.94
    start_codon              27865            83.57            3.00             0.03
    stop_codon               27861            83.57            3.00             0.03
    transcript               27873         62470.68         2241.26            20.75
    Total                   587869        274993.47          467.78            91.36

    The previous run (without RNAseq from flowers and flower buds)

    Type (3rd column)  Number  Size total (kb)  Size mean (bp)  % of the genome 
    cds                152346         36040.17          236.57            11.97
    exon               171367         43355.25          253.00            14.40
    gene                24376         52144.58         2139.18            17.32
    intron             124959         25315.38          202.59             8.41
    ncrna                7903         21370.51         2704.10             7.10
    ncrna_gene           6482         15411.90         2377.65             5.12
    start_codon         27378            82.11            3.00             0.03
    stop_codon          27380            82.12            3.00             0.03
    transcript          27387         61355.56         2240.32            20.38
    Total              569578        255157.58          447.98            84.77

## Extraction of protein and coding sequences using AGAT

``` sh
cd /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/annot_processing

# run the container
singularity run /storage/brno12-cerit/home/duchmil/SW/agat/agat_1.4.0--pl5321hdfd78af_0.sif

# Name of genome assembly
genome_assembly="Cardamine_glauca_CUNI_V1_2024_10_masked.fa"

# Protein sequences
agat_sp_extract_sequences.pl -g Cardamine_glauca_CUNI_V1_annotation_v1.1.gff -f ../genome_assembly/$genome_assembly -p -o Cardamine_glauca_CUNI_V1_annotation_v1.1_proteins.fasta
# CDS
agat_sp_extract_sequences.pl -g Cardamine_glauca_CUNI_V1_annotation_v1.1.gff -f ../genome_assembly/$genome_assembly -t cds -o Cardamine_glauca_CUNI_V1_annotation_v1.1_cds.fasta

# counting the number of sequences
grep -c ">" Cardamine_glauca_CUNI_V1_annotation_v1.1_proteins.fasta # 27873
grep -c ">" Cardamine_glauca_CUNI_V1_annotation_v1.1_cds.fasta # 27873

# counting internal stop codons
# converting from folded fasta to unfolded fasta for better counting and checking for internal stop codons
awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' < Cardamine_glauca_v1_annotation_proteins.fasta | grep \*[[:alpha:]] | wc -l
# 0
```

# Final files

Move the final files.

``` sh
cd /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/annot_processing

# compress the final files by gzip
gzip --keep Cardamine_glauca_CUNI_V1_annotation_v1.1.gff Cardamine_glauca_CUNI_V1_annotation_v1.1_proteins.fasta Cardamine_glauca_CUNI_V1_annotation_v1.1_cds.fasta

# directory for final files
mkdir -p /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/final_files

# move the final files
mv -t /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/final_files Cardamine_glauca_CUNI_V1_annotation_v1.1.gff.gz Cardamine_glauca_CUNI_V1_annotation_v1.1_proteins.fasta.gz Cardamine_glauca_CUNI_V1_annotation_v1.1_cds.fasta.gz


# compress genome assembly and add it to final files
cd /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/genome_assembly
genome_assembly="Cardamine_glauca_CUNI_V1_2024_10_masked.fa"

gzip --keep $genome_assembly
mv -t /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/final_files $genome_assembly.gz
```

# Reliability of the predicted genes

*This part of the script was added later according to script fo
Odontarrhena.*

This is one of the quality checks. I will produce a table that will show
how reliable the predicted genes are (if they have support in RNAseq
data or aligned proteins from other species).

## Intersect between gene models and aligned proteins

bedtools intersect
<https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html>

``` sh
cd /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06
mkdir -p intersects
cd intersects

module load bedtools2/2.30.0-gcc-10.2.1-5acjqve

ls ../protein_seqs_input/2_aligned/

# Intersects with proteins from single species
bedtools intersect -a ../annot_processing/Cardamine_glauca_CUNI_V1_annotation_v1.1.gff -b ../protein_seqs_input/2_aligned/A.thaliana.pep.gff -u -f 0.9 -r > braker_annotation_x_A.thaliana.pep.gff_intersect.tab

bedtools intersect -a ../annot_processing/Cardamine_glauca_CUNI_V1_annotation_v1.1.gff -b ../protein_seqs_input/2_aligned/A.lyrata.pep.gff -u -f 0.9 -r > braker_annotation_x_A.lyrata.pep.gff_intersect.tab

bedtools intersect -a ../annot_processing/Cardamine_glauca_CUNI_V1_annotation_v1.1.gff -b ../protein_seqs_input/2_aligned/B.rapa.pep.gff -u -f 0.9 -r > braker_annotation_x_B.rapa.pep.gff_intersect.tab

# Intersect with all proteins together
bedtools intersect -a ../annot_processing/Cardamine_glauca_CUNI_V1_annotation_v1.1.gff -b ../protein_seqs_input/2_aligned/*  -u -f 0.9 -r > braker_annotation_x_all.pep.gff_intersect.tab
```

This command finds the intersects between Braker output (`-a`) and at
least one of the aligned proteins (files in `-b`). It will report the
genes in `-a` only once (option `-u`). The overlap should be at least 90
% (`-f 0.9`) for both `-a` and `-b` (option `-r` as reciprocal).

## Intersect between gene models and assembled transcripts

Again bedtools intersect. This time we do not use the `-r` option - the
overlap should be at least 90 % of Braker output, but not necesarilly of
the Stringtie output. It is because the Stringtie tries to predict whole
transcripts including UTRs.

``` sh
cd /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/intersects

module load bedtools2/2.30.0-gcc-10.2.1-5acjqve

ls ../rnaseq/4_assembled_transcripts/

# Intersects with Stringtie assembled transcripts
bedtools intersect -a ../annot_processing/Cardamine_glauca_CUNI_V1_annotation_v1.1.gff -b ../rnaseq/4_assembled_transcripts/Assembled_transcripts.gtf -u -f 0.9 > braker_annotation_x_Assembled_transcripts.gtf_intersect.tab

# Number of genes predicted by Braker
grep -P -c "\tgene\t" ../annot_processing/Cardamine_glauca_CUNI_V1_annotation_v1.1.gff # 24828
# Number of transcripts assembled from RNAseq
grep -P -c "\ttranscript\t" ../rnaseq/4_assembled_transcripts/Assembled_transcripts.gtf # 45030
# Number of genes assembled from RNAseq
# (Gene features are not annotated separately by StringTie, but there is a field gene_id for transcripts which can be used.)
grep -P "\ttranscript\t" ../rnaseq/4_assembled_transcripts/Assembled_transcripts.gtf | sed -E 's/.*\tgene_id//' | sed -E 's/; transcript_id.*//' | sort -u | wc -l # 31004
# Intersect of genes predicted by Braker and transcripts assembled from RNAseq
grep -P -c "\tgene\t" braker_annotation_x_Assembled_transcripts.gtf_intersect.tab # 20404
```

**Note** There is a problem that some of the assembled transcripts have
long introns. If the whole gene predicted by Braker falls within intron
of assembled transcript, it is still reported as supported, even if no
reads map to this region. I don’t have any easy solution how to solve
this for now.

### Checking the outputs

Here I look at the number of genes that are reported to have intersect
with some aligned protein. There is a possibility that only some shorter
transcript (splicing variant) of the particular gene have overlap big
enough to be reported. Later I am putting into the reliability table
also these transcripts and their particular genes, so the numbers might
be slightly higher.

``` sh
grep -P -c "\tgene\t" ../annot_processing/Cardamine_glauca_CUNI_V1_annotation_v1.1.gff
# grep -P -c "\tmRNA\t" ../protein_seqs_input/2_aligned/A.thaliana.pep.gff # 47853
# grep -P -c "\tgene\t" braker_annotation_x_A.thaliana.pep.gff_intersect.tab # 16967

for FILE in braker_annotation_x_*.tab
do
COUNT=$(grep -P -c "\tgene\t" $FILE)
echo "$FILE $COUNT"
done
```

## Generating table of reliability of genes

Table will be generated from intersects of genes predicted by Braker and
aligned RNAseq data and proteins from other Brassicaceae species. It
will show which genes are supported by additional data and which not.

Using only protein coding genes for now. Also other genes could be used
in the future.

``` sh
cd /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/intersects

# checks
# grep -P "\tgene\t" braker_annotation_x_all.pep.gff_intersect.tab | cut -f 9 | sed 's/ID=//' | sort | head -n 20
# this will use genes belonging to all transcripts that were supported
grep -P "\ttranscript\t" braker_annotation_x_all.pep.gff_intersect.tab | cut -f 9 | sed 's/^.*Parent=//' | sort --unique | head -n 20
# Use also nc_RNA genes? Not now...
grep -P "\t(transcript|ncRNA)\t" braker_annotation_x_all.pep.gff_intersect.tab | cut -f 9 | sed 's/^.*Parent=//' | sed 's/;.*$//' | less

## Get the lists of genes from intersects
# This will extract genes belonging to all transcripts that were supported.
# It might be different from the case when I take just the "gene" features from GFF - there could be supported just one transcript and not the entire gene.
for FILE in braker_annotation_x_*_intersect.tab
do
NEW_FILE=$(echo $FILE | sed 's/.tab$/_genes.txt/')
# Extract transcripts and get the gene IDs from them.
grep -P "\ttranscript\t" $FILE | cut -f 9 | sed 's/^.*Parent=//' | sort --unique > $NEW_FILE
# number of genes extracted
echo $(wc -l $NEW_FILE)
done

## Get the list of genes from the annotation
# Just protein coding genes.
grep -P "\tgene\t" ../annot_processing/Cardamine_glauca_CUNI_V1_annotation_v1.1.gff | cut -f 9 | sed 's/ID=//' | sed 's/\r$//' | sort > protein_coding_genes.txt
wc -l protein_coding_genes.txt


## AWK script generating the table
# I will use also the list of all protein coding genes (to have a full list of genes) and then I will remove the respective column, where will be only ones and no zeros.
awk 'FNR==1{
    f++
    fname = FILENAME
    gsub(/^braker_annotation_x_/, "", fname)   # remove prefix
    gsub(/\.g[tf]f_intersect_genes\.txt$/, "", fname) # remove suffix
    files[f] = fname
}
{
    a[$1][f] = 1
}
END{
    # header
    printf "Gene"
    for(i=1;i<=f;i++) printf "\t%s", files[i]
    print ""

    # collect all gene names
    for(g in a) genes[++n] = g

    # sort gene names
    n = asort(genes)

    # print rows
    for(j=1;j<=n;j++){
        g = genes[j]
        printf "%s", g
        for(i=1;i<=f;i++) printf "\t%d", (a[g][i]?1:0)
        print ""
    }
}' protein_coding_genes.txt braker_annotation_x_Assembled_transcripts.gtf_intersect_genes.txt braker_annotation_x_A.lyrata.pep.gff_intersect_genes.txt braker_annotation_x_A.thaliana.pep.gff_intersect_genes.txt braker_annotation_x_B.rapa.pep.gff_intersect_genes.txt | 
     # remove the column with all genes ("protein_coding_genes.txt")
     cut -f2 --complement > Cardamine_glauca_CUNI_V1_annotation_v1.1_protein_coding_genes_support.tsv

# check - this should be number of gene features in final gff plus one for header
wc -l Cardamine_glauca_CUNI_V1_annotation_v1.1_protein_coding_genes_support.tsv

# Copying the reliability table to final files
cp -v Cardamine_glauca_CUNI_V1_annotation_v1.1_protein_coding_genes_support.tsv /storage/brno12-cerit/home/duchmil/annotations/Cardamine_glauca_2025_06/final_files/
```
