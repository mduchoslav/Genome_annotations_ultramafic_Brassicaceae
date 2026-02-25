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