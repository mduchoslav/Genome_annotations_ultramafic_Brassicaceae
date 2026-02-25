### Script for Metacentrum

#PBS -N trim_galore_for_RNAseq
#PBS -l select=1:ncpus=4:mem=32gb:scratch_local=10gb:brno=True
#PBS -l walltime=24:00:00 
#PBS -m ae

## trim_galore

# define a DATADIR variable: directory where the input files are taken from
DATADIR=/storage/brno12-cerit/home/duchmil/annotations/Aethionema_saxatile_2025_06/rnaseq/1_raw_reads
# directory for output
OUTDIR=/storage/brno12-cerit/home/duchmil/annotations/Aethionema_saxatile_2025_06/rnaseq/2_trimmed_reads

# append a line to a file "jobs_info.txt" containing the ID of the job, the hostname of node it is run on and the path to a scratch directory
# this information helps to find a scratch directory in case the job fails and you need to remove the scratch directory manually 
echo "$PBS_JOBID is running on node `hostname -f` in a scratch directory $SCRATCHDIR" >> $PBS_O_WORKDIR/jobs_info.txt

# move into data directory
cd $DATADIR

# running trim_galore
module load trim_galore/0.6.2_py3
trim_galore --version # version 0.6.2

## Trimming
# There is needed "sort" because the manual says:
# "Trim Galore! expects paired-end files to be supplied in a pairwise fashion, e.g. file1_1.fq file1_2.fq SRR2_1.fq.gz SRR2_2.fq.gz ... ."

# All together
# time find $DATADIR -type f,l \( -iname "*.fastq.gz" -o -iname "*.fq.gz" \) | sort | xargs trim_galore --paired --cores 4 -o $OUTDIR

# Macrogen and Novogene separately (they might have different adapters and Trim Galore uses the first 1 million sequences of the first file to detect adapter type)
time find $DATADIR -type f,l -iname "*Macrogen_*.fastq.gz" | sort | xargs trim_galore --paired --cores 4 -o $OUTDIR
time find $DATADIR -type f,l -iname "*.fq.gz" | sort | xargs trim_galore --paired --cores 4 -o $OUTDIR


# Remove trim_galore module
# The MultiQC otherwise does not work.
module remove trim_galore/0.6.2_py3

## FastQC and MultiQC after trim_galore

#### Files after trim_galore trimming end with *.fq.gz and not *.fastq.gz!

# report name
report_name="Aethionema_RNAseq_MultiQC_report_after_trimming"

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