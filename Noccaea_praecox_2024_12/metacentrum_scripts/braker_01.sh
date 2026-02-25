### Script for Metacentrum

#PBS -N braker_Noccaea_01
#PBS -l select=1:ncpus=16:mem=96gb:scratch_local=1000gb
#PBS -l walltime=48:00:00 
#PBS -m ae

# define a DATADIR variable: directory where the input files are taken from and where output will be copied to
DATADIR=/storage/brno12-cerit/home/duchmil/annotations/Noccaea_praecox_2024_12

# append a line to a file "jobs_info.txt" containing the ID of the job, the hostname of node it is run on and the path to a scratch directory
# this information helps to find a scratch directory in case the job fails and you need to remove the scratch directory manually 
echo "$PBS_JOBID is running on node `hostname -f` in a scratch directory $SCRATCHDIR" | ts '[%Y-%m-%d %H:%M:%S]' >> $PBS_O_WORKDIR/jobs_info.txt

# test if scratch directory is set
# if scratch directory is not set, issue error message and exit
test -n "$SCRATCHDIR" || { echo >&2 "Variable SCRATCHDIR is not set!"; exit 1; }

# copy files
cp -r $DATADIR/rnaseq/3_aligned_reads/RNAseq_trimmed_merged.bam $SCRATCHDIR || { echo >&2 "Error while copying input file(s)!"; exit 2; }
cp -r /storage/brno12-cerit/home/duchmil/annotations/OrthoDB_proteins/Viridiplantae.fa $SCRATCHDIR || { echo >&2 "Error while copying input file(s)!"; exit 2; }
cp -r $DATADIR/genome_assembly/Noccaea_polished_2024_10_Mahnaz.masked.fa $SCRATCHDIR || { echo >&2 "Error while copying input file(s)!"; exit 2; }


# move into scratch directory
cd $SCRATCHDIR 
mkdir results_braker_01

# running BRAKER
export BRAKER_SIF=/storage/brno12-cerit/home/duchmil/SW/braker_sw/braker3.sif

singularity exec -B ${PWD}:${PWD} ${BRAKER_SIF} braker.pl --bam=RNAseq_trimmed_merged.bam --genome=Noccaea_polished_2024_10_Mahnaz.masked.fa --prot_seq=Viridiplantae.fa --threads=16 --species=Noccaea_praecox --workingdir=$SCRATCHDIR/results_braker_01
# Note: Protein sequences shouldn't be compressed. It should be plain fasta.


# move the output to user's DATADIR or exit in case of failure
cp -r results_braker_01 $DATADIR/ || { echo >&2 "Result file(s) copying failed (with a code $?) !!"; exit 4; }

# clean the SCRATCH directory
clean_scratch