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
cp -v Assembled_transcripts.gtf $DATADIR/4_assembled_transcripts || { echo >&2 "Result file(s) copying failed (with a code $?) !!"; exit 4; }

echo "Output files copied." | ts '[%Y-%m-%d %H:%M:%S]'

# clean the SCRATCH directory
clean_scratch