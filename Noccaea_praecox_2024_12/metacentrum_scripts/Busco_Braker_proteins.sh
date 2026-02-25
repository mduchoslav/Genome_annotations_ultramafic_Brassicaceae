### Script for Metacentrum

#PBS -N Busco_Braker_proteins
#PBS -l select=1:ncpus=4:mem=8gb:scratch_local=1000gb
#PBS -l walltime=2:00:00 
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
cp -v $DATADIR/results_braker_01/braker.aa $SCRATCHDIR || { echo >&2 "Error while copying index file(s)!"; exit 2; }

echo "Input files copied." | ts '[%Y-%m-%d %H:%M:%S]'


# move into scratch directory
cd $SCRATCHDIR 


# activate Busco
source /storage/brno2/home/duchmil/SW/mambaforge/bin/activate busco_5_7_1

# run Busco
busco --in braker.aa --mode proteins --lineage_dataset brassicales_odb10 --cpu 4

echo "Main computation done." | ts '[%Y-%m-%d %H:%M:%S]'

# move the output to user's DATADIR or exit in case of failure
cp -vR * $DATADIR/busco_results || { echo >&2 "Result file(s) copying failed (with a code $?) !!"; exit 4; }

echo "Output files copied." | ts '[%Y-%m-%d %H:%M:%S]'

# clean the SCRATCH directory
clean_scratch