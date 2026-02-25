### Script for Metacentrum

#PBS -N Busco_proteins_CDS_assembly_transcripts
#PBS -l select=1:ncpus=4:mem=32gb:scratch_local=1000gb
#PBS -l walltime=4:00:00 
#PBS -m ae

# define a DATADIR variable: directory where the input files are taken from and where output will be copied to
DATADIR=/storage/brno12-cerit/home/duchmil/annotations/Odontarrhena_muralis_2025_10
# Name of genome assembly
genome_assembly="Odontarrhena_muralis_CUNI_V1_2025_05_masked.fa"

# append a line to a file "jobs_info.txt" containing the ID of the job, the hostname of node it is run on and the path to a scratch directory
# this information helps to find a scratch directory in case the job fails and you need to remove the scratch directory manually 
echo "$PBS_JOBID is running on node `hostname -f` in a scratch directory $SCRATCHDIR" | ts '[%Y-%m-%d %H:%M:%S]' >> $PBS_O_WORKDIR/jobs_info.txt

# test if scratch directory is set
# if scratch directory is not set, issue error message and exit
test -n "$SCRATCHDIR" || { echo >&2 "Variable SCRATCHDIR is not set!"; exit 1; }

# copy files
cp -v $DATADIR/results_braker_01/braker.aa $SCRATCHDIR || { echo >&2 "Error while copying input file(s)!"; exit 2; }
cp -v $DATADIR/results_braker_01/braker.codingseq $SCRATCHDIR || { echo >&2 "Error while copying input file(s)!"; exit 2; }
cp -v $DATADIR/genome_assembly/$genome_assembly $SCRATCHDIR || { echo >&2 "Error while copying input file(s)!"; exit 2; }
cp -v $DATADIR/rnaseq/4_assembled_transcripts/Assembled_transcripts.fasta $SCRATCHDIR || { echo >&2 "Error while copying input file(s)!"; exit 2; }



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

## Plot Busco results

cd $DATADIR/busco_results/

mkdir -p summaries_BUSCO

# copy summaries from all Busco results folders
find BUSCO* -name "short_summary.specific.*.txt" -exec cp {} summaries_BUSCO/ \;

# Rename some summaries (avoid dot in specific part of the name) so that Busco script will distinguish them
cd summaries_BUSCO/
mv short_summary.specific.brassicales_odb10.BUSCO_braker.aa.txt short_summary.specific.brassicales_odb10.BUSCO_braker_aa.txt
mv short_summary.specific.brassicales_odb10.BUSCO_braker.codingseq.txt short_summary.specific.brassicales_odb10.BUSCO_braker_codingseq.txt

# Use script to generate plot
generate_plot.py -wd .

echo "Plots generated." | ts '[%Y-%m-%d %H:%M:%S]'

# clean the SCRATCH directory
clean_scratch

# Resources: The script ran 2,5 h, with 73 % CPU time and 100% memory used.

