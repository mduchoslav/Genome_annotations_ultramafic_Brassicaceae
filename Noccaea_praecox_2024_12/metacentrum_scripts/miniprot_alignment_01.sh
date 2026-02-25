### Script for Metacentrum

#PBS -N miniprot_alignment_01
#PBS -l select=1:ncpus=16:mem=16gb:scratch_local=1000gb
#PBS -l walltime=2:00:00 
#PBS -m ae

# define a DATADIR variable: directory where the input files are taken from and where output will be copied to
DATADIR=/storage/brno12-cerit/home/duchmil/annotations/Noccaea_praecox_2024_12

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

cp -r $DATADIR/genome_assembly/Noccaea_polished_2024_10_Mahnaz.masked.fa $SCRATCHDIR || { echo >&2 "Error while copying input file(s)!"; exit 2; }

echo "Input files copied." | ts '[%Y-%m-%d %H:%M:%S]'


# move into scratch directory
cd $SCRATCHDIR 

# make index for miniprot
$MINIPROTDIR/miniprot -t16 -d miniprot_index.mpi Noccaea_polished_2024_10_Mahnaz.masked.fa

echo "Miniprot index prepared." | ts '[%Y-%m-%d %H:%M:%S]'

# running miniprot
$MINIPROTDIR/miniprot -Iut16 --gff miniprot_index.mpi Araport11_pep_20220914.gz > A.thaliana.pep.gff
echo "A. thaliana done." | ts '[%Y-%m-%d %H:%M:%S]'
$MINIPROTDIR/miniprot -Iut16 --gff miniprot_index.mpi Brapa_chiifu_v41_gene20230413.gff3.pep.fa.gz > B.rapa.pep.gff
echo "B. rapa done." | ts '[%Y-%m-%d %H:%M:%S]'
$MINIPROTDIR/miniprot -Iut16 --gff miniprot_index.mpi GCF_000004255.2_v.1.0_protein.faa.gz > A.lyrata.pep.gff
echo "A. lyrata done." | ts '[%Y-%m-%d %H:%M:%S]'




# move the output to user's DATADIR or exit in case of failure
cp -v *.gff $DATADIR/protein_seqs_input/2_aligned | ts '[%Y-%m-%d %H:%M:%S]' || { echo >&2 "Result file(s) copying failed (with a code $?) !!"; exit 4; }

# clean the SCRATCH directory
clean_scratch