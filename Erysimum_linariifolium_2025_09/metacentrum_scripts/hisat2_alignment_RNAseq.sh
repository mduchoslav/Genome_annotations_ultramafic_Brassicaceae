### Script for Metacentrum

#PBS -N hisat2_alignment_RNAseq
#PBS -l select=1:ncpus=16:mem=128gb:scratch_local=1000gb
#PBS -l walltime=24:00:00 
#PBS -m ae


# define a DATADIR variable: directory where the input files are taken from and where output will be copied to
DATADIR=/storage/brno12-cerit/home/duchmil/annotations/Erysimum_linariifolium_2025_09
# assembly name
genome_assembly="Erysimum_linariifolium_CUNI_V1_2024_09_masked.fa"

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

