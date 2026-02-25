### Script for Metacentrum

#PBS -N Download_from_ENA
#PBS -l select=1:ncpus=1:mem=8gb:scratch_local=10gb:brno=True
#PBS -l walltime=24:00:00 
#PBS -m ae
# I ask for machine in Brno, because it is close to storage and the data transfer might be faster.

cd /storage/brno12-cerit/home/duchmil/annotations/Noccaea_praecox_2024_12/rnaseq/1_raw_reads/Bocaj_et_al_2023

echo "Download start"  | ts '[%Y-%m-%d %H:%M:%S]'

wget -nv ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR249/066/SRR24947466/SRR24947466_1.fastq.gz
wget -nv ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR249/061/SRR24947461/SRR24947461_2.fastq.gz
wget -nv ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR249/061/SRR24947461/SRR24947461_1.fastq.gz
wget -nv ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR249/062/SRR24947462/SRR24947462_1.fastq.gz
wget -nv ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR249/065/SRR24947465/SRR24947465_1.fastq.gz
wget -nv ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR249/063/SRR24947463/SRR24947463_1.fastq.gz
wget -nv ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR249/064/SRR24947464/SRR24947464_1.fastq.gz
wget -nv ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR249/062/SRR24947462/SRR24947462_2.fastq.gz
wget -nv ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR249/066/SRR24947466/SRR24947466_2.fastq.gz
wget -nv ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR249/063/SRR24947463/SRR24947463_2.fastq.gz
wget -nv ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR249/064/SRR24947464/SRR24947464_2.fastq.gz
wget -nv ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR249/065/SRR24947465/SRR24947465_2.fastq.gz

echo "Download done"  | ts '[%Y-%m-%d %H:%M:%S]'

md5sum * > md5_to_check.txt

echo "md5sum done"  | ts '[%Y-%m-%d %H:%M:%S]'


# clean the SCRATCH directory
clean_scratch