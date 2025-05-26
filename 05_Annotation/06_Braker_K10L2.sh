#Load the modules
module load BRAKER/3.0.8-foss-2022a
module load GeneMark-ETP/1.0.0-GCCcore-11.3.0
module load BEDTools/2.30.0-GCC-12.2.0
module load gffread/0.12.7-GCCcore-11.3.0
module load StringTie/2.2.1-GCC-11.3.0

#Define the variables
OUT=""
REF=$OUT/RepeatMasker/CI66_K10L2_v1.fasta.masked
BAM=$OUT/K10L2_Hisat2

#Make the output directory 
cd $OUT
#mkdir $OUT/K10L2_Braker
cd $OUT/K10L2_Braker

#List the sorted bam files from the Hisat2 directory
LIST=$(ls $BAM/*.s.bam | tr '\n' ,)

#Download the orthodb protein file for plants

wget https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/Viridiplantae.fa.gz
gunzip Viridiplantae.fa.gz

#Run Braker to annotate the genome
braker.pl \
--genome $REF \
--bam $LIST \
--prot_seq $OUT/K10L2_Braker/Viridiplantae.fa \
--threads 40 \
--workingdir $OUT/K10L2_Braker/ \
--nocleanup



