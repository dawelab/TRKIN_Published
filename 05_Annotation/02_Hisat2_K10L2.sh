#Load Modules needed 
module load HISAT2/3n-20201216-gompi-2022a
module load SAMtools/1.17-GCC-12.2.0

#Set the file name 
THREADS=24
#K10L2 RNAseq reads come from doi: 10.1101/gad.340679.120
READDIR=path/to//RNAseq_K10L2
OUT=""

#Make the output directory and enter it
#mkdir $OUT/Hisat2
cd $OUT/K10L2_Hisat2

#List of all RNA seq identifiers (one per pair) from NCBI BioProject PRJNA285341 made in List.txt

#This pulls info from the array job
IT=$SLURM_ARRAY_TASK_ID
NAME=$(awk NR==${IT}'{print $1}' $OUT/K10L2_Hisat2/List.txt)

#Define the masked reference 
REF=$OUT/RepeatMasker/CI66_K10L2_v1.fasta.masked

#Build a hisat2 reference
#I ran this command separatly first to avoid having the reference re built every time
#hisat2-build $REF CI66_K10L2_v1.fasta.masked
#Align RNAseq reads to the masked reference
hisat2 -x $OUT/K10L2_Hisat2/CI66_K10L2_v1.fasta.masked -p $THREADS -1 $READDIR/${NAME}_R1_paired.fq -2 $READDIR/${NAME}_R2_paired.fq -S $OUT/K10L2_Hisat2/${NAME}.sam
#Convert the sam file output by hisat2 to a bam
samtools view -bS $OUT/K10L2_Hisat2/${NAME}.sam > $OUT/K10L2_Hisat2/${NAME}.bam
#sort the bam file
samtools sort $OUT/K10L2_Hisat2/${NAME}.bam -o $OUT/K10L2_Hisat2/${NAME}.s.bam


