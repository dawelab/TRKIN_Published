#Load Modules needed 
module load HISAT2/3n-20201216-gompi-2022a
module load SAMtools/1.17-GCC-12.2.0

#Set the file name 
THREADS=24
#Ab10 RNA seq reads come from https://doi.org/10.1186/s13059-020-02029-9
READDIR=path/to/RNAseq_Ab10
OUT=""

#Make the output directory and enter it
#mkdir $OUT/Hisat2
cd $OUT/Hisat2

#Manuallt list of all RNA seq identifiers (one per pair) from Bioproject PRJEB35367 made in List.txt

#This pulls info from the array job
IT=$SLURM_ARRAY_TASK_ID
NAME=$(awk NR==${IT}'{print $1}' $OUT/Hisat2/List.txt)

#Define the masked reference 
REF=$OUT/RepeatMasker/Ab10_HiFi_v2_corrected.fa.masked

#Build a hisat2 reference
#I ran this command separatly first to avoid having the reference re built every time
#hisat2-build $REF Ab10_HiFi_v2_corrected.fa.masked
#Align the RNAseq reads to the masked reference
hisat2 -x $OUT/Hisat2/Ab10_HiFi_v2_corrected.fa.masked -p $THREADS -1 $READDIR/${NAME}_R1_paired.fq -2 $READDIR/${NAME}_R2_paired.fq -S $OUT/Hisat2/${NAME}.sam
#Cenvert the sam file output by hisat2 to a bam file
samtools view -bS $OUT/Hisat2/${NAME}.sam > $OUT/Hisat2/${NAME}.bam
#Sorted the bam file
samtools sort $OUT/Hisat2/${NAME}.bam -o $OUT/Hisat2/${NAME}.s.bam


