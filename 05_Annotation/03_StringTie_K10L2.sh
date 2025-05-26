#Load the modules
module load StringTie/2.2.1-GCC-11.3.0
module load SAMtools/1.17-GCC-12.2.0

#Define the variables
OUTDIR=""
#mkdir $OUTDIR/Stringtie_K10L2
BAMDIR=$OUTDIR/K10L2_Hisat2

#Merge all the bam files 
BAM_LIST=$(ls $BAMDIR/*.s.bam)

#Merge all the bam files together
samtools merge -o $OUTDIR/Stringtie_K10L2/K10L2.AllReps.bam $BAM_LIST
#Sort the merged bam file
samtools sort $OUTDIR/Stringtie_K10L2/K10L2.AllReps.bam -o $OUTDIR/Stringtie_K10L2/K10L2.AllReps.s.bam

#Run the genome guided transcriptome assembly
stringtie -o $OUTDIR/Stringtie_K10L2/K10L2_Stringtie_Assembly.gtf $OUTDIR/Stringtie_K10L2/K10L2.AllReps.s.bam
