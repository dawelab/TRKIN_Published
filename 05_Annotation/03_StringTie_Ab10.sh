#Load the modules
module load StringTie/2.2.1-GCC-11.3.0
module load SAMtools/1.17-GCC-12.2.0

#Define the variables
OUTDIR=""
#mkdir $OUTDIR/Stringtie_Ab10
BAMDIR=$OUTDIR/Hisat2

#Merge all the bam files 
BAM_LIST=$(ls $BAMDIR/*.s.bam)

#Merge all the bam files together
samtools merge -o $OUTDIR/Stringtie_Ab10/B73Ab10_AllTissues.bam $BAM_LIST
#Sort the merged bam file
samtools sort $OUTDIR/Stringtie_Ab10/B73Ab10_AllTissues.bam -o $OUTDIR/Stringtie_Ab10/B73Ab10_AllTissues.s.bam

#Run the genome guided transcriptome assembly
stringtie -o $OUTDIR/Stringtie_Ab10/B73Ab10_Stringtie_Assembly.gtf $OUTDIR/Stringtie_Ab10/B73Ab10_AllTissues.s.bam
