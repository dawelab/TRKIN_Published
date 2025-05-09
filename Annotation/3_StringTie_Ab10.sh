#!/bin/bash
#SBATCH --job-name=Stringtie_Ab10
#SBATCH --output=Stringtie_Ab10.out
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mjb51923@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=200gb
#SBATCH --time=24:00:00

#Load the modules
module load StringTie/2.2.1-GCC-11.3.0
module load SAMtools/1.17-GCC-12.2.0

#Define the variables
OUTDIR=/scratch/mjb51923/TRKIN_CRISPR/out_paper
#mkdir $OUTDIR/Stringtie_Ab10
BAMDIR=/scratch/mjb51923/TRKIN_CRISPR/out_paper/Hisat2

#Merge all the bam files 
BAM_LIST=$(ls $BAMDIR/*.s.bam)

#Merge all the bam files together
samtools merge -o $OUTDIR/Stringtie_Ab10/B73Ab10_AllTissues.bam $BAM_LIST
samtools sort $OUTDIR/Stringtie_Ab10/B73Ab10_AllTissues.bam -o $OUTDIR/Stringtie_Ab10/B73Ab10_AllTissues.s.bam

#Run the genome guided assembly
stringtie -o $OUTDIR/Stringtie_Ab10/B73Ab10_Stringtie_Assembly.gtf $OUTDIR/Stringtie_Ab10/B73Ab10_AllTissues.s.bam
