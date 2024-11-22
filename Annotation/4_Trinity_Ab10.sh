#!/bin/bash
#SBATCH --job-name=Trinity_Ab10
#SBATCH --output=Trinity_Ab10.out
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mjb51923@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=300gb
#SBATCH --time=24:00:00

#Load the modules
module load Trinity/2.15.1-foss-2022a

#Define the variables
OUTDIR=/scratch/mjb51923/TRKIN_CRISPR/out_paper
#mkdir $OUTDIR/Trinity_Denovo_Ab10
READDIR=/scratch/mjb51923/raw_reads/RNA/Gapless_B73-Ab10I

#List the fastq files 
R1_LIST=$(ls $READDIR/*R1.fq.gz)
R2_LIST=$(ls $READDIR/*R2.fq.gz)

#Concatenate the fastq reads
#zcat $R1_LIST > $OUTDIR/Trinity_Denovo_Ab10/B73Ab10.AllTissues_R1.fq
#zcat $R2_LIST > $OUTDIR/Trinity_Denovo_Ab10/B73Ab10.AllTissues_R2.fq

#Denovo assemble the Ab10 reads
for i in "B73Ab10.AllTissues"
do
echo "########### Starting new line"
mkdir $OUTDIR/Trinity_Denovo_Ab10/trinity.$i
Trinity --seqType fq --output $OUTDIR/Trinity_Denovo_Ab10/trinity.$i --max_memory 270G --left $OUTDIR/Trinity_Denovo_Ab10/${i}_R1.fq --right $OUTDIR/Trinity_Denovo_Ab10/${i}_R2.fq --CPU 24
done