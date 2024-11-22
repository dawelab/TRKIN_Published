#!/bin/bash
#SBATCH --job-name=Trinity_K10L2
#SBATCH --output=Trinity_K10L2.out
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mjb51923@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=200gb
#SBATCH --time=24:00:00

#Load the modules
module load Trinity/2.15.1-foss-2022a

#Define the variables
OUTDIR=/scratch/mjb51923/TRKIN_CRISPR/out_paper
mkdir $OUTDIR/Trinity_Denovo_K10L2
READDIR=/scratch/mjb51923/raw_reads/RNA/K10L2

#List the fastq files 
R1_LIST=$(ls $READDIR/*R1.fq.gz)
R2_LIST=$(ls $READDIR/*R2.fq.gz)

#Concatenate the fastq reads
zcat $R1_LIST > $OUTDIR/Trinity_Denovo_K10L2/K10L2.AllReps_R1.fq
zcat $R2_LIST > $OUTDIR/Trinity_Denovo_K10L2/K10L2.AllReps_R2.fq

#Denovo assemble the K10L2 reads
for i in "K10L2.AllReps"
do
echo "########### Starting new line"
mkdir $OUTDIR/Trinity_Denovo_K10L2/trinity.$i
Trinity --seqType fq --output $OUTDIR/Trinity_Denovo_K10L2/trinity.$i --max_memory 170G --left $OUTDIR/Trinity_Denovo_K10L2/${i}_R1.fq --right $OUTDIR/Trinity_Denovo_K10L2/${i}_R2.fq --CPU 24
done