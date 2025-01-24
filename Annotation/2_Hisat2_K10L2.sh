#!/bin/bash
#SBATCH --job-name=Hisat2_K10L2
#SBATCH --output Hisat2_K10L2.%A-%a.out
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mjb51923@uga.edu
#SBATCH --ntasks=24
#SBATCH --mem=100gb
#SBATCH --time=24:00:00
#SBATCH --array=1-5

#Load Modules needed 
module load HISAT2/3n-20201216-gompi-2022a
module load SAMtools/1.17-GCC-12.2.0

#Set the file name 
THREADS=24
READDIR=/scratch/mjb51923/TRKIN_CRISPR/out_paper/Trimmed_RNAseq_K10L2
OUT=/scratch/mjb51923/TRKIN_CRISPR/out_paper

#Make the output directory and enter it
#mkdir $OUT/Hisat2
cd $OUT/K10L2_Hisat2

#List of all RNA seq identifiers (one per pair) from EBI E-MTAB-8641 made in List.txt

#This pulls info from the array job
IT=$SLURM_ARRAY_TASK_ID
NAME=$(awk NR==${IT}'{print $1}' $OUT/K10L2_Hisat2/List.txt)

#Define the masked reference 
REF=$OUT/RepeatMasker/CI66_K10L2_v1.fasta.masked

#Build a hisat2 reference
#I ran this command separatly first to avoid having the reference re built every time
#hisat2-build $REF CI66_K10L2_v1.fasta.masked
hisat2 -x $OUT/K10L2_Hisat2/CI66_K10L2_v1.fasta.masked -p $THREADS -1 $READDIR/${NAME}_R1_paired.fq -2 $READDIR/${NAME}_R2_paired.fq -S $OUT/K10L2_Hisat2/${NAME}.sam
samtools view -bS $OUT/K10L2_Hisat2/${NAME}.sam > $OUT/K10L2_Hisat2/${NAME}.bam
samtools sort $OUT/K10L2_Hisat2/${NAME}.bam -o $OUT/K10L2_Hisat2/${NAME}.s.bam


