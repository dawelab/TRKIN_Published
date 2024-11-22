#!/bin/bash
#SBATCH --job-name=Braker_HifiAb10
#SBATCH --output=Braker_HifiAb10.out
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mjb51923@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=200gb
#SBATCH --time=72:00:00

#Load the modules
module load BRAKER/3.0.8-foss-2022a
module load GeneMark-ETP/1.0.0-GCCcore-11.3.0
module load BEDTools/2.30.0-GCC-12.2.0
module load gffread/0.12.7-GCCcore-11.3.0
module load StringTie/2.2.1-GCC-11.3.0

#Define the variables
OUT=/scratch/mjb51923/TRKIN_CRISPR/out_paper
REF=$OUT/RepeatMasker/Ab10_HiFi_v2_corrected.fa.masked
BAM=$OUT/Hisat2

#Make the output directory 
#cd $OUT
#mkdir $OUT/Braker

#List the sorted bam files from the Hisat2 directory
LIST=$(ls $BAM/*.s.bam | tr '\n' ,)

#Download the orthodb protein file for plants
#wget https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/Viridiplantae.fa.gz
#gunzip Viridiplantae.fa.gz

#Run Braker. I can't add UTRs here for some reason. I can do it later
braker.pl \
--genome $REF \
--bam $LIST \
--prot_seq $OUT/Braker/Viridiplantae.fa \
--threads 40 \
--workingdir $OUT/Braker/ \
--nocleanup



