#!/bin/bash
#SBATCH --job-name=Mummer_N10vK10L2
#SBATCH --output=Mummer_N10vK10L2.out
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mjb51923@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=100gb
#SBATCH --time=72:00:00

#Load Modules
module load MUMmer/4.0.0rc1-GCCcore-11.3.0
module load SAMtools/1.10-GCC-8.3.0

#Define Variables
K10L2="/scratch/mjb51923/TRKIN_CRISPR/out_paper/RepeatMasker/CI66_K10L2_v1.fasta"
AB10='/scratch/mjb51923/ref_genomes/Ab10_HiFi_v2_corrected.fa'
N10='/scratch/mjb51923/ref_genomes/B73.PLATINUM.pseudomolecules-v1.fasta'
DIR="/scratch/mjb51923/TRKIN_CRISPR/out_paper/K10L2_Mummer"

#Make Region files 
printf "K10L2:2730186-31891546" > $DIR/CI66_K10L2.bed
printf "chr10:141115174-195055488" > $DIR/Ab10_HiFi_v2_corrected_Ab10Hap.bed
printf "chr10:141187279-152435371" > $DIR/N10_Ab10SharedRegion.bed

#Subset K10L2 reference to just chr 10 
samtools faidx $K10L2 -r $DIR/CI66_K10L2.bed -o $DIR/CI66_K10L2_v1.K10L2.fasta
samtools faidx $AB10 -r $DIR/Ab10_HiFi_v2_corrected_Ab10Hap.bed -o $DIR/Ab10_HiFi_v2_corrected_Ab10Hap.Ab10.fasta
samtools faidx $N10 -r $DIR/N10_Ab10SharedRegion.bed -o $DIR/B73.PLATINUM.pseudomolecules-v1.N10SharedRegion.fasta

#-t is threads, -l is minimum length generate maximal unique matches, -maxmatch has it compute all maximal matches not just unique ones, -b computes both forward and reverse complements, I believe that
cd $DIR
nucmer --maxmatch -p N10_v_K10L2_nucmer -t 30 -l 300 $DIR/B73.PLATINUM.pseudomolecules-v1.N10SharedRegion.fasta $DIR/CI66_K10L2_v1.K10L2.fasta

#This outputs the cords
show-coords -c -l -r -T N10_v_K10L2_nucmer.delta > N10_v_K10L2_nucmer.coords