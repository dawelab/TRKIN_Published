#!/bin/bash
#SBATCH --job-name=Mummer_OGAb10vHiFiAb10
#SBATCH --output=Mummer_OGAb10vHiFiAb10.out
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
DIR="/scratch/mjb51923/TRKIN_CRISPR/out_paper/K10L2_Mummer"

#This subsets the Ab10 original assembly
printf "chr10:140964412-195026473" > $DIR/B73Ab10_Ab10Hap.bed
samtools faidx  /scratch/mjb51923/ref_genomes/Zm-B73_AB10-REFERENCE-NAM-1.0.fa -r $DIR/B73Ab10_Ab10Hap.bed > $DIR/Zm-B73_AB10-REFERENCE-NAM-1.0.Ab10Hap.fa

#-t is threads, -l is minimum length generate maximal unique matches, -maxmatch has it compute all maximal matches not just unique ones, -b computes both forward and reverse complements, I believe that
cd $DIR
nucmer --maxmatch -p OGAb10_v_HiFiAb10_nucmer -t 30 -l 300 $DIR/Zm-B73_AB10-REFERENCE-NAM-1.0.Ab10Hap.fa $DIR/Ab10_HiFi_v2_corrected_Ab10Hap.Ab10.fasta

#This outputs the cords
show-coords -c -l -r -T OGAb10_v_HiFiAb10_nucmer.delta > OGAb10_v_HiFiAb10_nucmer.coords
