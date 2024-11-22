#!/bin/bash
#SBATCH --job-name=Mummer_HiFiAb10vHiFiAb10
#SBATCH --output=Mummer_HiFiAb10vHiFiAb10.out
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mjb51923@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=100gb
#SBATCH --time=72:00:00

#Define Variables
DIR="/scratch/mjb51923/TRKIN_CRISPR/out_paper/K10L2_Mummer"

#Load Modules
module load MUMmer/4.0.0rc1-GCCcore-11.3.0

#-t is threads, -l is minimum length generate maximal unique matches, -maxmatch has it computer all maximal matches not just unique ones, -b computes both forward and reverse complements, I believe that --nosimplify ignore all hits with the exact same start coordinate in both sequences  (mums)
cd $DIR
nucmer --nosimplify --maxmatch -p Ab10_v_Ab10_nucmer -t 30 -l 300 $DIR/Ab10_HiFi_v2_corrected_Ab10Hap.Ab10.fasta $DIR/Ab10_HiFi_v2_corrected_Ab10Hap.Ab10.fasta

#This outputs the cords
show-coords -c -l -r -T Ab10_v_Ab10_nucmer.delta > Ab10_v_Ab10_nucmer.coords