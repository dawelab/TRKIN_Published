#!/bin/bash
#SBATCH --job-name=AGAT_Ab10
#SBATCH --output=AGAT_Ab10.out
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mjb51923@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50gb
#SBATCH --time=10:00:00

#Load the modules 
module load AGAT/1.1.0

#Define variables
DIR=/scratch/mjb51923/TRKIN_CRISPR/out_paper/AGAT

#Extract the CDS sequences
agat_sp_extract_sequences.pl -g /scratch/mjb51923/TRKIN_CRISPR/out_paper/Liftoff/Ab10_HiFi_v2_corrected.gene.v2.sorted.gff3 -f /scratch/mjb51923/TRKIN_CRISPR/out_paper/RepeatMasker/Ab10_HiFi_v2_corrected.fa.masked -t cds -o $DIR/HiFiAb10.v2.cds.fasta

#Extract the mRNA sequence 
agat_sp_extract_sequences.pl -g /scratch/mjb51923/TRKIN_CRISPR/out_paper/Liftoff/Ab10_HiFi_v2_corrected.gene.v2.sorted.gff3 -f /scratch/mjb51923/TRKIN_CRISPR/out_paper/RepeatMasker/Ab10_HiFi_v2_corrected.fa.masked -t exon -merge -o $DIR/HiFiAb10.v2.mRNA.fasta
