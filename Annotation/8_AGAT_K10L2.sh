#!/bin/bash
#SBATCH --job-name=AGAT_K10L2
#SBATCH --output=AGAT_K10L2.out
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
agat_sp_extract_sequences.pl -g /scratch/mjb51923/TRKIN_CRISPR/out_paper/Liftoff/CI66_K10L2_v1.gene.v2.sorted.gff3 -f /scratch/mjb51923/TRKIN_CRISPR/out_paper/RepeatMasker/CI66_K10L2_v1.fasta.masked -t cds -o $DIR/CI66_K10L2.v2.cds.fasta

#Extract the mRNA sequences
agat_sp_extract_sequences.pl -g /scratch/mjb51923/TRKIN_CRISPR/out_paper/Liftoff/CI66_K10L2_v1.gene.v2.sorted.gff3 -f /scratch/mjb51923/TRKIN_CRISPR/out_paper/RepeatMasker/CI66_K10L2_v1.fasta.masked -t exon --merge -o $DIR/CI66_K10L2.v2.mRNA.fasta
