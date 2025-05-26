#!/bin/bash
#SBATCH --job-name=BLAST_CI66_K10L2
#SBATCH --output=BLAST_CI66_K10L2.out
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mjb51923@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=5gb
#SBATCH --time=5:00:00

#load the modules necessary from the cluster
module load BLAST+/2.13.0-gompi-2022a

#Define Variables, variables in all caps 
DIR="/scratch/mjb51923/TRKIN_CRISPR/out_paper"
REF=$DIR/RepeatMasker/CI66_K10L2_v1.fasta.masked

#make blast database
makeblastdb -in $REF -parse_seqids -dbtype nucl

#Blast various relevant features to the CI66_K10L2 genome and convert them to a bed for IGV visualization
blastn -num_threads 30 -task "blastn" -outfmt 6 -db $REF -query /scratch/mjb51923/annotations/TR1.fasta -out $DIR/BLAST/BLAST_TR1_v_CI66_K10L2.out
awk '{print $2, $9, $10, $3}' $DIR/BLAST/BLAST_TR1_v_CI66_K10L2.out > $DIR/BLAST/BLAST_TR1_v_CI66_K10L2.bed

blastn -num_threads 30 -task "blastn" -outfmt 6 -db $REF -query /scratch/mjb51923/annotations/knob180.fasta -out $DIR/BLAST/BLAST_knob180_v_CI66_K10L2.out
awk '{print $2, $9, $10, $3}' $DIR/BLAST/BLAST_knob180_v_CI66_K10L2.out > $DIR/BLAST/BLAST_knob180_v_CI66_K10L2.bed

blastn -num_threads 30 -task "blastn" -outfmt 6 -db $REF -query /scratch/mjb51923/annotations/trkin_cds_1.fa -out $DIR/BLAST/BLAST_trkin_v_CI66_K10L2.out
awk '{print $2, $9, $10, $3}' $DIR/BLAST/BLAST_trkin_v_CI66_K10L2.out > $DIR/BLAST/BLAST_trkin_v_CI66_K10L2.bed

blastn -num_threads 30 -task "blastn" -outfmt 6 -db $REF -query /scratch/mjb51923/annotations/Zm00001eb434490_T001.fasta -out $DIR/BLAST/BLAST_sr2_v_CI66_K10L2.out
awk '{print $2, $9, $10, $3}' $DIR/BLAST/BLAST_sr2_v_CI66_K10L2.out > $DIR/BLAST/BLAST_sr2_v_CI66_K10L2.bed

blastn -num_threads 30 -task "blastn" -outfmt 6 -db $REF -query /scratch/mjb51923/annotations/Zm00001eb429330.fasta -out $DIR/BLAST/BLAST_r1_v_CI66_K10L2.out
awk '{print $2, $9, $10, $3}' $DIR/BLAST/BLAST_r1_v_CI66_K10L2.out > $DIR/BLAST/BLAST_r1_v_CI66_K10L2.bed
