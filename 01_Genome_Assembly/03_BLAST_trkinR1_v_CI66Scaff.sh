#!/bin/bash
#SBATCH --partition=batch 
#SBATCH -J BLAST_trkin_v_CI66Scaffolds
#SBATCH --output BLAST_trkin_v_CI66Scaffolds.out
#SBATCH --mem=5GB
#SBATCH --time=2:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=30
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=END

#load the modules necessary from the cluster
module load BLAST+/2.13.0-gompi-2022a


#Define Variables, variables in all caps 
DIR="/scratch/mjb51923/CI66_Assembly/out"
NAME="CI66_rq99.asm"

#make blast database for the B73Ab10 reference genome
makeblastdb -in $DIR/${NAME}.bp.p_ctg.gfa.fasta -parse_seqids -dbtype nucl

#blast trkin against the scaffolds 
blastn -num_threads 30 -task "blastn" -outfmt 6 -db $DIR/${NAME}.bp.p_ctg.gfa.fasta -query /scratch/mjb51923/annotations/trkin_cds_1.fa -out $DIR/BLAST_trkin_v_CI66Scaff.out

#Blast colored1 against the CI66 scaffolds
blastn -num_threads 30 -task "blastn" -outfmt 6 -db $DIR/${NAME}.bp.p_ctg.gfa.fasta -query /scratch/mjb51923/annotations/Zm00001eb429330.fasta -out $DIR/BLAST_r1_v_CI66Scaff.out
