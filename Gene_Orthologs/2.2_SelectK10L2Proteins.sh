#!/bin/bash
#SBATCH --partition=batch
#SBATCH -J Seqkit_K10L2Pro
#SBATCH --output Seqkit_K10L2Pro.out
#SBATCH --mem=50000
#SBATCH --time=80:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=1
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=END

module load seqtk/1.3-GCC-11.3.0
module load SAMtools/1.11-GCC-10.2.0

DIR="/scratch/mjb51923/TRKIN_CRISPR/out_paper/OrthoFinder"
PROT="/scratch/mjb51923/TRKIN_CRISPR/out_paper/Entap/K10L2/entap_outfiles_nr_update/transcriptomes/CI66_K10L2_final.fasta"

while read i;
do
NAME=$i"\."
LIST=$(grep $NAME $PROT | tr ' ' '\n' | sed 's/>//g')
for j in $LIST;
do
    echo $j >> $DIR/K10L2GeneNamesSpec.txt
done
done < $DIR/CI66_K10L2.K10L2hapGeneNames.txt

seqtk subseq $PROT $DIR/K10L2GeneNamesSpec.txt > $DIR/CI66_K10L2.K10L2hapProtein.fasta
samtools faidx $DIR/CI66_K10L2.K10L2hapProtein.fasta


