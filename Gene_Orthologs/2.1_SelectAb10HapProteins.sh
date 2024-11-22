#!/bin/bash
#SBATCH --partition=batch
#SBATCH -J Seqkit_Ab10HapPro
#SBATCH --output Seqkit_Ab10HapPro.out
#SBATCH --mem=50000
#SBATCH --time=80:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=1
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=END

module load seqtk/1.3-GCC-11.3.0
module load SAMtools/1.11-GCC-10.2.0

DIR="/scratch/mjb51923/TRKIN_CRISPR/out_paper/OrthoFinder"
PROT="/scratch/mjb51923/TRKIN_CRISPR/out_paper/Entap/Ab10/entap_outfiles_v2/transcriptomes/HiFiAb10_final.fasta"

while read i
do
NAME=$i"\."
LIST=$(grep $NAME $PROT | tr ' ' '\n' | sed 's/>//g')
for j in $LIST
do
    echo $j >> $DIR/Ab10hapGeneNamesSpec.txt
done
done < $DIR/HiFiAb10.Ab10hapGeneNames.txt

seqtk subseq $PROT $DIR/Ab10hapGeneNamesSpec.txt > $DIR/HiFiAb10.Ab10hapProtein.fasta
samtools faidx $DIR/HiFiAb10.Ab10hapProtein.fasta


