module load seqtk/1.3-GCC-11.3.0
module load SAMtools/1.11-GCC-10.2.0

DIR=""
#This file was generated in 05.12
PROT="/path/to/Entap/K10L2/entap_outfiles_nr_update/transcriptomes/CI66_K10L2_final.fasta"

#Select the protein sequences associated with genes on the non-shared region of the K10L2 haplotype

#This loop produces a list of all the protein names
while read i;
do
NAME=$i"\."
LIST=$(grep $NAME $PROT | tr ' ' '\n' | sed 's/>//g')
for j in $LIST;
do
    echo $j >> $DIR/K10L2GeneNamesSpec.txt
done
done < $DIR/CI66_K10L2.K10L2hapGeneNames.txt

#This actually selects those fasta sequences 
seqtk subseq $PROT $DIR/K10L2GeneNamesSpec.txt > $DIR/CI66_K10L2.K10L2hapProtein.fasta
#This indexes the protein sequence fasta
samtools faidx $DIR/CI66_K10L2.K10L2hapProtein.fasta


