module load seqtk/1.3-GCC-11.3.0
module load SAMtools/1.11-GCC-10.2.0

DIR=""
#This file was generated in 05.12
PROT="path/to/Entap/Ab10/entap_outfiles_v2/transcriptomes/HiFiAb10_final.fasta"

#Select the protein sequences associated with genes on the non-shared region of the Ab10 haplotype

#This loop produces a list of all the protein names
while read i
do
NAME=$i"\."
LIST=$(grep $NAME $PROT | tr ' ' '\n' | sed 's/>//g')
for j in $LIST
do
    echo $j >> $DIR/Ab10hapGeneNamesSpec.txt
done
done < $DIR/HiFiAb10.Ab10hapGeneNames.txt

#This actually selects those fasta sequences 
seqtk subseq $PROT $DIR/Ab10hapGeneNamesSpec.txt > $DIR/HiFiAb10.Ab10hapProtein.fasta
#This indexes the protein sequence fasta
samtools faidx $DIR/HiFiAb10.Ab10hapProtein.fasta


