module load seqtk/1.3-GCC-11.3.0

DIR=""
#This file was generated in 05.12
PROT="/path/to/Entap/Ab10/entap_outfiles/transcriptomes/HiFiAb10_final.fasta"

#This selects the fasta sequence for the longest isoforms identified in the previous step
seqtk subseq $PROT $DIR/HiFiAb10.Ab10hapProtein.fasta.Ab10hapGeneNamesLongest.txt > $DIR/HiFiAb10.Ab10hapProtein.LongestIsoform.fasta

