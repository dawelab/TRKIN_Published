module load seqtk/1.3-GCC-11.3.0

DIR=""
#This file was generated in 05.12
PROT="/path/to/Entap/K10L2/entap_outfiles_nr_update/transcriptomes/CI66_K10L2_final.fasta"

#This selects the fasta sequence for the longest isoforms identified in the previous step
seqtk subseq $PROT $DIR/CI66_K10L2.K10L2hapProtein.fasta.K10L2hapGeneNamesLongest.txt > $DIR/CI66_K10L2.K10L2hapProtein.LongestIsoform.fasta

