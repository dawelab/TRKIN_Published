module load seqtk/1.3-GCC-11.3.0

DIR="/scratch/mjb51923/TRKIN_CRISPR/out_paper/OrthoFinder"
PROT="/scratch/mjb51923/TRKIN_CRISPR/out_paper/Entap/K10L2/entap_outfiles_nr_update/transcriptomes/CI66_K10L2_final.fasta"

seqtk subseq $PROT $DIR/CI66_K10L2.K10L2hapProtein.fasta.K10L2hapGeneNamesLongest.txt > $DIR/CI66_K10L2.K10L2hapProtein.LongestIsoform.fasta

