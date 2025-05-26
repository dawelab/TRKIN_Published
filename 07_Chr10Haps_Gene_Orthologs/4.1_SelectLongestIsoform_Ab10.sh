module load seqtk/1.3-GCC-11.3.0

DIR="/scratch/mjb51923/TRKIN_CRISPR/out_paper/OrthoFinder"
PROT="/scratch/mjb51923/TRKIN_CRISPR/out_paper/Entap/Ab10/entap_outfiles/transcriptomes/HiFiAb10_final.fasta"

seqtk subseq $PROT $DIR/HiFiAb10.Ab10hapProtein.fasta.Ab10hapGeneNamesLongest.txt > $DIR/HiFiAb10.Ab10hapProtein.LongestIsoform.fasta

