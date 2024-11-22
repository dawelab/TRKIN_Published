module load seqtk/1.3-GCC-11.3.0

DIR="/scratch/mjb51923/TRKIN_CRISPR/out_paper/OrthoFinder"
PROT="/scratch/mjb51923/annotations/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.protein.fa"

seqtk subseq $PROT $DIR/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.protei.GeneNamesLongest.txt> $DIR/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.Protein.LongestIsoform.fasta

