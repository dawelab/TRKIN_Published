module load seqtk/1.3-GCC-11.3.0

DIR=""
#This file is assocaited with 
#This file is from https://www.maizegdb.org/genome/assembly/Zm-B73-REFERENCE-NAM-5.0
PROT="Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.protein.fa"

#This selects the fasta sequence for the longest isoforms identified in the previous step
seqtk subseq $PROT $DIR/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.protei.GeneNamesLongest.txt> $DIR/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.Protein.LongestIsoform.fasta

