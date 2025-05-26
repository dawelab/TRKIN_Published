#Load the modules 
module load AGAT/1.1.0

#Define variables
OUTDIR=""
DIR=$OUTDIR/AGAT

#Extract the CDS sequences
agat_sp_extract_sequences.pl -g $OUTDIR/Liftoff/CI66_K10L2_v1.gene.v2.sorted.gff3 -f  $OUTDIR/RepeatMasker/CI66_K10L2_v1.fasta.masked -t cds -o $DIR/CI66_K10L2.v2.cds.fasta

#Extract the mRNA sequences
agat_sp_extract_sequences.pl -g $OUTDIR/Liftoff/CI66_K10L2_v1.gene.v2.sorted.gff3 -f  $OUTDIR/RepeatMasker/CI66_K10L2_v1.fasta.masked -t exon --merge -o $DIR/CI66_K10L2.v2.mRNA.fasta
