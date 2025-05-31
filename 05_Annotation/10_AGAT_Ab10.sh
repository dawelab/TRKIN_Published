#Load the modules 
module load AGAT/1.1.0

#Define variables
OUTDIR=""
DIR=$OUTDIR/AGAT

#Extract the CDS sequences
agat_sp_extract_sequences.pl -g $OUTDIR/Liftoff/B73_Ab10_HiFi_v2.gene.sorted.gff3 -f $OUTDIR/RepeatMasker/B73_Ab10_HiFi_v2.fa.masked -t cds -o $DIR/B73_Ab10_HiFi_v2.cds.fasta

#Extract the mRNA sequence 
agat_sp_extract_sequences.pl -g $OUTDIR/Liftoff/B73_Ab10_HiFi_v2.gene.sorted.gff3 -f $OUTDIR/RepeatMasker/B73_Ab10_HiFi_v2.fa.masked -t exon -merge -o $DIR/B73_Ab10_HiFi_v2.mRNA.fasta
