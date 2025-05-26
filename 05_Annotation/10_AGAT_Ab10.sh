#Load the modules 
module load AGAT/1.1.0

#Define variables
OUTDIR=""
DIR=$OUTDIR/AGAT

#Extract the CDS sequences
agat_sp_extract_sequences.pl -g $OUTDIR/Liftoff/Ab10_HiFi_v2_corrected.gene.v2.sorted.gff3 -f $OUTDIR/RepeatMasker/Ab10_HiFi_v2_corrected.fa.masked -t cds -o $DIR/HiFiAb10.v2.cds.fasta

#Extract the mRNA sequence 
agat_sp_extract_sequences.pl -g $OUTDIR/Liftoff/Ab10_HiFi_v2_corrected.gene.v2.sorted.gff3 -f $OUTDIR/RepeatMasker/Ab10_HiFi_v2_corrected.fa.masked -t exon -merge -o $DIR/HiFiAb10.v2.mRNA.fasta
