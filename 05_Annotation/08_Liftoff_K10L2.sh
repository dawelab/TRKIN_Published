#Load modules 
module load Liftoff/1.6.3
module load SAMtools/1.18-GCC-12.3.0

#Define Variables
OUTDIR=""
DIR="$OUTDIR/Liftoff"
cd $DIR

#Define reference genomes
K10L2=/path/to/CI66_K10L2_v1.fasta
Ab10=/path/to/B73_Ab10_HiFi_v2.fa

#This subsets the K10L2 genome to the trkin specific regions (The region between the two TR-1 knobs)
samtools faidx $K10L2 K10L2:7429619-25502283 > $DIR/CI66_K10L2_v1.trkinreg.fa

#I manually edited the K10L2 trkin region file to name the region K10L2 instead of K10L2:7429619-25502283

#This isolates the trkin region of each gene annotation file
bedtools intersect -a $OUTDIR/Pasa_K10L2/PASA_updates.2.gff3 -b DIR/CI66_K10L2_v1.trkinreg.fa -o $DIR/Pasa_K10L2/PASA_updates.3.trkinreg.gff3

#Lift over the annotations from the original Ab10 genome to the HiFi genome 
liftoff -a 0.4 -s 0.4 -copies -sc 0.7 -flank 0.1 -g $OUTDIR/Pasa_Ab10/PASA_updates.3.trkinreg.gff3 -o $DIR/CI66_K10L2_v1.trkinreg.Ab10Liftoff.gff3 $K10L2 $Ab10 

#Identify the liftedoff genes that didn't appear in the original annotations
bedtools intersect -a $DIR/CI66_K10L2_v1.trkinreg.Ab10Liftoff.gff3 -b OUTDIR/Pasa_K10L2/PASA_updates.2.gff3 -v > $DIR/CI66_K10L2_v1.trkinreg.Ab10liftoff.newonly.gff3

#I downloaded this file in R and subset it to only the trkin specific regions. This script is under 08_FilterAndRename_Liftoff.R

#This combines both gff3 files to create a new gff3 file with the new genes
cat $DIR/CI66_K10L2_v1.trkinreg.Ab10liftoff.newonly.fix.gff3 $OUTDIR/Pasa_K10L2/PASA_updates.2.gff3 > $DIR/CI66_K10L2_v1.gene.gff3
#This sorts the new gff3 file
bedtools sort -i $DIR/CI66_K10L2_v1.gene.gff3 > $DIR/CI66_K10L2_v1.gene.sorted.gff3
