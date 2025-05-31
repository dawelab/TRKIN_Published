#Load modules 
module load Liftoff/1.6.3
module load SAMtools/1.18-GCC-12.3.0
module load BEDTools/2.31.0-GCC-12.3.0

#Define Variables
$OUTDIR=""
DIR=$OUTDIR"/Liftoff"
cd $DIR

#Define reference genomes
Ab10=/path/to/B73_Ab10_HiFi_v2.fa
K10L2=/path/to/CI66_K10L2_v1.fasta

#This subsets the K10L2 and Ab10 genomes to their respective trkin specific regions
samtools faidx $Ab10 chr10:142472000-153145000 > $DIR/B73_Ab10_HiFi_v2.trkinreg.fa
samtools faidx $K10L2 K10L2:7429619-25502283 > $DIR/CI66_K10L2_v1.trkinreg.fa

#I manually edited the K10L2 trkin region file to name the region K10L2 instead of K10L2:7429619-25502283 and the same for Ab10

#This isolates the trkin region of each gene annotation file
bedtools intersect -a $OUTDIR/Pasa_Ab10/PASA_updates.2.gff3 -b $DIR/B73_Ab10_HiFi_v2.trkinreg.fa -o $DIR/Pasa_Ab10/PASA_updates.3.trkinreg.gff3 
bedtools intersect -a $OUTDIR/Pasa_K10L2/PASA_updates.2.gff3 -b DIR/CI66_K10L2_v1.trkinreg.fa -o $DIR/Pasa_K10L2/PASA_updates.3.trkinreg.gff3 

#Lift over the annotations from K10L2 to the new Ab10 genome 
liftoff -a 0.4 -s 0.4 -copies -sc 0.7 -flank 0.1 -g $DIR/Pasa_Ab10/PASA_updates.3.trkinreg.gff3 -o $DIR/B73_Ab10_HiFi_v2.trkinreg.K10L2liftoff.gff3 $Ab10 $K10L2

#Identify the liftedoff genes that didn't appear in the original annotations
bedtools intersect -a $DIR/B73_Ab10_HiFi_v2.trkinreg.K10L2liftoff.gff3 -b $DIR/Pasa_Ab10/PASA_updates.3.trkinreg.gff3  -v > $DIR/B73_Ab10_HiFi_v2.trkinreg.K10L2liftoff.newonly.gff3

#I downloaded this file in R and subset it to only the trkin specific regions These scripts are under 09_FilterAndRename_Liftoff.R

#This combines both gff3 files to create a new gff3 file with the new genes
cat $DIR/B73_Ab10_HiFi_v2.trkinreg.K10L2liftoff.newonly.gff3 $OUTDIR/Pasa_Ab10/PASA_updates.2.gff3 > $DIR/B73_Ab10_HiFi_v2.gene.gff3
#This sorts the new gff3 file
bedtools sort -i $$DIR/B73_Ab10_HiFi_v2.gene.gff3 > $DIR/B73_Ab10_HiFi_v2.gene.sorted.gff3
