#!/bin/bash
#SBATCH --partition=batch
#SBATCH -J Liftoff_HiFiAb10
#SBATCH --output Liftoff_HiFiAb10.out
#SBATCH --time=12:00:00
#SBATCH --mem=100GB
#SBATCH	--nodes=1
#SBATCH	--ntasks=12
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=END

#Load modules 
module load Liftoff/1.6.3
module load SAMtools/1.18-GCC-12.3.0
module load BEDTools/2.31.0-GCC-12.3.0

#Define Variables
DIR="/scratch/mjb51923/TRKIN_CRISPR/out_paper/Liftoff"
cd $DIR
#This subsets the K10L2 and Ab10 genomes to their respective trkin specific regions
#samtools faidx /scratch/mjb51923/ref_genomes/Ab10_HiFi_v2_corrected.fa chr10:142472000-153145000 > $DIR/Ab10_HiFi_v2_corrected.trkinreg.fa
#samtools faidx /scratch/mjb51923/TRKIN_CRISPR/out_paper/RepeatMasker/CI66_K10L2_v1.fasta K10L2:7429619-25502283 > $DIR/CI66_K10L2_v1.trkinreg.fa

#I manually edited the K10L2 trkin region file to name the region K10L2 instead of K10L2:7429619-25502283 and the same for Ab10

Ab10=/scratch/mjb51923/ref_genomes/Ab10_HiFi_v2_corrected.fa
K10L2=/scratch/mjb51923/TRKIN_CRISPR/out_paper/RepeatMasker/CI66_K10L2_v1.fasta

#Lift over the annotations from the original Ab10 genome to the HiFi genome 
#liftoff -a 0.4 -s 0.4 -copies -sc 0.7 -flank 0.1 -g $DIR/pasa_db_K10L2.gene_structures_post_PASA_updates.3815155.trkinreg.gff3 -o $DIR/Ab10_HiFi_v2_corrected.trkinreg.K10L2liftoff.gff3 $Ab10 $K10L2

#Identify the liftedoff genes that didn't appear in the original annotations
bedtools intersect -a $DIR/Ab10_HiFi_v2_corrected.trkinreg.K10L2liftoff.gff3 -b $DIR/pasa_db_Ab10.gene_structures_post_PASA_updates.280839.trkinreg.gff3 -v > $DIR/Ab10_HiFi_v2_corrected.trkinreg.K10L2liftoff.newonly.gff3

#I downloaded this file in R and subset it to only the trkin specific regions
#This combines both gff3 files to create a new gff3 file with the new genes
cat $DIR/Ab10_HiFi_v2_corrected.trkinreg.K10L2liftoff.newonly.fix.gff3  /scratch/mjb51923/TRKIN_CRISPR/out_paper/Pasa_Ab10/pasa_db_Ab10.gene_structures_post_PASA_updates.280839.gff3 > $DIR/Ab10_HiFi_v2_corrected.gene.v2.gff3
#This sorts the new gff3 file
bedtools sort -i $DIR/Ab10_HiFi_v2_corrected.gene.v2.gff3 > $DIR/Ab10_HiFi_v2_corrected.gene.v2.sorted.gff3