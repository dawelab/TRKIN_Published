#!/bin/bash
#SBATCH --job-name=OrthoFinder_Ab10_K10L2
#SBATCH --output OrthoFinder_Ab10_K10L2.out
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mjb51923@uga.edu
#SBATCH --ntasks=40
#SBATCH --mem=200gb
#SBATCH --time=24:00:00

#load the modules necessary from the cluster
module load OrthoFinder/2.5.5-foss-2023a

#Define Variables 
DIR="/scratch/mjb51923/TRKIN_CRISPR/out_paper/OrthoFinder"
ODIRNAME="Ab10_K10L2_proteomes"
PROT1=$DIR/HiFiAb10.Ab10hapProtein.LongestIsoform.fasta
PROT2=$DIR/CI66_K10L2.K10L2hapProtein.LongestIsoform.fasta

#make a directory with the proteomes of the two lines
mkdir $DIR/$ODIRNAME
cp $PROT1 $DIR/$ODIRNAME
cp $PROT2 $DIR/$ODIRNAME

#These are the files in the folder $DIR/Ab10_N10_proteomes: Zm-B73_AB10-REFERENCE-NAM-1.0_Zm00043a.1.evd.protein.Ab10HapLongest.fa Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.protein.Longest.fa
orthofinder -t 40 -f $DIR/$ODIRNAME

ODIRNAME="Ab10_N10_proteomes"
PROT1=$DIR/HiFiAb10.Ab10hapProtein.LongestIsoform.fasta
PROT2=$DIR/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.Protein.LongestIsoform.fasta

#make a directory with the proteomes of the two lines
mkdir $DIR/$ODIRNAME
cp $PROT1 $DIR/$ODIRNAME
cp $PROT2 $DIR/$ODIRNAME

#These are the files in the folder $DIR/Ab10_N10_proteomes: Zm-B73_AB10-REFERENCE-NAM-1.0_Zm00043a.1.evd.protein.Ab10HapLongest.fa Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.protein.Longest.fa
orthofinder -t 40 -f $DIR/$ODIRNAME

ODIRNAME="K10L2_N10_proteomes"
PROT1=$DIR/CI66_K10L2.K10L2hapProtein.LongestIsoform.fasta
PROT2=$DIR/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.Protein.LongestIsoform.fasta

#make a directory with the proteomes of the two lines
mkdir $DIR/$ODIRNAME
cp $PROT1 $DIR/$ODIRNAME
cp $PROT2 $DIR/$ODIRNAME

#These are the files in the folder $DIR/Ab10_N10_proteomes: Zm-B73_AB10-REFERENCE-NAM-1.0_Zm00043a.1.evd.protein.Ab10HapLongest.fa Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.protein.Longest.fa
orthofinder -t 40 -f $DIR/$ODIRNAME
