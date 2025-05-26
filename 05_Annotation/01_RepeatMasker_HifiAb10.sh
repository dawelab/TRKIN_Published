#!/bin/bash
#SBATCH --job-name=RepeatMasker_HifiAb10
#SBATCH --output RepeatMasker_HifiAb10.out
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mjb51923@uga.edu
#SBATCH --ntasks=40
#SBATCH --mem=200gb
#SBATCH --time=24:00:00

#Load Modules needed 
module load RepeatMasker/4.1.5-foss-2022a

#Set the file name 
REF=/scratch/mjb51923/ref_genomes/Ab10_HiFi_v2_corrected.fa
OUT=/scratch/mjb51923/TRKIN_CRISPR/out_paper

#Make the output folder
#mkdir $OUT/RepeatMasker
#cd $OUT/RepeatMasker

#Download the curated maize TE library
#git clone https://github.com/oushujun/MTEC.git
#mv MTEC/maizeTE02052020 ./maizeTE02052020.fasta
#rm -r MTEC


#Mask the repeats using both of these repeat repositories. I am choosing to use RepBase as well because Ab10 and K10L2 are 
#weird and may not be fully represented by other maize lines. I don't think it's really necessary to use repeat modeler
#given how welll curated the maize library is
#This runs repeat asker as maize on the reference using blast (-e), 40 threads (-pa), the custom librarby made above (-lib),
#outputs into the designated output directory (-dir), soft masks (-small), and a TE gff file (-gff)
RepeatMasker -e ncbi -pa 40 -lib $OUT/RepeatMasker/maizeTE02052020.fasta -dir $OUT/RepeatMasker -small -xsmall -gff $REF


