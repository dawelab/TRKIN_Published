#!/bin/bash
#SBATCH --partition=batch
#SBATCH -J BLAST_Df_Repeats_All
#SBATCH --output /scratch/mjb51923/Ab10_FT_Mapping/scripts/BLAST_Df_Repeats_All.out
#SBATCH --mem=100000
#SBATCH --time=100:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks-per-node=30
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=END

#load the modules necessary from the cluster
module load BLAST/2.2.26-Linux_x86_64
module load seqtk/1.2-foss-2019b
module load BEDTools/2.30.0-GCC-8.3.0

#Define the variables
DIR="/scratch/mjb51923/Ab10_FT_Mapping/out/Df_Seq/Repeat_Abundance"
REF_DIR="/scratch/mjb51923/annotations"
READ_DIR="/scratch/mjb51923/raw_reads/Ab10"

#Index the references
#formatdb -p F -o T -i $REF_DIR/knob180.fasta
#formatdb -p F -o T -i $REF_DIR/TR1.fasta
#formatdb -p F -o T -i $REF_DIR/CentC.fasta

#Convert the Sequence to .fasta
#for x in "W23-Ab10_R1_001.fastq.gz" \
#"W23-DfL_R1_001.fastq.gz" \
#"W23-DfK_R1_001.fastq.gz"
#do
#	SAMP=$(echo $x | cut -d "_" -f 1)
#	#Convert .fastq reads to .fasta so BLAST can process them 
#	seqtk seq -a $READ_DIR/$x > $DIR/${SAMP}_R1.fasta
#done

#for x in "W23-Ab10_R2_001.fastq.gz" \
#"W23-DfL_R2_001.fastq.gz" \
#"W23-DfK_R2_001.fastq.gz"
#do
# 	SAMP=$(echo $x | cut -d "_" -f 1)
# 	#Convert .fastq reads to .fasta so BLAST can process them 
# 	seqtk seq -a $READ_DIR/$x > $DIR/${SAMP}_R2.fasta
# done

#BLAST all the reads for all of my W23 lines to each repeat. These are nested for loops, but bash will not accept them with it indented to show this
for i in "knob180.fasta" \
"TR1.fasta" \
"CentC.fasta"
do 
	REF=$(echo $i | cut -d "." -f 1)
for x in "W23-Ab10_R1_001.fastq.gz" \
"W23-DfL_R1_001.fastq.gz" \
"W23-DfK_R1_001.fastq.gz"
do
	SAMP=$(echo $x | cut -d "_" -f 1)
	#BLAST the reads
	blastall -b 5000 -F F -p blastn -d $REF_DIR/$i -i $DIR/${SAMP}_R1.fasta -o $DIR/"BLAST_"$REF"_"$SAMP"_R1.2.out" -m 8
done
done

#BLAST all the reads for all of my W23 lines to each repeat. These are nested for loops, but bash will not accept them with it indented to show this
for i in "knob180.fasta" \
"TR1.fasta" \
"CentC.fasta"
do 
	REF=$(echo $i | cut -d "." -f 1)
for x in "W23-Ab10_R2_001.fastq.gz" \
"W23-DfL_R2_001.fastq.gz" \
"W23-DfK_R2_001.fastq.gz"
do
	SAMP=$(echo $x | cut -d "_" -f 1)
	#BLAST the reads
	blastall -b 5000 -F F -p blastn -d $REF_DIR/$i -i $DIR/${SAMP}_R2.fasta -o $DIR/"BLAST_"$REF"_"$SAMP"_R2.2.out" -m 8
done
done

#This section combines the R1 and R2 hit files for each 
for i in "knob180.fasta" \
"TR1.fasta" \
"CentC.fasta"
do 	
	REF=$(echo $i | cut -d "." -f 1)
for x in "W23-Ab10_R2_001.fastq.gz" \
"W23-DfL_R2_001.fastq.gz" \
"W23-DfK_R2_001.fastq.gz"
do
	SAMP=$(echo $x | cut -d "_" -f 1)
	#merge the two hit files together
	#cat $DIR/"BLAST_"$REF"_"$SAMP"_R1.2.out" $DIR/"BLAST_"$REF"_"$SAMP"_R2.2.out" >> $DIR/"BLAST_"$REF"_"$SAMP"_ALL.2.out"
	#awk 'BEGIN{FS="\t"} {print $1, $7, $8, $4, $11}' OFS='\t' $DIR/"BLAST_"$REF"_"$SAMP"_ALL.2.out" > $DIR/"BLAST_"$REF"_"$SAMP"_ALL.2.bed"
	sort -k 1,1 -k2,2n $DIR/"BLAST_"$REF"_"$SAMP"_ALL.2.bed" > $DIR/"BLAST_"$REF"_"$SAMP"_ALL.2.s.bed"
	bedtools merge -i $DIR/"BLAST_"$REF"_"$SAMP"_ALL.2.s.bed" > $DIR/"BLAST_"$REF"_"$SAMP"_ALL.merge.bed" 

done
done		