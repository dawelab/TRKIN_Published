#!/usr/bin/bash
#SBATCH --partition=highmem_p 
#SBATCH -J HISAT2_Ab10I
#SBATCH --output HISAT2_Ab10I.out
#SBATCH --mem=500GB
#SBATCH --time=48:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=19
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END

#load modules
module load HISAT2/3n-20201216-gompi-2022a
module load SAMtools/1.16.1-GCC-11.3.0

#Define Variables
READ_DIR="/scratch/mjb51923/TRKIN_CRISPR/out/Ab10I_Trimmed"
OUT_DIR="/scratch/mjb51923/TRKIN_CRISPR/out/RNA_Aln"
REF="/scratch/mjb51923/ref_genomes/Zm-B73_AB10-REFERENCE-NAM-1.0.fa"

#Index the reference
#hisat2-build $REF $REF

for i in Ab10IMMR_1 \
Ab10IMMR_2 \
Ab10IMMR_3 
do
hisat2 -p 12 -x $REF -1 $READ_DIR/$i"_1_paired.fastq" -2 $READ_DIR/$i"_2_paired.fastq" -S $OUT_DIR/$i"_v_B73Ab10.bam"
samtools sort $OUT_DIR/$i"_v_B73Ab10.bam" -o $OUT_DIR/$i"_v_B73Ab10.s.bam"
samtools index $OUT_DIR/$i"_v_B73Ab10.s.bam"
done

for i in N10_1 \
N10_2 \
N10_3
do
hisat2 -p 12 -x $REF -1 $READ_DIR/$i"_1_paired.fastq" -2 $READ_DIR/$i"_2_paired.fastq" -S $OUT_DIR/$i"_Ab10Sib_v_B73Ab10.bam"
samtools sort $OUT_DIR/$i"_Ab10Sib_v_B73Ab10.bam" -o $OUT_DIR/$i"_Ab10Sib_v_B73Ab10.s.bam"
samtools index $OUT_DIR/$i"_Ab10Sib_v_B73Ab10.s.bam"
done
