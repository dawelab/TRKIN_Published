#!/usr/bin/bash
#SBATCH --partition=highmem_p 
#SBATCH -J HISAT2_K10L2
#SBATCH --output HISAT2_K10L2.out
#SBATCH --mem=500GB
#SBATCH --time=48:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=19
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END

#load modules
module load HISAT2/2.1.0-foss-2019b
module load SAMtools/1.9-GCC-8.3.0

#Define Variables
READ_DIR="/scratch/mjb51923/TRKIN_CRISPR/out/K10L2_Trimmed"
OUT_DIR="/scratch/mjb51923/TRKIN_CRISPR/out/RNA_Aln"
REF="/scratch/mjb51923/ref_genomes/Zm-B73_AB10-REFERENCE-NAM-1.0.fa"

#Index the reference
hisat2-build $REF $REF

for i in K10L2_1 \
K10L2_2 \
K10L2_3 \
K10L2_4 \
K10L2_5 \
N10_1 \
N10_2 \
N10_3
do
hisat2 -p 12 -x $REF -1 $OUT_DIR/$i"_1_paired.fastq" -2 $OUT_DIR/$i"_2_paired.fastq" -S $OUT_DIR/$i"_v_B73Ab10.bam"
samtools sort $OUT_DIR/$i"_v_B73Ab10.bam" -o $OUT_DIR/$i"_v_B73Ab10.s.bam"
samtools index $OUT_DIR/$i"_v_B73Ab10.s.bam"
done
