#!/usr/bin/bash
#SBATCH --partition=batch
#SBATCH -J HISAT2_B73_Ab10
#SBATCH --output HISAT2_B73_Ab10.%A-%a.out
#SBATCH --mem=300GB
#SBATCH --time=72:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=19
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --array=1-20

#load modules
module load HISAT2/3n-20201216-gompi-2022a
module load SAMtools/1.16.1-GCC-11.3.0

#Define Variables
READ_DIR="/scratch/mjb51923/raw_reads/RNA/Gapless_B73-Ab10I"
OUT_DIR="/scratch/mjb51923/TRKIN_CRISPR/out_paper/B73_Ab10_RNA"
REF=/scratch/mjb51923/ref_genomes/Ab10_HiFi_v2_corrected.fa
#This pulls info from the array job
FILE=$(awk NR==${SLURM_ARRAY_TASK_ID}'{print $1}' /scratch/mjb51923/TRKIN_CRISPR/TRKIN_Published/trkin1and2_Expression/Tissue_List.txt)

#Build the reference
#hisat2-build $REF $REF

#Align the reads
hisat2 -p 12 -x $REF -1 $READ_DIR/$FILE"_R1.fq" -2 $READ_DIR/$FILE"_R2.fq" -S $OUT_DIR/$FILE"_v_B73Ab10v2.bam"
#sort the output
samtools sort $OUT_DIR/$FILE"_v_B73Ab10v2.bam" -o $OUT_DIR/$FILE"_v_B73Ab10v2.s.bam"
#index the output
samtools index $OUT_DIR/$FILE"_v_B73Ab10v2.s.bam"
#filter the output to only reads with no mismatches
samtools view -b -e '[NM]==0' -O BAM -o $OUT_DIR/$FILE"_v_B73Ab10v2.nomismatch.s.bam" $OUT_DIR/$FILE"_v_B73Ab10v2.s.bam"
#index the output
samtools index $OUT_DIR/$FILE"_v_B73Ab10v2.nomismatch.s.bam"
#filter the output to a mapQ of 20 with no mismatches
samtools view -b -q 20 -e '[NM]==0' -O BAM -o $OUT_DIR/$FILE"_v_B73Ab10v2.MAPQ20nomismatch.s.bam" $OUT_DIR/$FILE"_v_B73Ab10v2.s.bam"
#index the output
samtools index $OUT_DIR/$FILE"_v_B73Ab10v2.MAPQ20nomismatch.s.bam"
#filter the output to a mapQ of 10 with no mismatches
samtools view -b -q 10 -e '[NM]==0' -O BAM -o $OUT_DIR/$FILE"_v_B73Ab10v2.MAPQ10nomismatch.s.bam" $OUT_DIR/$FILE"_v_B73Ab10v2.s.bam"
#output only chromosome 10 of the raw unfiltered file
samtools view -b $OUT_DIR/$FILE"_v_B73Ab10v2.s.bam" "chr10" > $OUT_DIR/$FILE"_v_B73Ab10v2.chr10.s.bam"

