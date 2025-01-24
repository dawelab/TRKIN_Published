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
OUT_DIR="/scratch/mjb51923/TRKIN_CRISPR/out/B73_Ab10_RNA"
REF="/scratch/mjb51923/ref_genomes/Zm-B73_AB10-REFERENCE-NAM-1.0.fa"
#This pulls info from the array job
FILE=$(awk NR==${SLURM_ARRAY_TASK_ID}'{print $1}' $OUT_DIR/Tissue_List.txt)

#Build the reference
#hisat2-build $REF $REF

#Align the reads
hisat2 -p 12 -x $REF -1 $READ_DIR/$FILE"_R1.fq.gz" -2 $READ_DIR/$FILE"_R2.fq.gz" -S $OUT_DIR/$FILE"_v_B73Ab10.bam"
#sort the output
samtools sort $OUT_DIR/$FILE"_v_B73Ab10.bam" -o $OUT_DIR/$FILE"_v_B73Ab10.s.bam"
#index the output
samtools index $OUT_DIR/$FILE"_v_B73Ab10.s.bam"
#filter the output to only reads with no mismatches
samtools view -b -e '[NM]==0' -O BAM -o $OUT_DIR/$FILE"_v_B73Ab10.nomismatch.s.bam" $OUT_DIR/$FILE"_v_B73Ab10.s.bam"
#index the output
samtools index $OUT_DIR/$FILE"_v_B73Ab10.nomismatch.s.bam"
#filter the output to a mapQ of 20 with no mismatches
samtools view -b -q 20 -e '[NM]==0' -O BAM -o $OUT_DIR/$FILE"_v_B73Ab10.MAPQ20nomismatch.s.bam" $OUT_DIR/$FILE"_v_B73Ab10.s.bam"
#index the output
samtools index $OUT_DIR/$FILE"_v_B73Ab10.MAPQ20nomismatch.s.bam"
#filter the output to a mapQ of 20 with no mismatches
samtools view -b -q 10 -e '[NM]==0' -O BAM -o $OUT_DIR/$FILE"_v_B73Ab10.MAPQ1nomismatch.s.bam" $OUT_DIR/$FILE"_v_B73Ab10.s.bam"

samtools view -b $OUT_DIR/$FILE"_v_B73Ab10.s.bam" "chr10" > $OUT_DIR/$FILE"_v_B73Ab10.chr10.s.bam"
