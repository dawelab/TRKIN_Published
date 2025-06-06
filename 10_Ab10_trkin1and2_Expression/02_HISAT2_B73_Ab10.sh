#load modules
module load HISAT2/3n-20201216-gompi-2022a
module load SAMtools/1.16.1-GCC-11.3.0

#Define Variables
READ_DIR="/path/to/reads/from/step1"
OUT_DIR="/"
REF=B73_Ab10_HiFi_v2.fa
#This pulls info from the array job
FILE=$(awk NR==${SLURM_ARRAY_TASK_ID}'{print $1}' Tissue_List.txt)

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

