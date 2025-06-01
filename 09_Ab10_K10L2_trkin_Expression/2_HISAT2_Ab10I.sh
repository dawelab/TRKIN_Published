#load modules
module load HISAT2/3n-20201216-gompi-2022a
module load SAMtools/1.16.1-GCC-11.3.0

#From doi 10.1101/gad.340679.120
RNA_DIR="/path/to/reads"
OUT_DIR=""
#From doi 10.1186/s13059-020-02029-9
REF="Zm-B73_AB10-REFERENCE-NAM-1.0.fa"

#Index the reference
hisat2-build $REF $REF

#Align the Ab10 and N10 sequences to the reference, sort and index the bam files
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
