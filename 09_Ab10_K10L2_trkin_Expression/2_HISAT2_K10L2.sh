#load modules
module load HISAT2/2.1.0-foss-2019b
module load SAMtools/1.9-GCC-8.3.0

#From doi 10.1101/gad.340679.120
RNA_DIR="/path/to/reads"
OUT_DIR=""
#From doi 10.1186/s13059-020-02029-9
REF="Zm-B73_AB10-REFERENCE-NAM-1.0.fa"

#Index the reference
hisat2-build $REF $REF

#Align the K10L2 and N10 sequences to the reference, sort and index the bam files
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
