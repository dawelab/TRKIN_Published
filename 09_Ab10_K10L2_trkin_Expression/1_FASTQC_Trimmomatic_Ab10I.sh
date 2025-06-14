module load FastQC/0.11.9-Java-11
module load Trimmomatic/0.39-Java-13

#From doi 10.1101/gad.340679.120
RNA_DIR="/path/to/reads"
OUT_DIR=""


#This checks data quality
for i in Ab10IMMR_1_1 \
Ab10IMMR_1_2 \
Ab10IMMR_2_1 \
Ab10IMMR_2_2 \
Ab10IMMR_3_1 \
Ab10IMMR_3_2 \
N10_1_1 \
N10_1_2 \
N10_2_1 \
N10_2_2 \
N10_3_1 \
N10_3_2
do
fastqc $RNA_DIR/$i.fastq --extract -o $OUT_DIR
done

#This trims adapters and low quality sequences
for i in Ab10IMMR_1 \
Ab10IMMR_2 \
Ab10IMMR_3 \
N10_1 \
N10_2 \
N10_3
do
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 12 $RNA_DIR/$i"_1.fastq" $RNA_DIR/$i"_2.fastq" \
$OUT_DIR/$i"_1_paired.fastq" $OUT_DIR/$i"_1_unpaired.fastq" \
$OUT_DIR/$i"_2_paired.fastq" $OUT_DIR/$i"_2_unpaired.fastq" \
ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE-2.fa:2:30:10:1:keepBothReads LEADING:3 TRAILING:3 MINLEN:36
done

#This rechecks quality 
for i in Ab10IMMR_1_1_paired.fastq \
Ab10IMMR_1_2_paired.fastq \
Ab10IMMR_2_1_paired.fastq \
Ab10IMMR_2_2_paired.fastq \
Ab10IMMR_3_1_paired.fastq \
Ab10IMMR_3_2_paired.fastq \
N10_1_1_paired.fastq \
N10_1_2_paired.fastq \
N10_2_1_paired.fastq \
N10_2_2_paired.fastq \
N10_3_1_paired.fastq \
N10_3_2_paired.fastq
do
fastqc $OUT_DIR/$i --extract -o $OUT_DIR
done
