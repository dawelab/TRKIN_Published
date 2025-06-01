module load FastQC/0.11.9-Java-11
module load Trimmomatic/0.39-Java-13

#From doi 10.1101/gad.340679.120
RNA_DIR="/path/to/reads"
OUT_DIR="/scratch/mjb51923/TRKIN_CRISPR/out/K10L2_Trimmed"


#This checks data quality
# for i in K10L2_1_1.fastq \
# K10L2_1_2.fastq \
# K10L2_2_1.fastq \
# K10L2_2_2.fastq \
# K10L2_3_1.fastq \
# K10L2_3_2.fastq \
# K10L2_4_1.fastq \
# K10L2_4_2.fastq \
# K10L2_5_1.fastq \
# K10L2_5_2.fastq \
# N10_1_1.fastq \
# N10_1_2.fastq \
# N10_2_1.fastq \
# N10_2_2.fastq \
# N10_3_1.fastq \
# N10_3_2.fastq
# do
# fastqc $RNA_DIR/$i --extract -o $OUT_DIR
# done

#This trims adapters and low quality sequences
for i in K10L2_1 \
K10L2_2 \
K10L2_3 \
K10L2_4 \
K10L2_5 \
N10_1 \
N10_2 \
N10_3
do
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 12 $RNA_DIR/$i"_1.fastq" $RNA_DIR/$i"_2.fastq" \
$OUT_DIR/$i"_1_paired.fastq" $OUT_DIR/$i"_1_unpaired.fastq"
$OUT_DIR/$i"_2_paired.fastq" $OUT_DIR/$i"_2_unpaired.fastq" \
ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE-2.fa:2:30:10:1:keepBothReads LEADING:3 TRAILING:3 MINLEN:36
done

#This rechecks quality 

for i in K10L2_1_1_paired.fastq \
K10L2_1_2_paired.fastq \
K10L2_2_1_paired.fastq \
K10L2_2_2_paired.fastq \
K10L2_3_1_paired.fastq \
K10L2_3_2_paired.fastq \
K10L2_4_1_paired.fastq \
K10L2_4_2_paired.fastq \
K10L2_5_1_paired.fastq \
K10L2_5_2_paired.fastq \
N10_1_1_paired.fastq \
N10_1_2_paired.fastq \
N10_2_1_paired.fastq \
N10_2_2_paired.fastq \
N10_3_1_paired.fastq \
N10_3_2_paired.fastq
do
fastqc $OUT_DIR/$i --extract -o $OUT_DIR
done
