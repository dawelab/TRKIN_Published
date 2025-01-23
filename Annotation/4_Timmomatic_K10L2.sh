#!/usr/bin/bash
#SBATCH --partition=batch 
#SBATCH -J FASTQC_Trimmomatic_K10L2
#SBATCH --output FASTQC_Trimmomatic_K10L2.%A-%a.out
#SBATCH --mem=100GB
#SBATCH --time=10:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=12
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --array=1-10

module load FastQC/0.11.9-Java-11
module load Trimmomatic/0.39-Java-13

RNA_DIR=/scratch/mjb51923/raw_reads/RNA/K10L2
OUT_DIR="/scratch/mjb51923/TRKIN_CRISPR/out_paper/Trimmed_RNAseq_K10L2"
#qmkdir $OUT_DIR
i=$(awk NR==${SLURM_ARRAY_TASK_ID}'{print $1}' /scratch/mjb51923/TRKIN_CRISPR/TRKIN_Published/Annotation/4_K10L2_RNASeq_ReadFiles.txt)

#Unzips the file
gunzip $RNA_DIR/$i"_R1.fq.gz"
gunzip $RNA_DIR/$i"_R2.fq.gz"

#This checks data quality
fastqc $RNA_DIR/$i"_R1.fq" --extract -o $OUT_DIR
fastqc $RNA_DIR/$i"_R2.fq" --extract -o $OUT_DIR

#This trims adapters and low quality sequences
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 12 $RNA_DIR/$i"_R1.fq" $RNA_DIR/$i"_R2.fq" \
$OUT_DIR/$i"_R1_paired.fq" $OUT_DIR/$i"_R1_unpaired.fq" \
$OUT_DIR/$i"_R2_paired.fq" $OUT_DIR/$i"_R2_unpaired.fq" \
ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE-2.fa:2:30:10:1:keepBothReads LEADING:3 TRAILING:3 MINLEN:36

#This rechecks quality 
fastqc $OUT_DIR/$i"_R1_paired.fq" --extract -o $OUT_DIR
fastqc $OUT_DIR/$i"_R2_paired.fq" --extract -o $OUT_DIR
