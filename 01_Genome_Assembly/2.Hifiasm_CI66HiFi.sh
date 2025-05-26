#!/usr/bin/bash
#SBATCH --partition=batch
#SBATCH -J Hifiasm_CI66HiFi
#SBATCH --output Hifiasm_CI66HiFi.out
#SBATCH --mem=300GB
#SBATCH --time=96:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=20
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END

module load hifiasm/0.19.6-GCCcore-11.3.0

OUT_DIR="/scratch/mjb51923/CI66_Assembly/out"
READS="/scratch/mjb51923/raw_reads/CI66_HiFi/CI66.HiFi.rq99.fastq"

cd $OUT_DIR

#write commands just write more files for checks if you need them later. -t is threads, -u disables post-joining, which can lead to misassemblies. -l0 indicates that this is a homozygote
hifiasm --write-ec --write-paf -u 0 -o CI66_rq99.asm -l0 -t 20 $READS 
