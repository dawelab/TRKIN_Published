#!/usr/bin/bash
#SBATCH --partition=batch
#SBATCH -J Filt_BamToFastq_CI66HiFi
#SBATCH --output Filt_BamToFastq_CI66HiFi.out
#SBATCH --mem=50GB
#SBATCH --time=5:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=1
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END

module load BEDTools/2.30.0-GCC-12.2.0

module load BamTools/2.5.2-GCC-11.2.0

cd /scratch/mjb51923/raw_reads/CI66_HiFi

bamtools filter -in m84082_240109_064404_s1.hifi_reads.bc2079.bam -out m84082_240109_064404_s1.hifi_reads.bc2079.rq99.bam -tag "rq":">=0.99" 

bedtools bamtofastq -i m84082_240109_064404_s1.hifi_reads.bc2079.rq99.bam -fq CI66.HiFi.rq99.fastq
