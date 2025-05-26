#Load modules
module load BEDTools/2.30.0-GCC-12.2.0
module load BamTools/2.5.2-GCC-11.2.0

#enter directory with CI66 Raw Reads located on NCBI SRA under bioproject PRJNA1254310
cd /path/to/CI66_HiFi_Reads 

#Filter the reads
bamtools filter -in m84082_240109_064404_s1.hifi_reads.bc2079.bam -out m84082_240109_064404_s1.hifi_reads.bc2079.rq99.bam -tag "rq":">=0.99" 

#Convert the bam file to a fastq
bedtools bamtofastq -i m84082_240109_064404_s1.hifi_reads.bc2079.rq99.bam -fq CI66.HiFi.rq99.fastq
