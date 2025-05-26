t#Load the modile
module load hifiasm/0.19.6-GCCcore-11.3.0

OUTDIR=""
READS="/path/to/rawreads/CI66.HiFi.rq99.fastq"

#Enter output directory
cd $OUT_DIR

#Perform assemble
#write commands just write more files for checks if you need them later. -t is threads, -u disables post-joining, which can lead to misassemblies. -l0 indicates that this is a homozygote
hifiasm --write-ec --write-paf -u 0 -o CI66_rq99.asm -l0 -t 20 $READS 
