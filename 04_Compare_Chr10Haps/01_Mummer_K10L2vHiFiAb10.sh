#Load Modules
module load MUMmer/4.0.0rc1-GCCcore-11.3.0
module load SAMtools/1.10-GCC-8.3.0

#Define Variables
DIR=""

#-t is threads, -l is minimum length generate maximal unique matches, -maxmatch has it compute all maximal matches not just unique ones, -b computes both forward and reverse complements, I believe that
cd $DIR
nucmer --maxmatch -p K10L2_v_Ab10_nucmer -t 30 -l 300 $DIR/CI66_K10L2_v1.K10L2.fasta $DIR/Ab10_HiFi_v2_corrected_Ab10Hap.Ab10.fasta

#This outputs the cords
show-coords -c -l -r -T K10L2_v_Ab10_nucmer.delta > K10L2_v_Ab10_nucmer.coords
