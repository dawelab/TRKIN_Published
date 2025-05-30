#Define Variables
DIR=""

#Load Modules
module load MUMmer/4.0.0rc1-GCCcore-11.3.0
module load SAMtools/1.10-GCC-8.3.0

#-t is threads, -l is minimum length generate maximal unique matches, -maxmatch has it computer all maximal matches not just unique ones, -b computes both forward and reverse complements, I believe that --nosimplify ignore all hits with the exact same start coordinate in both sequences  (mums)
cd $DIR
nucmer --nosimplify --maxmatch -p K10L2_v_K10L2_nucmer -t 30 -l 300 $DIR/CI66_K10L2_v1.K10L2.fasta $DIR/CI66_K10L2_v1.K10L2.fasta

#This outputs the cords
show-coords -c -l -r -T K10L2_v_K10L2_nucmer.delta > K10L2_v_K10L2_nucmer.coords
