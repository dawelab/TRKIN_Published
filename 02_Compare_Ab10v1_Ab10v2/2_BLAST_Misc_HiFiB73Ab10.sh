###All of this was run on the command line 

#load the modules necessary from the cluster
module load BLAST+/2.13.0-gompi-2022a

#Define Variables, variables in all caps 
DIR=""
REF="/path/to/B73_Ab10_HiFi_v2.fa"

#make blast database for the B73Ab10 reference genome
#This was run in a previous script 
makeblastdb -in $REF -parse_seqids -dbtype nucl

blastn -num_threads 30 -task "blastn" -outfmt 6 -db $REF -query /scratch/mjb51923/annotations/TR1.fasta -out $DIR/BLAST_TR1_v_HiFiAb10.out
awk '{print $2, $9, $10, $3}' $DIR/BLAST_TR1_v_HiFiAb10.out > $DIR/BLAST_TR1_v_HiFiAb10.bed

blastn -num_threads 30 -task "blastn" -outfmt 6 -db $REF -query /scratch/mjb51923/annotations/knob180.fasta -out $DIR/BLAST_knob180_v_HiFiAb10.out
awk '{print $2, $9, $10, $3}' $DIR/BLAST_knob180_v_HiFiAb10.out > $DIR/BLAST_knob180_v_HiFiAb10.bed

blastn -num_threads 30 -task "blastn" -outfmt 6 -db $REF -query /scratch/mjb51923/annotations/trkin_cds_1.fa -out $DIR/BLAST_trkin_v_HiFiAb10.out
awk '{print $2, $9, $10, $3}' $DIR/BLAST_trkin_v_HiFiAb10.out > $DIR/BLAST_trkin_v_HiFiAb10.bed

blastn -num_threads 30 -task "blastn" -outfmt 6 -db $REF -query /scratch/mjb51923/annotations/Zm00001eb429330.fasta -out $DIR/BLAST_r1_v_HiFiAb10.out
awk '{print $2, $9, $10, $3}' $DIR/BLAST_r1_v_HiFiAb10.out > $DIR/BLAST_r1_v_HiFiAb10.bed

blastn -num_threads 30 -task "blastn" -outfmt 6 -db $REF -query /scratch/mjb51923/ref_genomes/kindr_E9_CDS.fasta -out $DIR/BLAST_kindr_v_HiFiAb10.out
awk '{print $2, $9, $10, $3}' $DIR/BLAST_kindr_v_HiFiAb10.out > $DIR/BLAST_kindr_v_HiFiAb10.bed

