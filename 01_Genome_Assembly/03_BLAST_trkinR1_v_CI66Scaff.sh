#load the modules necessary from the cluster
module load BLAST+/2.13.0-gompi-2022a


#Define Variables, variables in all caps 
DIR=""
NAME="CI66_rq99.asm"

#make blast database for the CI66 K10L2 scaffolds.
makeblastdb -in $DIR/${NAME}.bp.p_ctg.gfa.fasta -parse_seqids -dbtype nucl

#Identify scaffolds containing the K10L2 haplotype

#blast trkin against the scaffolds 
blastn -num_threads 30 -task "blastn" -outfmt 6 -db $DIR/${NAME}.bp.p_ctg.gfa.fasta -query /scratch/mjb51923/annotations/trkin_cds_1.fa -out $DIR/BLAST_trkin_v_CI66Scaff.out

#Blast colored1 against the CI66 scaffolds
blastn -num_threads 30 -task "blastn" -outfmt 6 -db $DIR/${NAME}.bp.p_ctg.gfa.fasta -query /scratch/mjb51923/annotations/Zm00001eb429330.fasta -out $DIR/BLAST_r1_v_CI66Scaff.out
