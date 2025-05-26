#load the modules necessary from the cluster
module load BLAST+/2.13.0-gompi-2022a


#Define Variables, variables in all caps 
OUTDIR=""
NAME="CI66_rq99.asm"

#make blast database for the CI66 K10L2 scaffolds.
makeblastdb -in $OUTDIR/${NAME}.bp.p_ctg.gfa.fasta -parse_seqids -dbtype nucl

#Identify scaffolds containing the K10L2 haplotype

#blast trkin against the scaffolds. trkin cds can be found under NCBI GenBank: MT459824.1
blastn -num_threads 30 -task "blastn" -outfmt 6 -db $OUTDIR/${NAME}.bp.p_ctg.gfa.fasta -query trkin_cds_1.fa -out $OUTDIR/BLAST_trkin_v_CI66Scaff.out

#Blast colored1 against the CI66 scaffolds. The colored 1 gene can be found by searching the gene ID Zm00001eb429330 on https://www.maizegdb.org/
blastn -num_threads 30 -task "blastn" -outfmt 6 -db $OUTDIR/${NAME}.bp.p_ctg.gfa.fasta -query Zm00001eb429330.fasta -out $OUTDIR/BLAST_r1_v_CI66Scaff.out
