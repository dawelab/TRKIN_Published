module load AGAT/1.1.0
#Convert the Stringtie gtf to a gff file 
gunzip /scratch/mjb51923/TRKIN_CRISPR/out_paper/Liftoff/Ab10_HiFi_v2_corrected.gene.v2.gff3.gz
agat_convert_sp_gff2gtf.pl --gff /scratch/mjb51923/TRKIN_CRISPR/out_paper/Liftoff/Ab10_HiFi_v2_corrected.gene.v2.gff3 -o /scratch/mjb51923/TRKIN_CRISPR/out_paper/Liftoff/Ab10_HiFi_v2_corrected.gene.v2.gtf
