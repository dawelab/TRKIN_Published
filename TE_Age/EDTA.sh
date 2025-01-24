module load EDTA/2.2.0

REF=/scratch/mjb51923/ref_genomes/Ab10_HiFi_v2_corrected.fa
TE_ANNO=/scratch/mjb51923/TRKIN_CRISPR/out_paper/RepeatMasker/Ab10_HiFi_v2_corrected.fa.ori.out
#Remove the scaffolds
grep "chr" $REF > /scratch/mjb51923/ref_genomes/Ab10_HiFi_v2_corrected.noscaff.fa

#Average mutation rate for maize from https://doi.org/10.1038/s41467-017-02063-5
EDTA.pl --genome /scratch/mjb51923/ref_genomes/Ab10_HiFi_v2_corrected.noscaff.fa --species Maize --step anno --rmout $TE_ANNO -u 0.0000000302