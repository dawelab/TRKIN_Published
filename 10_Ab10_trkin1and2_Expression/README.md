# These scripts calculate the transcripts per million (TPM) of Ab10 trkin1 and Ab10 trkin2 
1. Download the necessary RNAseq data
2. Align the RNA seq data to B73_Ab10_v2
3. Convert the gff file generated in step 05 to a gtf for compatibility with feature counts
4. Count the number of reads overlapping  each exon
5. Convert feature count to transcripts per million, isolate only the exons that are different between trkin1 and trkin2, and plot the results. 
