# These scripts compare the expression of trkin between Ab10-I and K10L2
1. These scripts assess the fastq quality, trim the adaptors, and then reasses the trimmed fastq quality for the Ab10-I and K10L2 RNAseq reads.
2. These scripts align the trimmed RNAseq reads to the Ab10_v1 assembly (Liu et al. 2020) for Ab10-I and K10L2 using HiSat2.
3. These scripts run the R package feature count to determine the raw expression of annotated genes from the v1 annotation (Liu et al. 2020).
4. This calculates TPM and plots the data. 
