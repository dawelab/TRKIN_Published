# These scripts annotate both the Ab10 and K10L2 genomes. 
## In all cases the Ab10 and K10L2 annotations are independent of each other, 
1. These scripts use RepeatMasker and an existing publically available maize TE and repeat library to annotate and mask TEs and simple repeats (https://github.com/oushujun/MTEC)
2. These scripts take all the available RNAseq data for Ab10 and K10L2 respectivly and map them to their respective genomes using HiSat2.
3. These scirpts use StringTie2 to generate a genome guided transcriptome assembly for Ab10 and K10L2 using all available RNAseq data.
4. These scripts use Trinity to generate a denoo transcriptome assembly for Ab10 and K10L2 using all available RNAseq data.
5. These scripts use Braker 3 to annotate both genomes.
6. These sciripts use Pasa to generate a comprehensize transcriptome from the StringTie and Trinity transcriptome assemblies and then polish the Braker annotation 4 times.
7. These scripts use liftoff on the K10L2 haplotype and trkin bearing interknob region on Ab10 to update the annotations in each region.
8. These scripts use AGAT to extract CDS and cDNA sequence from the final annotation.
9. These scripts use EnTAP to produce functional annotations for all annotated genes.
   
