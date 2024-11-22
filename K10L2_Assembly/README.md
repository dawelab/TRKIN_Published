# These scripts assemble the K10L2 haplotype from PacBio HiFi data derived from CI66 DNA. 
1. Filter the PacBio HiFi reads and convert them to a fastq formet. 
2. Use Hifiasm to assemble raw reads into contigs
3. Identify contigs with homology to the colored1 (R1) and trkin genes.
4. I manually generated this agp file to assemble the colored1 (R1) bearing contig and the trkin bearing contig into a larger contig with an N gap between them.
5. Assemble the K10L2 contig accoring to the agp described in 4.
6. Identify relavant features of the K10L2 haplotype using BLAST. 
