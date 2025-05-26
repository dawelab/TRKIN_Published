#module load R/4.4.1-foss-2022b
#install.packages("tidyverse")
library(tidyverse)

setwd("/scratch/mjb51923/TRKIN_CRISPR/out_paper/OrthoFinder")

data <- readLines("/scratch/mjb51923/TRKIN_CRISPR/out_paper/Liftoff/CI66_K10L2_v1.gene.v2.sorted.gff3")
filtered_data <- data[!grepl("^#", data)]

write.table(filtered_data, "CI66_K10L2.genes.edit.gff3", row.names = FALSE, quote = FALSE)

#Removed the header to make this file compatible with R 
K10L2_GFF <- read.delim("CI66_K10L2.genes.edit.gff3", header = FALSE)
colnames(K10L2_GFF) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

#The R1 gene is at 2730186 in the K10L2HiFI genome

#This selects only genes on the K10L2 haplotype
K10L2_GFF <- subset(K10L2_GFF, feature == "gene" & seqname == "K10L2" & start >= 2730186)

#This separates out the gene ID 
K10L2_GFF_GENE_temp1 <- separate(K10L2_GFF, col=attribute, into= c("ID", "biotype", "logic_name"), sep= ';')
K10L2_GFF_GENE_temp2 <- separate(K10L2_GFF_GENE_temp1, col= ID, into=c("Trash", "ID"), sep="=")
K10L2_GFF_GENE_temp3 <- K10L2_GFF_GENE_temp2[,c("seqname", "start", "end", "ID")]
colnames(K10L2_GFF_GENE_temp3) <- c("K10L2_seqname", "K10L2_start", "K10L2_end", "K10L2_ID")
K10L2_GFF_GENE_temp3$K10L2_ID <-  gsub("gene:", "", K10L2_GFF_GENE_temp3$K10L2_ID)
K10L2_GFF_GENE <- K10L2_GFF_GENE_temp3

#This writes the gene names on the K10L2 haplotype
write.table(K10L2_GFF_GENE$K10L2_ID, "CI66_K10L2.K10L2hapGeneNames.txt", row.names = FALSE, quote = FALSE)
