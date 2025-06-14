#module load R/4.4.1-foss-2022b

library(tidyr)
library(dplyr)

setwd("")

#This file comes from indexing the protein fasta in 07.2.2
K10L2HapProt <- read.delim("CI66_K10L2.K10L2hapProtein.fasta.fai", header=FALSE)
names(K10L2HapProt) <- c("name", "length", "offset", "linebases", "linewidth")

#This separates out the name field
K10L2HapProt <- separate(K10L2HapProt, col = name, into=c("gene", "isoform"), sep="\\.", extra = "merge")

#This identifies the longest isoform for each gene
Long <- K10L2HapProt %>% group_by(K10L2HapProt$gene) %>% arrange(-length) %>% slice(1)
Long$ID <- paste(Long$gene, Long$isoform, sep=".")

#Print a list of gene names and numbers of the longest isoforms
write.table(Long$ID, file = "CI66_K10L2.K10L2hapProtein.fasta.K10L2hapGeneNamesLongest.txt", sep="\n", row.names = FALSE, col.names = FALSE, quote = FALSE)
