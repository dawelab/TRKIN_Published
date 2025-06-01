#module load R/4.4.1-foss-2022b

library(tidyr)
library(dplyr)

setwd("")

#This file is from https://www.maizegdb.org/genome/assembly/Zm-B73-REFERENCE-NAM-5.0
B73Prot <- read.delim("Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.protein.fa.fai", header=FALSE)
names(B73Prot) <- c("name", "length", "offset", "linebases", "linewidth")

#This separates out the name field
B73Prot <- separate(B73Prot, col = name, into=c("gene", "isoform"), sep="_")

#This isolates the longest isoform for each gene
Long <- B73Prot %>% group_by(B73Prot$gene) %>% arrange(-length) %>% slice(1)
Long$ID <- paste(Long$gene, Long$isoform, sep="_")

#Print a list of gene names and numbers of the longest isoforms
write.table(Long$ID, file = "Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.protei.GeneNamesLongest.txt", sep="\n", row.names = FALSE, col.names = FALSE, quote = FALSE)
