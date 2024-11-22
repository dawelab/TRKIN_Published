#module load R/4.4.1-foss-2022b

library(tidyr)
library(dplyr)

setwd("/scratch/mjb51923/TRKIN_CRISPR/out_paper/OrthoFinder")

Ab10HapProt <- read.delim("HiFiAb10.Ab10hapProtein.fasta.fai", header=FALSE)

names(Ab10HapProt) <- c("name", "length", "offset", "linebases", "linewidth")


Ab10HapProt <- separate(Ab10HapProt, col = name, into=c("gene", "isoform"), sep="\\.", extra = "merge")

genes <- unique(Ab10HapProt$gene)

Long <- Ab10HapProt %>% group_by(Ab10HapProt$gene) %>% arrange(-length) %>% slice(1)
length(Long$gene)

Long$ID <- paste(Long$gene, Long$isoform, sep=".")

#Print a list of gene names and numbers of the longest isoforms
write.table(Long$ID, file = "HiFiAb10.Ab10hapProtein.fasta.Ab10hapGeneNamesLongest.txt", sep="\n", row.names = FALSE, col.names = FALSE, quote = FALSE)
