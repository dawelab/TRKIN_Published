library(tidyverse)
library(topGO)
library(lintr)
library(lattice)
library(dplyr)

setwd("")

#I manually removed all commented out lines for compatibility with R

This comes from 05.08
Ab10_GFF <- read.delim("B73_Ab10_HiFi_v2.gene.gff3", header = FALSE)
colnames(Ab10_GFF) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

#This section alters the Ab10 GFF to be compatible with the OrthoFinder Output
Ab10_GFF_GENE <- subset(Ab10_GFF, feature == "gene")
Ab10_GFF_GENE_temp1 <- separate(Ab10_GFF_GENE, col=attribute, into= c("ID", "biotype", "logic_name"), sep= ';')
Ab10_GFF_GENE_temp1$biotype <-  gsub("Name=", "", Ab10_GFF_GENE_temp1$biotype)
Ab10_GFF_GENE <- Ab10_GFF_GENE_temp1

#remove the scaffolds
Ab10_GFF_GENE_f1 <- subset(Ab10_GFF_GENE, seqname == "chr1" | seqname == "chr2" | seqname == "chr3" | seqname == "chr4" | seqname == "chr5" | seqname == "chr6" | seqname == "chr7" | seqname == "chr8" | seqname == "chr9" | seqname == "chr10")

#Select only chr10
Ab10_GFF_GENE_f2 <- subset(Ab10_GFF_GENE_f1, seqname == "chr10")

#This extracts the Ab10 Specific region genes
Ab10_SPEC <- subset(Ab10_GFF_GENE_f2, start >= 141187279 & ((start >= 142472000 & end <= 153145000) | (start >= 168235374)))

#This builds the custom gene list
geneNames <- Ab10_GFF_GENE_f1$biotype
geneList <- factor(as.integer(geneNames %in% Ab10_SPEC$biotype))
names(geneList) <- geneNames
str(geneList)

#This builds the custom GO term mapping 
#This file is derived from step 05.12
Ab10_entap <- read.delim("Ab10_entap_results.tsv")
Ab10_entap$GO_ALL <- paste(Ab10_entap$EggNOG.GO.Biological, Ab10_entap$EggNOG.GO.Cellular, Ab10_entap$EggNOG.GO.Molecular, sep = ",")

Ab10_GO <- Ab10_entap[,c("Query.Sequence" ,"GO_ALL")]
Ab10_GO$GO_ALL <- gsub("NA,", "", Ab10_GO$GO_ALL)
write.table(Ab10_GO, file="Ab10_GOMapping.tbl", quote=FALSE, sep='\t', row.names = FALSE, col.names = FALSE)

#This builds the topgo object for the Biological, Cellular, and Molecular categories
#node size requires 10 genes must be assigned that GO to include it
DATA_BP <- new("topGOdata", description="Ab10", ontology = "BP", allGenes = geneList, nodeSize = 10, annot = annFUN.file, file = "Ab10_GOMapping.tbl")
resultFisherBP <- runTest(DATA_BP, algorithm = "classic", statistic = "fisher")
BPRes <- GenTable(DATA_BP, classicFisher = resultFisherBP,
                  ranksOf = "classicFisher", topNodes = 10)
BPRes$Type <- "Biological Process"

DATA_MF <- new("topGOdata", description="Ab10", ontology = "MF", allGenes = geneList, nodeSize = 10, annot = annFUN.file, file = "Ab10_GOMapping.tbl")
resultFisherMF <- runTest(DATA_MF, algorithm = "classic", statistic = "fisher")
MFRes <- GenTable(DATA_MF, classicFisher = resultFisherMF,
                  ranksOf = "classicFisher", topNodes = 10)
MFRes$Type <- "Molecular Function"

DATA_CC <- new("topGOdata", description="Ab10", ontology = "CC", allGenes = geneList, nodeSize = 10, annot = annFUN.file, file = "Ab10_GOMapping.tbl")
resultFisherCC <- runTest(DATA_CC, algorithm = "classic", statistic = "fisher")
CCRes <- GenTable(DATA_CC, classicFisher = resultFisherCC,
                  ranksOf = "classicFisher", topNodes = 10)
CCRes$Type <- "Cellular Component"

#This brings together the results of all the Goterm erichment analyses
AllRes <- rbind(BPRes, MFRes, CCRes)
AllRes$classicFisher <- as.numeric(AllRes$classicFisher)

#This writes out the GO term enrichment 
write.csv(AllRes, "Ab10_GOTermEnrichment.csv", quote = FALSE, row.names = FALSE)

