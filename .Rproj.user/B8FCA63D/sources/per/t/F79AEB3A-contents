library(tidyverse)
library(topGO)
library(lintr)
library(lattice)
library(dplyr)

setwd("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/Paper/TRKIN_Paper")

K10L2_GFF <- read.delim("~/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/OrthoFinder/CI66_K10L2_v1.gene.v2.gff3", header = FALSE)
#This drops an unnecessary column at the top
K10L2_GFF <- K10L2_GFF[-c(1),]
colnames(K10L2_GFF) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

#This section alters the K10L2 GFF to be compatible with the OrthoFinder Output
K10L2_GFF_GENE <- subset(K10L2_GFF, feature == "gene")

K10L2_GFF_GENE_temp1 <- separate(K10L2_GFF_GENE, col=attribute, into= c("ID", "biotype", "logic_name"), sep= ';')

K10L2_GFF_GENE_temp1$biotype <-  gsub("Name=", "", K10L2_GFF_GENE_temp1$biotype)

K10L2_GFF_GENE <- K10L2_GFF_GENE_temp1

K10L2_GFF_GENE_f1 <- subset(K10L2_GFF_GENE, seqname == "K10L2")

#This extracts the K10L2 Specific region genes
K10L2_SPEC<- subset(K10L2_GFF_GENE_f1, start >= 7429619 & end <= 25502283)

#This builds the custom gene list
geneNames <- K10L2_GFF_GENE_f1$biotype
geneList <- factor(as.integer(geneNames %in% K10L2_SPEC$biotype))
names(geneList) <- geneNames
str(geneList)

#This builds the custom GO term mapping 
K10L2_entap <- read.delim("~/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/GOTermEnrichment/K10L2_entap_results.tsv")
K10L2_entap$GO_ALL <- paste(K10L2_entap$EggNOG.GO.Biological, K10L2_entap$EggNOG.GO.Cellular, K10L2_entap$EggNOG.GO.Molecular, sep = ",")

K10L2_GO <- K10L2_entap[,c("Query.Sequence" ,"GO_ALL")]
K10L2_GO$GO_ALL <- gsub("NA,", "", K10L2_GO$GO_ALL)
write.table(K10L2_GO, file="K10L2_GOMapping.tbl", quote=FALSE, sep='\t', row.names = FALSE, col.names = FALSE)

#This builds the topgo object
#node size says 10 genes must be assigned that GO to include it
DATA_BP <- new("topGOdata", description="K10L2", ontology = "BP", allGenes = geneList, nodeSize = 10, annot = annFUN.file, file = "K10L2_GOMapping.tbl")
resultFisherBP <- runTest(DATA_BP, algorithm = "classic", statistic = "fisher")
BPRes <- GenTable(DATA_BP, classicFisher = resultFisherBP,
                  ranksOf = "classicFisher", topNodes = 10)
BPRes$Type <- "Biological Process"

DATA_MF <- new("topGOdata", description="K10L2", ontology = "MF", allGenes = geneList, nodeSize = 10, annot = annFUN.file, file = "K10L2_GOMapping.tbl")
resultFisherMF <- runTest(DATA_MF, algorithm = "classic", statistic = "fisher")
MFRes <- GenTable(DATA_MF, classicFisher = resultFisherMF,
                  ranksOf = "classicFisher", topNodes = 10)
MFRes$Type <- "Molecular Function"

DATA_CC <- new("topGOdata", description="K10L2", ontology = "CC", allGenes = geneList, nodeSize = 10, annot = annFUN.file, file = "K10L2_GOMapping.tbl")
resultFisherCC <- runTest(DATA_CC, algorithm = "classic", statistic = "fisher")
CCRes <- GenTable(DATA_CC, classicFisher = resultFisherCC,
                  ranksOf = "classicFisher", topNodes = 10)
CCRes$Type <- "Cellular Component"

AllRes <- rbind(BPRes, MFRes, CCRes)

AllRes$classicFisher <- as.numeric(AllRes$classicFisher)
write.csv(AllRes, "K10L2_GOTermEnrichment.csv", quote = FALSE, row.names = FALSE)

