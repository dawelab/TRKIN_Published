library(tidyverse)
library(topGO)
library(lintr)
library(lattice)
library(dplyr)

setwd("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/Paper/TRKIN_Paper")

Ab10_GFF <- read.delim("~/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/OrthoFinder/Ab10_HiFi_v2_corrected.gene.v2.gff3", header = FALSE)
colnames(Ab10_GFF) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

#This section alters the Ab10 GFF to be compatible with the OrthoFinder Output
Ab10_GFF_GENE <- subset(Ab10_GFF, feature == "gene")

Ab10_GFF_GENE_temp1 <- separate(Ab10_GFF_GENE, col=attribute, into= c("ID", "biotype", "logic_name"), sep= ';')

Ab10_GFF_GENE_temp1$biotype <-  gsub("Name=", "", Ab10_GFF_GENE_temp1$biotype)

Ab10_GFF_GENE <- Ab10_GFF_GENE_temp1

#remove the scaffolds as I'm not sure where they come from so they might obscure the analysis
Ab10_GFF_GENE_f1 <- subset(Ab10_GFF_GENE, seqname == "chr1" | seqname == "chr2" | seqname == "chr3" | seqname == "chr4" | seqname == "chr5" | seqname == "chr6" | seqname == "chr7" | seqname == "chr8" | seqname == "chr9" | seqname == "chr10")

Ab10_GFF_GENE_f2 <- subset(Ab10_GFF_GENE_f1, seqname == "chr10")

#This extracts the Ab10 Specific region genes
Ab10_SPEC <- subset(Ab10_GFF_GENE_f2, start >= 141187279 & ((start >= 142472000 & end <= 153145000) | (start >= 168235374)))

#This drops 8 copies of the 9 copies of nrpd2/e2 as they seem to be skewing the analysis 
#Ab10_SPEC <- subset(Ab10_SPEC, biotype != "g5851.t1" & biotype != "g5853.t1" & biotype != "g5854.t1" & biotype != "g5856.t1" & biotype != "g5858.t1" & biotype != "g5860.t1" & biotype != "g5862.t1" & biotype != "g5863.t1" & biotype != "g5880.t1" & biotype != "g5879.t1" & biotype != "g5874.t1")

#This builds the custom gene list
geneNames <- Ab10_GFF_GENE_f1$biotype
geneList <- factor(as.integer(geneNames %in% Ab10_SPEC$biotype))
names(geneList) <- geneNames
str(geneList)

#This builds the custom GO term mapping 
Ab10_entap <- read.delim("~/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/GOTermEnrichment/Ab10_entap_results.tsv")
Ab10_entap$GO_ALL <- paste(Ab10_entap$EggNOG.GO.Biological, Ab10_entap$EggNOG.GO.Cellular, Ab10_entap$EggNOG.GO.Molecular, sep = ",")

Ab10_GO <- Ab10_entap[,c("Query.Sequence" ,"GO_ALL")]
Ab10_GO$GO_ALL <- gsub("NA,", "", Ab10_GO$GO_ALL)
write.table(Ab10_GO, file="Ab10_GOMapping.tbl", quote=FALSE, sep='\t', row.names = FALSE, col.names = FALSE)

#This builds the topgo object
#node size says 10 genes must be assigned that GO to include it
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

AllRes <- rbind(BPRes, MFRes, CCRes)

AllRes$classicFisher <- as.numeric(AllRes$classicFisher)
write.csv(AllRes, "Ab10_GOTermEnrichment.csv", quote = FALSE, row.names = FALSE)

