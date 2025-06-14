library(topGO)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(grid)
library(egg)
library(tidyverse)

setwd("")


Ab10_GFF <- read.delim("B73_Ab10_HiFi_v2.gene.gff3", header = FALSE)
#This drops an unnecessary column at the top
Ab10_GFF <- Ab10_GFF[-c(1),]
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
Ab10_SPEC<- subset(Ab10_GFF_GENE_f2, start >= 141187279 & ((start >= 142472000 & end <= 153145000) | (start >= 168235374)))

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
DATA_MF <- new("topGOdata", description="Ab10", ontology = "MF", allGenes = geneList, nodeSize = 10, annot = annFUN.file, file = "Ab10_GOMapping.tbl")
DATA_CC <- new("topGOdata", description="Ab10", ontology = "CC", allGenes = geneList, nodeSize = 10, annot = annFUN.file, file = "Ab10_GOMapping.tbl")

#This loads in the 10 most enriched GO terms for each category
ENR <- read.csv("1.2_Ab10_GOTermEnrichment.csv")
ENR$Type <- gsub("Biological Process", "DATA_BP", ENR$Type)
ENR$Type <- gsub("Cellular Component", "DATA_CC", ENR$Type)
ENR$Type <- gsub("Molecular Function", "DATA_MF", ENR$Type)

#This loads in the entap results, this file is derived from step 05.12
Ab10_entap <- read.delim("Ab10_entap_results.tsv")
Ab10_entap$GO_ALL <- paste(Ab10_entap$EggNOG.GO.Biological, Ab10_entap$EggNOG.GO.Cellular, Ab10_entap$EggNOG.GO.Molecular, sep = ",")

#This merges the entap results and GFF file
MERGE <- merge(Ab10_SPEC, Ab10_entap, by.x="biotype", by.y="Query.Sequence")
MERGE_GO_ALL <- data.frame("chr"=c(), "start"=c(), "end"=c(), "term"=c(), "enrich"=c(), "p"=c(), "logp"=c(), "type"=c() )

#This function identifies the genes associated with each GO term
GO2GENES <- function(GO, DATA, i) {
  #Extract a list of genes associated with this term 
  LISTGO <- genesInTerm(get(DATA), whichGO = GO)
  #Get the GO Description
  TERM <- ENR[ENR$GO.ID %in% c(GO), 2]
  #Extract gene IDs with the GO terms
  MERGE_GO <- MERGE[MERGE$biotype %in% LISTGO[1][[1]],c("seqname", "start", "end")]
  #Assign the GO term 
  MERGE_GO$Term <- TERM
  #Add the fold enrichment
  MERGE_GO$enrich <- rep(paste(ENR[i,"Annotated"]/ENR[i,"Expected"]), nrow(MERGE_GO))
  #Add the p value
  MERGE_GO$p <- ENR[i,"classicFisher"]
  #Add the p value
  MERGE_GO$logp <- -log(ENR[i,"classicFisher"])
  #Add the axis label color
  if(DATA=="DATA_BP"){MERGE_GO$labcol ="Biological\n Process"}
  if(DATA=="DATA_MF"){MERGE_GO$labcol ="Molecular\n Function"}
  if(DATA=="DATA_CC"){MERGE_GO$labcol ="Cellular\n Component"}
  #Define the column names
  names(MERGE_GO) <- c("chr", "start", "end", "term", "enrich", "p", "logp", "type")
  #Add these columns to the final data frame
  MERGE_GO_ALL <<- rbind(MERGE_GO_ALL, MERGE_GO)
}

for(i in 1:length(ENR$GO.ID)) {
  GO2GENES(ENR[i,1], ENR[i,7], i)
}

#This sorts the data by order of enrichment
MERGE_GO_ALL <- MERGE_GO_ALL[order(MERGE_GO_ALL$enrich, decreasing = TRUE),]
#This defines it as a factor to maintain the ordering
MERGE_GO_ALL$term <- factor(MERGE_GO_ALL$term, levels=rev(unique(MERGE_GO_ALL$term)))

#This plots the GO term enrichment
pdf("Ab10_GOterm_Enrichment.pdf", height=10.5, width=7)
p <- ggplot() +
  geom_point(data=MERGE_GO_ALL, aes(x=term, y=start, color=as.numeric(enrich), size=logp), shape=4) +
  scale_colour_viridis_c(option = "turbo" , direction=1, begin = 0, end=0.9) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  labs(x="Enriched GO Term", y="Abnormal Chromosome 10", color ="Fold\n Enrich.", size ="-log(p)") +
  ylim(141115174, 195055488) + 
  facet_grid(vars(type), scales = "free_x") +
  theme_dark() +
  theme(legend.position = "bottom", ncol=1, nrow=3, axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title = element_text(size=19), axis.text.x = element_text(size=14), strip.text.y.right = element_text(size = 15)) +
  #####################Ab10 Annotations
  ####TR1 Knob1
  annotate("rect", ymin=142472000, ymax=146699300, xmin=-1, xmax=-2, alpha=0.5, fill="deepskyblue") +
  annotate("rect", ymin=150656000, ymax=153145000, xmin=-1, xmax=-2, alpha=0.5, fill="deepskyblue" ) +
  annotate("rect", ymin=157485200, ymax=159356550, xmin=-1, xmax=-2, alpha=0.5, fill="deepskyblue" ) +
  #####Knob 18
  annotate("rect", ymin=174433450, ymax=182846100, xmin=-1, xmax=-2, alpha=0.5, fill="darkorange3" ) +
  ######Shared Region
  annotate("rect", ymin=141115174, ymax=142472000, xmin=-1, xmax=-2, alpha=0.5, fill="darkgoldenrod1" ) +
  annotate("rect", ymin=152050000, ymax=156350000, xmin=-1, xmax=-2, alpha=0.5, fill="darkgoldenrod1" ) +
  annotate("rect", ymin=158250000, ymax=166820000, xmin=-1, xmax=-2, alpha=0.5, fill="darkgoldenrod1" ) +
  #######trkin
  annotate("rect", ymin=148964528, ymax=149082763, xmin=-1, xmax=-2, fill="blue" ) +
  annotate("rect", ymin=150358622, ymax=150457752, xmin=-1, xmax=-2, fill="blue" ) +
  #######kindr
  annotate("rect", ymin=189326066, ymax=190330226, xmin=-1, xmax=-2, alpha=0.5, fill="hotpink" ) +
  coord_cartesian(clip="off")

p
dev.off()



