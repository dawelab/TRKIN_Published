library(topGO)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(grid)
library(egg)

setwd("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/Paper/TRKIN_Paper")


Ab10_GFF <- read.delim("~/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/OrthoFinder/HiFiAb10.genes.edit.gff3", header = FALSE)
#This drops an unnecessary column at the top
Ab10_GFF <- Ab10_GFF[-c(1),]
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
Ab10_SPEC<- subset(Ab10_GFF_GENE_f2, start >= 141187279 & ((start >= 142472000 & end <= 153145000) | (start >= 168235374)))

#This drops 8 copies of the 9 copies of nrpd2/e2 as they seem to be skewing the analysis 
Ab10_SPEC <- subset(Ab10_SPEC, biotype != "g5851.t1" & biotype != "g5853.t1" & biotype != "g5854.t1" & biotype != "g5856.t1" & biotype != "g5858.t1" & biotype != "g5860.t1" & biotype != "g5862.t1" & biotype != "g5863.t1" & biotype != "g5880.t1" & biotype != "g5879.t1" & biotype != "g5874.t1")

#This drops the 9 copies of 

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
DATA_MF <- new("topGOdata", description="Ab10", ontology = "MF", allGenes = geneList, nodeSize = 10, annot = annFUN.file, file = "Ab10_GOMapping.tbl")
DATA_CC <- new("topGOdata", description="Ab10", ontology = "CC", allGenes = geneList, nodeSize = 10, annot = annFUN.file, file = "Ab10_GOMapping.tbl")

#This loads in the 10 most enriched GO terms for each category
ENR <- read.csv("Ab10_GOTermEnrichment_NoDup.csv")
ENR$Type <- gsub("Biological Process", "DATA_BP", ENR$Type)
ENR$Type <- gsub("Cellular Component", "DATA_CC", ENR$Type)
ENR$Type <- gsub("Molecular Function", "DATA_MF", ENR$Type)

#This loads in the entap results 
Ab10_entap <- read.delim("~/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/GOTermEnrichment/Ab10_entap_results.tsv")
Ab10_entap$GO_ALL <- paste(Ab10_entap$EggNOG.GO.Biological, Ab10_entap$EggNOG.GO.Cellular, Ab10_entap$EggNOG.GO.Molecular, sep = ",")

#This merges the entap results and GFF file
MERGE <- merge(Ab10_SPEC, Ab10_entap, by.x="biotype", by.y="Query.Sequence")

MERGE_GO_ALL <- data.frame("chr"=c(), "start"=c(), "end"=c(), "term"=c(), "enrich"=c(), "p"=c(), "logp"=c(), "type"=c() )


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

pdf("Ab10_GOterm_Enrichment_NoDup.pdf", height=7, width=10)
p <- ggplot() +
  geom_point(data=MERGE_GO_ALL, aes(x=start, y=term, color=as.numeric(enrich), size=logp), shape=4) +
  scale_colour_viridis_c(option = "turbo" , direction=1, breaks = c(255.6, 281), labels = c("255","281"), begin = 0, end=0.9) +
  labs(y="Enriched GO Term", x= "Abnormal Chromosome 10", color = "Fold Enrichment", size = "-log(p)") +
  xlim(141115174, 195055488) + 
  facet_grid(vars(type), scales = "free_y") +
  theme(legend.position = "bottom", ncol=1, nrow=3, axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title = element_text(size=20), axis.text = element_text(size=10), strip.text.y.right = element_text(size = 15)) +
  #####################Ab10 Annotations
  ####TR1 Knob1
  annotate("rect", xmin=142472000, xmax=146699300, ymin=-1, ymax=-2, alpha=0.5, fill="deepskyblue") +
  annotate("rect", xmin=150656000, xmax=153145000, ymin=-1, ymax=-2, alpha=0.5, fill="deepskyblue" ) +
  annotate("rect", xmin=157485200, xmax=159356550, ymin=-1, ymax=-2, alpha=0.5, fill="deepskyblue" ) +
  #####Knob 180
  annotate("rect", xmin=174433450, xmax=182846100, ymin=-1, ymax=-2, alpha=0.5, fill="darkorange3" ) +
  ######Shared Region
  annotate("rect", xmin=141115174, xmax=142472000, ymin=-1, ymax=-2, alpha=0.5, fill="darkgoldenrod1" ) +
  annotate("rect", xmin=152050000, xmax=156350000, ymin=-1, ymax=-2, alpha=0.5, fill="darkgoldenrod1" ) +
  annotate("rect", xmin=158250000, xmax=166820000, ymin=-1, ymax=-2, alpha=0.5, fill="darkgoldenrod1" ) +
  #######trkin
  annotate("rect", xmin=148964528, xmax=149082763, ymin=-1, ymax=-2, fill="blue" ) +
  annotate("rect", xmin=150358622, xmax=150457752, ymin=-1, ymax=-2, fill="blue" ) +
  #######kindr
  annotate("rect", xmin=189326066, xmax=190330226, ymin=-1, ymax=-2, alpha=0.5, fill="hotpink" ) +
  coord_cartesian(clip="off")

p
dev.off()

