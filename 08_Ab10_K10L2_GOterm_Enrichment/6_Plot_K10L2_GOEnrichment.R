library(topGO)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(grid)
library(egg)

setwd("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/Paper/TRKIN_Published/GOterm_Enrichment")

K10L2_GFF <- read.delim("~/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/OrthoFinder/CI66_K10L2_v1.gene.v2.gff3", header = FALSE)
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
DATA_MF <- new("topGOdata", description="K10L2", ontology = "MF", allGenes = geneList, nodeSize = 10, annot = annFUN.file, file = "K10L2_GOMapping.tbl")
DATA_CC <- new("topGOdata", description="K10L2", ontology = "CC", allGenes = geneList, nodeSize = 10, annot = annFUN.file, file = "K10L2_GOMapping.tbl")

#This loads in the 10 most enriched GO terms for each category
ENR <- read.csv("3.2_K10L2_GOTermEnrichment.csv")
ENR$Type <- gsub("Biological Process", "DATA_BP", ENR$Type)
ENR$Type <- gsub("Cellular Component", "DATA_CC", ENR$Type)
ENR$Type <- gsub("Molecular Function", "DATA_MF", ENR$Type)

#This filters anything with a p below 0.01
ENR <- subset(ENR, classicFisher <= 0.01)

#This loads in the entap results 
K10L2_entap <- read.delim("~/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/GOTermEnrichment/K10L2_entap_results.tsv")
K10L2_entap$GO_ALL <- paste(K10L2_entap$EggNOG.GO.Biological, K10L2_entap$EggNOG.GO.Cellular, K10L2_entap$EggNOG.GO.Molecular, sep = ",")

#This merges the entap results and GFF file
MERGE <- merge(K10L2_SPEC, K10L2_entap, by.x="biotype", by.y="Query.Sequence")

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
  if(DATA=="DATA_BP"){MERGE_GO$type ="Biological\n Process"}
  if(DATA=="DATA_MF"){MERGE_GO$type ="Molecular\n Function"}
  if(DATA=="DATA_CC"){MERGE_GO$type ="Cellular\n Component"}
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


pdf("K10L2_GOterm_Enrichment.pdf", height=10, width=5)
p <- ggplot() +
  geom_point(data=MERGE_GO_ALL, aes(y=start, x=term, color=as.numeric(enrich), size=logp), shape=4) +
  scale_colour_viridis_c(option = "turbo" , direction=1, breaks = c(24.4, 26.19), labels = c("24","26"), begin = 0, end=0.9) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  labs(x="Enriched GO Term", y= "K10L2", color = "Fold\n Enrich.", size = "-log(p)") +
  ylim(2730186, 31891546) + 
  facet_grid(vars(type), scales = "free_x") +
  theme_dark() +
  theme(legend.position = "right", ncol=1, nrow=3, axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title = element_text(size=19), axis.text = element_text(size=14), strip.text.x.right = element_text(size = 15)) +
  #####################K10L2 Annotations
  ####TR1 Knob1
  annotate("rect", ymin=7429619, ymax=13957822, xmin=-1, xmax=-2, alpha=0.5, fill="deepskyblue") +
  annotate("rect", ymin=15726572, ymax=15898240, xmin=-1, xmax=-2, alpha=0.5, fill="deepskyblue" ) +
  annotate("rect", ymin=16787371, ymax=25024178, xmin=-1, xmax=-2, alpha=0.5, fill="deepskyblue" ) +
  annotate("rect", ymin=25498094, ymax=25502283, xmin=-1, xmax=-2, alpha=0.5, fill="deepskyblue" ) +
  ######Shared Region
  annotate("rect", ymin=2730186, ymax=7429619, xmin=-1, xmax=-2, alpha=0.5, fill="darkgoldenrod1" ) +
  annotate("rect", ymin=25502283, ymax=31891546, xmin=-1, xmax=-2, alpha=0.5, fill="darkgoldenrod1" ) +
  #######trkin
  annotate("rect", ymin=16328801, ymax=16357145, xmin=-1, xmax=-2, fill="blue" ) +
  coord_cartesian(clip="off")

p
dev.off()

