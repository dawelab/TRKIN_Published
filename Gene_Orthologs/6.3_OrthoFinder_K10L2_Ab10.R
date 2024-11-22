library(tidyverse)
library(reshape2)
library(splitstackshape)
library(karyoploteR)
library(readxl)
library(ggplot2)
library(pafr)
library(Rsamtools)

setwd("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/Paper/TRKIN_Paper")

ORTHO <- read.delim("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/OrthoFinder/HiFiAb10.Ab10hapProtein.LongestIsoform__v__CI66_K10L2.K10L2hapProtein.LongestIsoform.tsv")

K10L2_GFF <- read.delim("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/OrthoFinder/CI66_K10L2_v1.gene.v2.gff3", header = FALSE)
colnames(K10L2_GFF) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

Ab10_GFF <- read.delim("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/OrthoFinder/Ab10_HiFi_v2_corrected.gene.v2.gff3", header = FALSE)
colnames(Ab10_GFF) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")


#This section alters the K10L2 GFF to be compatible with the OrthoFinder Output
K10L2_GFF_GENE <- subset(K10L2_GFF, feature == "gene")

K10L2_GFF_GENE_temp1 <- separate(K10L2_GFF_GENE, col=attribute, into= c("ID", "biotype", "logic_name"), sep= ';')

K10L2_GFF_GENE_temp2 <- separate(K10L2_GFF_GENE_temp1, col= ID, into=c("Trash", "ID"), sep="=")

K10L2_GFF_GENE_temp3 <- K10L2_GFF_GENE_temp2[,c("seqname", "start", "end", "ID")]
colnames(K10L2_GFF_GENE_temp3) <- c("K10L2_seqname", "K10L2_start", "K10L2_end", "K10L2_ID")

K10L2_GFF_GENE <- K10L2_GFF_GENE_temp3

#This section alters the Ab10 GFF to be compatible with the OrthoFinder Output
Ab10_GFF_GENE <- subset(Ab10_GFF, feature == "gene")

Ab10_GFF_GENE_temp1 <- separate(Ab10_GFF_GENE, col=attribute, into= c("ID", "biotype", "logic_name"), sep= ';')

Ab10_GFF_GENE_temp2 <- separate(Ab10_GFF_GENE_temp1, col= ID, into=c("Trash", "ID"), sep="=")

Ab10_GFF_GENE_temp3 <- Ab10_GFF_GENE_temp2[,c("seqname", "start", "end", "ID")]
colnames(Ab10_GFF_GENE_temp3) <- c("Ab10_seqname", "Ab10_start", "Ab10_end", "Ab10_ID")

Ab10_GFF_GENE_temp3$Ab10_ID <-  gsub("gene:", "", Ab10_GFF_GENE_temp3$Ab10_ID)

Ab10_GFF_GENE <- Ab10_GFF_GENE_temp3

#Determine which genes might be Potential TEs 
Ab10_GFF_EXON <- subset(Ab10_GFF, feature == "exon")
Ab10_GFF_EXON <- separate(Ab10_GFF_EXON, col = "attribute" , into = c("geneID", "Ab10_Iso"), sep=";")
Ab10_GFF_EXON$geneID <- gsub("ID=", "", Ab10_GFF_EXON$geneID)
Ab10_GFF_EXON <- separate(Ab10_GFF_EXON, col = "geneID" , into = c("geneID", "TRASH"), sep="\\.")
Ab10_GFF_EXON_BED <- Ab10_GFF_EXON[,c("seqname","start", "end", "geneID")]
colnames(Ab10_GFF_EXON_BED) <- c("chr","start", "end", "geneID")
write.table(Ab10_GFF_EXON_BED, "Ab10_GFF_EXON_BED.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')


K10L2_GFF_EXON <- subset(K10L2_GFF, feature == "exon")
K10L2_GFF_EXON <- separate(K10L2_GFF_EXON, col = "attribute" , into = c("geneID", "K10L2_Iso"), sep=";")
K10L2_GFF_EXON$geneID <- gsub("ID=", "", K10L2_GFF_EXON$geneID)
K10L2_GFF_EXON <- separate(K10L2_GFF_EXON, col = "geneID" , into = c("geneID", "TRASH"), sep="\\.")
K10L2_GFF_EXON_BED <- K10L2_GFF_EXON[,c("seqname","start", "end", "geneID")]
colnames(K10L2_GFF_EXON_BED) <- c("chr","start", "end", "geneID")
write.table(K10L2_GFF_EXON_BED, "K10L2_GFF_EXON_BED.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')

#This section downloads the TE annotations and removes any overlapping orthologs
Ab10_TE <- read.delim("~/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/OrthoFinder/Ab10_HiFi_v2_corrected.TE.gff", header=FALSE)
colnames(Ab10_TE) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
Ab10_TE_Ab10 <- subset(Ab10_TE, seqname == "chr10")
Ab10_TE_Ab10_BED <- Ab10_TE_Ab10[,c("seqname", "start", "end")]
Ab10_TE_Ab10_BED <- subset(Ab10_TE_Ab10_BED, seqname == "chr10")
write.table(Ab10_TE_Ab10_BED, "Ab10_TE_Ab10.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')

K10L2_TE <- read.delim("~/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/OrthoFinder/CI66_K10L2_v1.TE.gff", header=FALSE)
colnames(K10L2_TE) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
K10L2_TE_K10L2 <- subset(K10L2_TE, seqname == "K10L2")
K10L2_TE_K10L2_BED <- K10L2_TE_K10L2[,c("seqname", "start", "end")]
K10L2_TE_Ab10_BED <- subset(K10L2_TE_K10L2_BED, seqname == "K10L2")
write.table(K10L2_TE_K10L2_BED, "K10L2_TE_K10L2.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')

#I took these files and ran them with bedtools intersect to find gene exons overlapping with annotated TEs so that I can filter them out
FALSE_GENE_Ab10 <- read.delim("~/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/OrthoFinder/Ab10_TE_Genes.bed", header = FALSE)
FALSE_GENE_K10L2 <- read.delim("~/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/OrthoFinder/K10L2_TE_Genes.bed", header = FALSE)

FALSE_GENE_Ab10 <- unique(FALSE_GENE_Ab10$V4)
FALSE_GENE_K10L2 <- unique(FALSE_GENE_K10L2$V4)

#This section converts the file so that there is only one gene in each line
ORTHO_MELT_temp1 <- cSplit(ORTHO, "CI66_K10L2.K10L2hapProtein.LongestIsoform", sep = ",", direction = "long")

ORTHO_MELT_temp2 <- cSplit(ORTHO_MELT_temp1, "HiFiAb10.Ab10hapProtein.LongestIsoform", sep = ",", direction = "long")

ORTHO_MELT_temp3 <- separate(ORTHO_MELT_temp2, col = CI66_K10L2.K10L2hapProtein.LongestIsoform, into = c("K10L2_ID", "K10L2_Iso"))

ORTHO_MELT_temp4 <- separate(ORTHO_MELT_temp3, col = HiFiAb10.Ab10hapProtein.LongestIsoform , into = c("Ab10_ID", "Ab10_Iso"))

ORTHO_MELT <- ORTHO_MELT_temp4

#Drop any genes with exons that overlap annotated TEs
ORTHO_MELT_filt1 <- ORTHO_MELT[-c(ORTHO_MELT$Ab10_ID %in% FALSE_GENE_Ab10),]
ORTHO_MELT_filt1 <- ORTHO_MELT[-c(ORTHO_MELT$K10L2_ID %in% FALSE_GENE_K10L2),]

#This section merges the GFF files with the OrthoFinder Files
DATA_temp1 <- merge(ORTHO_MELT, K10L2_GFF_GENE)
DATA_temp2 <- merge(DATA_temp1, Ab10_GFF_GENE)
DATA <- DATA_temp2

#This alters the data for KaryoploteR
GENOME <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/OrthoFinder/Ab10K10L2_KaryoploteR_genome.csv")

GENOME_Ab10 <- subset(GENOME, Chr == "Ab10" | Chr == "K10L2")

#This makes the Ab10 and K10L2 genome
Ab10.genome <- toGRanges(GENOME_Ab10)

#################################################### Orthofinder Plots

#This loads in the links from orthofinder
OrthoLink_10_temp1 <- subset(DATA, DATA$K10L2_seqname == "K10L2" & DATA$Ab10_seqname == "chr10" & DATA$K10L2_start >= 2730186 & DATA$Ab10_start >= 141115174)

#This adds a column to make it compatible with the entap results
OrthoLink_10_temp1$Query.Sequence <- paste(OrthoLink_10_temp1$Ab10_ID, OrthoLink_10_temp1$Ab10_Iso, sep="." )

#This extracts the links between the specific regions of K10L2 and Ab10
SPEC <- subset(OrthoLink_10_temp1, (K10L2_start >= 7429619 & K10L2_end <= 25502283) & ((Ab10_start >= 142472000 & Ab10_end <= 153145000) | (Ab10_start >= 168235374)))

################
#UPDATE THIS WHEN NEW ENTAP IS DONE
################

#This loads in the entap results for the Ab10 annotation
Ab10_entap <- read.delim("~/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/OrthoFinder/Ab10_entap_results.tsv")

#This merges the SGE specific links and the entap results
Ab10_entap_SPEC <- merge(SPEC, Ab10_entap, by="Query.Sequence")

#This writes out the file
write.csv(Ab10_entap_SPEC, "K10L2Ab10Specific_Orthologs.csv", row.names = FALSE, quote = FALSE)

#This defines the known regions for plotting 
Uninverted <- subset(OrthoLink_10_temp1, Ab10_start < 143000000)
Uninverted_start <- Uninverted[,c("K10L2_seqname", "K10L2_start", "K10L2_end")]
Uninverted_end <- Uninverted[,c("Ab10_seqname", "Ab10_start", "Ab10_end")]
Uninverted_start$K10L2_seqname <- "K10L2"
Uninverted_end$Ab10_seqname <- "Ab10"
Uninverted_start_range <- toGRanges(Uninverted_start)
Uninverted_end_range <- toGRanges(Uninverted_end)

Inv1 <- subset(OrthoLink_10_temp1, Ab10_start > 151000000 & Ab10_start < 157960000 )
Inv1_start <- Inv1[,c("K10L2_seqname", "K10L2_start", "K10L2_end")]
Inv1_end <- Inv1[,c("Ab10_seqname", "Ab10_start", "Ab10_end")]
Inv1_start$K10L2_seqname <- "K10L2"
Inv1_end$Ab10_seqname <- "Ab10"
Inv1_start_range <- toGRanges(Inv1_start)
Inv1_end_range <- toGRanges(Inv1_end)

Inv2 <- subset(OrthoLink_10_temp1, Ab10_start > 158960000 & Ab10_start < 167968901)
Inv2_start <- Inv2[,c("K10L2_seqname", "K10L2_start", "K10L2_end")]
Inv2_end <- Inv2[,c("Ab10_seqname", "Ab10_start", "Ab10_end")]
Inv2_start$K10L2_seqname <- "K10L2"
Inv2_end$Ab10_seqname <- "Ab10"
Inv2_start_range <- toGRanges(Inv2_start)
Inv2_end_range <- toGRanges(Inv2_end)

Uninverted2 <- subset(OrthoLink_10_temp1, Ab10_start >= 167968911 & Ab10_end <= 168235374)
Uninverted2_start <- Uninverted2[,c("K10L2_seqname", "K10L2_start", "K10L2_end")]
Uninverted2_end <- Uninverted2[,c("Ab10_seqname", "Ab10_start", "Ab10_end")]
Uninverted2_start$K10L2_seqname <- "K10L2"
Uninverted2_end$Ab10_seqname <- "Ab10"
Uninverted2_start_range <- toGRanges(Uninverted2_start)
Uninverted2_end_range <- toGRanges(Uninverted2_end)

All <- OrthoLink_10_temp1
All_start <- All[,c("K10L2_seqname", "K10L2_start", "K10L2_end")]
All_end <- All[,c("Ab10_seqname", "Ab10_start", "Ab10_end")]
All_start$K10L2_seqname <- "K10L2"
All_end$Ab10_seqname <- "Ab10"
All_start_range <- toGRanges(All_start)
All_end_range <- toGRanges(All_end)

#This alters plot parameters, this is actually just the default still
pp <- getDefaultPlotParams(plot.type=1)
pp$ideogramheight <- 1
pp$data1height <- 200
pp$data2height <- 100
pp$bottommargin <- 200

#This does the plotting
png(file="K10L2_Ab10_Orthologes_KaryoploteR.png", width = 480*4, height = 480*3, units = "px", pointsize = 24, bg = "white")
kp <- plotKaryotype(genome = Ab10.genome, chromosomes = c("K10L2","Ab10"), plot.type=1, plot.params=pp)
kpRect(kp, chr="K10L2", x0=7429619, x1=13957822, y0=0.01, y1=0.1, col = "deepskyblue")
kpRect(kp, chr="K10L2", x0=15726572, x1=15898240, y0=0.01, y1=0.1, col = "deepskyblue")
kpRect(kp, chr="K10L2", x0=16787371, x1=25024178, y0=0.01, y1=0.1, col = "deepskyblue")
kpRect(kp, chr="K10L2", x0=25498094, x1=25502283, y0=0.01, y1=0.1, col = "deepskyblue")
kpRect(kp, chr="K10L2", x0=16328801, x1=16357145, y0=0.01, y1=0.1, col = "blue")
kpRect(kp, chr="K10L2", x0=2729416, x1=7429619, y0=0.01, y1=0.1, col = "darkgoldenrod1")
kpRect(kp, chr="K10L2", x0=25502283, x1=31891546, y0=0.01, y1=0.1, col = "darkgoldenrod1")
kpRect(kp, chr="Ab10", x0=141115174, x1=142472000, y0=-0.25, y1=-0.35, col = "darkgoldenrod1")
kpRect(kp, chr="Ab10", x0=152050000, x1=156880000, y0=-0.25, y1=-0.35, col = "darkgoldenrod1")
kpRect(kp, chr="Ab10", x0=158250000, x1=167721000, y0=-0.25, y1=-0.35, col = "darkgoldenrod1")
kpRect(kp, chr="Ab10", x0=142472000, x1=146699300, y0=-0.25, y1=-0.35, col = "deepskyblue")
kpRect(kp, chr="Ab10", x0=150656000, x1=153145000, y0=-0.25, y1=-0.35, col = "deepskyblue")
kpRect(kp, chr="Ab10", x0=157485200, x1=159356550, y0=-0.25, y1=-0.35, col = "deepskyblue")
kpRect(kp, chr="Ab10", x0=174433450, x1=182846100, y0=-0.25, y1=-0.35, col = "darkorange3")
kpRect(kp, chr="Ab10", x0=148964528, x1=149082763, y0=-0.25, y1=-0.35, col = "blue", border= "blue")
kpRect(kp, chr="Ab10", x0=150358622, x1=150457752, y0=-0.25, y1=-0.35, col = "blue", border= "blue")
kpRect(kp, chr="Ab10", x0=189326066, x1=190330226, y0=-0.25, y1=-0.35, col = "hotpink")
kpPlotLinks(kp, data=All_start_range, data2=All_end_range, r0=-.5, r1 = -.25, y= 1.6, col = "purple" , border = "purple")
kpPlotLinks(kp, data=Uninverted_start_range, data2=Uninverted_end_range, r0=-.5, r1 = -.25, y= 1.6, col = "palegreen", border = "palegreen" )
kpPlotLinks(kp, data=Inv1_start_range, data2=Inv1_end_range, r0=-.5, r1 = -.25, y= 1.6, col = "seagreen1", border = "seagreen1")
kpPlotLinks(kp, data=Inv2_start_range, data2=Inv2_end_range, r0=-.5, r1 = -.25, y= 1.6, col = "aquamarine", border = "aquamarine")
kpPlotLinks(kp, data=Uninverted2_start_range, data2=Uninverted2_end_range, r0=-.5, r1 = -.25, y= 1.6, col = "aquamarine4", border = "aquamarine4" )
dev.off()

