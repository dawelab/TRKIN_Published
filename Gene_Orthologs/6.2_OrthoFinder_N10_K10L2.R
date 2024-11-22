library(tidyverse)
library(reshape2)
library(splitstackshape)
library(karyoploteR)
library(readxl)
library(ggplot2)
library(pafr)
library(Rsamtools)

setwd("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/Paper/TRKIN_Paper")

ORTHO <- read.delim("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/OrthoFinder/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.Protein.LongestIsoform__v__CI66_K10L2.K10L2hapProtein.LongestIsoform.tsv")

B73_GFF <- read.delim("/Volumes/Transcend/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.noheader.gff3", header = FALSE)
colnames(B73_GFF) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

K10L2_GFF <- read.delim("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/OrthoFinder/CI66_K10L2_v1.gene.v2.gff3", header = FALSE)
colnames(K10L2_GFF) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")


#This section alters the B73 GFF to be compatible with the OrthoFinder Output
B73_GFF_GENE <- subset(B73_GFF, feature == "gene")

B73_GFF_GENE_temp1 <- separate(B73_GFF_GENE, col=attribute, into= c("ID", "biotype", "logic_name"), sep= ';')

B73_GFF_GENE_temp2 <- separate(B73_GFF_GENE_temp1, col= ID, into=c("Trash", "ID"), sep="=")

B73_GFF_GENE_temp3 <- B73_GFF_GENE_temp2[,c("seqname", "start", "end", "ID")]
colnames(B73_GFF_GENE_temp3) <- c("B73_seqname", "B73_start", "B73_end", "B73_ID")

B73_GFF_GENE <- B73_GFF_GENE_temp3

#This section alters the K10L2 GFF to be compatible with the OrthoFinder Output
K10L2_GFF_GENE <- subset(K10L2_GFF, feature == "gene")

K10L2_GFF_GENE_temp1 <- separate(K10L2_GFF_GENE, col=attribute, into= c("ID", "biotype", "logic_name"), sep= ';')

K10L2_GFF_GENE_temp2 <- separate(K10L2_GFF_GENE_temp1, col= ID, into=c("Trash", "ID"), sep="=")

K10L2_GFF_GENE_temp3 <- K10L2_GFF_GENE_temp2[,c("seqname", "start", "end", "ID")]
colnames(K10L2_GFF_GENE_temp3) <- c("K10L2_seqname", "K10L2_start", "K10L2_end", "K10L2_ID")

K10L2_GFF_GENE_temp3$K10L2_ID <-  gsub("gene:", "", K10L2_GFF_GENE_temp3$K10L2_ID)

K10L2_GFF_GENE <- K10L2_GFF_GENE_temp3

#This section converts the file so that there is only one gene in each file 
ORTHO_MELT_temp1 <- cSplit(ORTHO, "Zm.B73.REFERENCE.NAM.5.0_Zm00001eb.1.Protein.LongestIsoform", sep = ",", direction = "long")

ORTHO_MELT_temp2 <- cSplit(ORTHO_MELT_temp1, "CI66_K10L2.K10L2hapProtein.LongestIsoform", sep = ",", direction = "long")

ORTHO_MELT_temp3 <- separate(ORTHO_MELT_temp2, col = Zm.B73.REFERENCE.NAM.5.0_Zm00001eb.1.Protein.LongestIsoform, into = c("B73_ID", "B73_Iso"))

ORTHO_MELT_temp4 <- separate(ORTHO_MELT_temp3, col = CI66_K10L2.K10L2hapProtein.LongestIsoform , into = c("K10L2_ID", "K10L2_Iso"))

ORTHO_MELT <- ORTHO_MELT_temp4

#This section merges the GFF files with the OrthoFinder Files
DATA_temp1 <- merge(ORTHO_MELT, B73_GFF_GENE)
DATA_temp2 <- merge(DATA_temp1, K10L2_GFF_GENE)
DATA <- DATA_temp2

#This alters the data for KaryoploteR
GENOME <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/OrthoFinder/K10L2N10_KaryoploteR_genome.csv")

GENOME_K10L2 <- subset(GENOME, Chr == "K10L2" | Chr == "N10")

#This makes the K10L2 and N10 genome
K10L2.genome <- toGRanges(GENOME_K10L2)

#################################################### Orthofinder Plots

#This loads in the links from orthofinder
OrthoLink_10_temp1 <- subset(DATA, DATA$B73_seqname == "chr10" & DATA$K10L2_seqname == "K10L2")

All <- OrthoLink_10_temp1
All_start <- All[,c("B73_seqname", "B73_start", "B73_end")]
All_end <- All[,c("K10L2_seqname", "K10L2_start", "K10L2_end")]
All_start$B73_seqname <- "N10"
All_end$K10L2_seqname <- "K10L2"
All_start_range <- toGRanges(All_start)
All_end_range <- toGRanges(All_end)

Uninverted <- subset(OrthoLink_10_temp1, B73_start > 141210513-1000 & B73_end < 142645291+1000)
Uninverted_start <- Uninverted[,c("K10L2_seqname", "K10L2_start", "K10L2_end")]
Uninverted_end <- Uninverted[,c("B73_seqname", "B73_start", "B73_end")]
Uninverted_start$K10L2_seqname <- "K10L2"
Uninverted_end$B73_seqname <- "N10"
Uninverted_start_range <- toGRanges(Uninverted_start)
Uninverted_end_range <- toGRanges(Uninverted_end)

Inv1 <- subset(OrthoLink_10_temp1, B73_start > 142645291 & B73_start < 145909912)
Inv1_start <- Inv1[,c("K10L2_seqname", "K10L2_start", "K10L2_end")]
Inv1_end <- Inv1[,c("B73_seqname", "B73_start", "B73_end")]
Inv1_start$K10L2_seqname <- "K10L2"
Inv1_end$B73_seqname <- "N10"
Inv1_start_range <- toGRanges(Inv1_start)
Inv1_end_range <- toGRanges(Inv1_end)

Inv2 <- subset(OrthoLink_10_temp1, B73_start > 145909912 & B73_start < 152287382)
Inv2_start <- Inv2[,c("K10L2_seqname", "K10L2_start", "K10L2_end")]
Inv2_end <- Inv2[,c("B73_seqname", "B73_start", "B73_end")]
Inv2_start$K10L2_seqname <- "K10L2"
Inv2_end$B73_seqname <- "N10"
Inv2_start_range <- toGRanges(Inv2_start)
Inv2_end_range <- toGRanges(Inv2_end)

Uninverted2 <- subset(OrthoLink_10_temp1, B73_start >= 151806944 & B73_end <= 151949983)
Uninverted2_start <- Uninverted2[,c("K10L2_seqname", "K10L2_start", "K10L2_end")]
Uninverted2_end <- Uninverted2[,c("B73_seqname", "B73_start", "B73_end")]
Uninverted2_start$K10L2_seqname <- "K10L2"
Uninverted2_end$B73_seqname <- "N10"
Uninverted2_start_range <- toGRanges(Uninverted2_start)
Uninverted2_end_range <- toGRanges(Uninverted2_end)

#This alters plot parameters, this is actually just the default still
pp <- getDefaultPlotParams(plot.type=1)
pp$ideogramheight <- 1
pp$data1height <- 200
pp$data2height <- 100
pp$bottommargin <- 200

#This does the plotting
png(file="N10_K10L2_Orthologes_KaryoploteR.png", width = 480*4, height = 480*3, units = "px", pointsize = 24, bg = "white")
kp <- plotKaryotype(genome = K10L2.genome, chromosomes = c("N10", "K10L2"), plot.type=1, plot.params=pp)
kpRect(kp, chr="N10", x0=141187279, x1=152435371, y0=0.01, y1=0.1, col = "darkgoldenrod1")
kpRect(kp, chr="K10L2", x0=7429619, x1=13957822, y0=-0.25, y1=-0.35, col = "deepskyblue")
kpRect(kp, chr="K10L2", x0=15726572, x1=15898240, y0=-0.25, y1=-0.35, col = "deepskyblue")
kpRect(kp, chr="K10L2", x0=16787371, x1=25024178, y0=-0.25, y1=-0.35, col = "deepskyblue")
kpRect(kp, chr="K10L2", x0=25498094, x1=25502283, y0=-0.25, y1=-0.35, col = "deepskyblue")
kpRect(kp, chr="K10L2", x0=16328801, x1=16357145, y0=-0.25, y1=-0.35, col = "blue")
kpRect(kp, chr="K10L2", x0=2729416, x1=7429619, y0=-0.25, y1=-0.35, col = "darkgoldenrod1")
kpRect(kp, chr="K10L2", x0=25502283, x1=31891546, y0=-0.25, y1=-0.35, col = "darkgoldenrod1")
kpPlotLinks(kp, data=Uninverted_start_range, data2=Uninverted_end_range, r0=-.5, r1 = -.25, y= 1.6, col = "palegreen", border = "palegreen" )
kpPlotLinks(kp, data=Inv1_start_range, data2=Inv1_end_range, r0=-.5, r1 = -.25, y= 1.6, col = "seagreen1", border = "seagreen1")
kpPlotLinks(kp, data=Inv2_start_range, data2=Inv2_end_range, r0=-.5, r1 = -.25, y= 1.6, col = "aquamarine", border = "aquamarine")
kpPlotLinks(kp, data=Uninverted2_start_range, data2=Uninverted2_end_range, r0=-.5, r1 = -.25, y= 1.6, col = "aquamarine4", border = "aquamarine4" )
dev.off()

