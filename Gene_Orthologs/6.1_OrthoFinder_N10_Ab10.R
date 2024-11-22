library(tidyverse)
library(reshape2)
library(splitstackshape)
library(karyoploteR)
library(readxl)
library(ggplot2)
library(pafr)
library(Rsamtools)

setwd("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/Paper/TRKIN_Paper")

ORTHO <- read.delim("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/OrthoFinder/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.Protein.LongestIsoform__v__HiFiAb10.Ab10hapProtein.LongestIsoform.tsv")

B73_GFF <- read.delim("/Volumes/Transcend/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.noheader.gff3", header = FALSE)
colnames(B73_GFF) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

Ab10_GFF <- read.delim("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/OrthoFinder/Ab10_HiFi_v2_corrected.gene.v2.gff3", header = FALSE)
colnames(Ab10_GFF) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")


#This section alters the B73 GFF to be compatible with the OrthoFinder Output
B73_GFF_GENE <- subset(B73_GFF, feature == "gene")

B73_GFF_GENE_temp1 <- separate(B73_GFF_GENE, col=attribute, into= c("ID", "biotype", "logic_name"), sep= ';')

B73_GFF_GENE_temp2 <- separate(B73_GFF_GENE_temp1, col= ID, into=c("Trash", "ID"), sep="=")

B73_GFF_GENE_temp3 <- B73_GFF_GENE_temp2[,c("seqname", "start", "end", "ID")]
colnames(B73_GFF_GENE_temp3) <- c("B73_seqname", "B73_start", "B73_end", "B73_ID")

B73_GFF_GENE <- B73_GFF_GENE_temp3

#This section alters the Ab10 GFF to be compatible with the OrthoFinder Output
Ab10_GFF_GENE <- subset(Ab10_GFF, feature == "gene")

Ab10_GFF_GENE_temp1 <- separate(Ab10_GFF_GENE, col=attribute, into= c("ID", "biotype", "logic_name"), sep= ';')

Ab10_GFF_GENE_temp2 <- separate(Ab10_GFF_GENE_temp1, col= ID, into=c("Trash", "ID"), sep="=")

Ab10_GFF_GENE_temp3 <- Ab10_GFF_GENE_temp2[,c("seqname", "start", "end", "ID")]
colnames(Ab10_GFF_GENE_temp3) <- c("Ab10_seqname", "Ab10_start", "Ab10_end", "Ab10_ID")

Ab10_GFF_GENE_temp3$Ab10_ID <-  gsub("gene:", "", Ab10_GFF_GENE_temp3$Ab10_ID)

Ab10_GFF_GENE <- Ab10_GFF_GENE_temp3

#This section converts the file so that there is only one gene in each file 
ORTHO_MELT_temp1 <- cSplit(ORTHO, "Zm.B73.REFERENCE.NAM.5.0_Zm00001eb.1.Protein.LongestIsoform", sep = ",", direction = "long")

ORTHO_MELT_temp2 <- cSplit(ORTHO_MELT_temp1, "HiFiAb10.Ab10hapProtein.LongestIsoform", sep = ",", direction = "long")

ORTHO_MELT_temp3 <- separate(ORTHO_MELT_temp2, col = Zm.B73.REFERENCE.NAM.5.0_Zm00001eb.1.Protein.LongestIsoform, into = c("B73_ID", "B73_Iso"))

ORTHO_MELT_temp4 <- separate(ORTHO_MELT_temp3, col = HiFiAb10.Ab10hapProtein.LongestIsoform , into = c("Ab10_ID", "Ab10_Iso"))

ORTHO_MELT <- ORTHO_MELT_temp4


#This section merges the GFF files with the OrthoFinder Files
DATA_temp1 <- merge(ORTHO_MELT, B73_GFF_GENE)
DATA_temp2 <- merge(DATA_temp1, Ab10_GFF_GENE)
DATA <- DATA_temp2

#This alters the data for KaryoploteR
GENOME <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/OrthoFinder/Ab10N10_KaryoploteR_genome.csv")

GENOME_Ab10 <- subset(GENOME, Chr == "Ab10" | Chr == "N10")

#This makes the Ab10 and N10 genome
Ab10.genome <- toGRanges(GENOME_Ab10)

#################################################### Orthofinder Plots

#This loads in the links from orthofinder
OrthoLink_10_temp1 <- subset(DATA, DATA$B73_seqname == "chr10" & DATA$Ab10_seqname == "chr10")
#This adds a column to make it compatible with the entap results
OrthoLink_10_temp1$Query.Sequence <- paste(OrthoLink_10_temp1$Ab10_ID, OrthoLink_10_temp1$Ab10_Iso, sep="." )

#This extracts the links between N10 and the Ab10 Specific region of Ab10
Ab10_SPEC<- subset(OrthoLink_10_temp1, B73_start >= 141187279 & ((Ab10_start >= 142472000 & Ab10_end <= 153145000) | (Ab10_start >= 168235374)))

################
#UPDATE THIS WHEN NEW ENTAP IS DONE
################

#This loads in the entap results for the Ab10 annotation
Ab10_entap <- read.delim("~/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/OrthoFinder/Ab10_entap_results.tsv")

#This merges the weird links and the entap results
Ab10_entap_SPEC<- merge(Ab10_SPEC, Ab10_entap, by="Query.Sequence")

#This writes this out
write.csv(Ab10_entap_SPEC, "N10Ab10Specific_Orthologs.csv", row.names = FALSE, quote = FALSE)

#These are the only other genes within the rpd2 repeats that are highlighted by the above file g5850.t1, g5852.t1, g5855.t1, g5857.t1, g5859.t1, g5861.t1 they are unconvincing with only one having very weak similarity to an uncharacterized protein

#This creates the links color
Uninverted <- subset(OrthoLink_10_temp1, Ab10_start < 143000000)
Uninverted_start <- Uninverted[,c("B73_seqname", "B73_start", "B73_end")]
Uninverted_end <- Uninverted[,c("Ab10_seqname", "Ab10_start", "Ab10_end")]
Uninverted_start$B73_seqname <- "N10"
Uninverted_end$Ab10_seqname <- "Ab10"
Uninverted_start_range <- toGRanges(Uninverted_start)
Uninverted_end_range <- toGRanges(Uninverted_end)

Inv1 <- subset(OrthoLink_10_temp1, Ab10_start > 151000000 & Ab10_start < 157960000 )
Inv1_start <- Inv1[,c("B73_seqname", "B73_start", "B73_end")]
Inv1_end <- Inv1[,c("Ab10_seqname", "Ab10_start", "Ab10_end")]
Inv1_start$B73_seqname <- "N10"
Inv1_end$Ab10_seqname <- "Ab10"
Inv1_start_range <- toGRanges(Inv1_start)
Inv1_end_range <- toGRanges(Inv1_end)

Inv2 <- subset(OrthoLink_10_temp1, Ab10_start > 158960000 & Ab10_start < 167968901)
Inv2_start <- Inv2[,c("B73_seqname", "B73_start", "B73_end")]
Inv2_end <- Inv2[,c("Ab10_seqname", "Ab10_start", "Ab10_end")]
Inv2_start$B73_seqname <- "N10"
Inv2_end$Ab10_seqname <- "Ab10"
Inv2_start_range <- toGRanges(Inv2_start)
Inv2_end_range <- toGRanges(Inv2_end)

Uninverted2 <- subset(OrthoLink_10_temp1, Ab10_start >= 167968911 & Ab10_end <= 168235374)
Uninverted2_start <- Uninverted2[,c("B73_seqname", "B73_start", "B73_end")]
Uninverted2_end <- Uninverted2[,c("Ab10_seqname", "Ab10_start", "Ab10_end")]
Uninverted2_start$B73_seqname <- "N10"
Uninverted2_end$Ab10_seqname <- "Ab10"
Uninverted2_start_range <- toGRanges(Uninverted2_start)
Uninverted2_end_range <- toGRanges(Uninverted2_end)

All <- OrthoLink_10_temp1
All_start <- All[,c("B73_seqname", "B73_start", "B73_end")]
All_end <- All[,c("Ab10_seqname", "Ab10_start", "Ab10_end")]
All_start$B73_seqname <- "N10"
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

png(file="N10_Ab10_Orthologes_KaryoploteR.png", width = 480*4, height = 480*3, units = "px", pointsize = 24, bg = "white")
kp <- plotKaryotype(genome = Ab10.genome, chromosomes = c("Ab10", "N10"), plot.type=1, plot.params=pp)
kpRect(kp, chr="N10", x0=141187279, x1=152435371, y0=-0.25, y1=-0.35, col = "darkgoldenrod1")
kpRect(kp, chr="Ab10", x0=141115174, x1=142472000, y0=0.01, y1=0.1, col = "darkgoldenrod1")
kpRect(kp, chr="Ab10", x0=152050000, x1=156880000, y0=0.01, y1=0.1, col = "darkgoldenrod1")
kpRect(kp, chr="Ab10", x0=158250000, x1=167721000, y0=0.01, y1=0.1, col = "darkgoldenrod1")
kpRect(kp, chr="Ab10", x0=142472000, x1=146699300, y0=0.01, y1=0.1, col = "deepskyblue")
kpRect(kp, chr="Ab10", x0=150656000, x1=153145000, y0=0.01, y1=0.1, col = "deepskyblue")
kpRect(kp, chr="Ab10", x0=157485200, x1=159356550, y0=0.01, y1=0.1, col = "deepskyblue")
kpRect(kp, chr="Ab10", x0=174433450, x1=182846100, y0=0.01, y1=0.1, col = "darkorange3")
kpRect(kp, chr="Ab10", x0=148964528, x1=149082763, y0=0.01, y1=0.1, col = "blue", border= "blue")
kpRect(kp, chr="Ab10", x0=150358622, x1=150457752, y0=0.01, y1=0.1, col = "blue", border= "blue")
kpRect(kp, chr="Ab10", x0=189326066, x1=190330226, y0=0.01, y1=0.1, col = "hotpink")
kpPlotLinks(kp, data=All_start_range, data2=All_end_range, r0=-.5, r1 = -.25, y= 1.6, col = "purple" , border = "purple")
kpPlotLinks(kp, data=Uninverted_start_range, data2=Uninverted_end_range, r0=-.5, r1 = -.25, y= 1.6, col = "palegreen", border = "palegreen" )
kpPlotLinks(kp, data=Inv1_start_range, data2=Inv1_end_range, r0=-.5, r1 = -.25, y= 1.6, col = "seagreen1", border = "seagreen1")
kpPlotLinks(kp, data=Inv2_start_range, data2=Inv2_end_range, r0=-.5, r1 = -.25, y= 1.6, col = "aquamarine", border = "aquamarine")
kpPlotLinks(kp, data=Uninverted2_start_range, data2=Uninverted2_end_range, r0=-.5, r1 = -.25, y= 1.6, col = "aquamarine4", border = "aquamarine4" )
dev.off()



