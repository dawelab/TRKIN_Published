library(ggplot2)
library(ggpubr)
library(readxl)
library(stringr)
library(tidyverse)

setwd("")

#I manually merged the TE and gene trkin region annotation files for Ab10trkin1, Ab10trkin2, and K10L2trkin edited these and isolated the TE superfamily in excel

#Load the All feautes file with both TE and gene annotations
Ab10trkin1_AllFeatures <- read.csv("Ab10trkin1_AllFeatures.csv")
#Set the start of the gene to 0 for easier plotting
Ab10trkin1_AllFeatures$start <- Ab10trkin1_AllFeatures$start-148968414
Ab10trkin1_AllFeatures$end <- Ab10trkin1_AllFeatures$end-148968414
#Set the y coordinate for plotting
Ab10trkin1_AllFeatures$y <- 0

#Separate genes and TE annotations for plotting
Ab10trkin1_AllFeatures_TE <- subset(Ab10trkin1_AllFeatures, repeattype == "dispersed_repeat")
Ab10trkin1_AllFeatures_Ex <- subset(Ab10trkin1_AllFeatures, repeattype != "dispersed_repeat")

#Load the All feautes file with both TE and gene annotations
Ab10trkin2_AllFeatures <- read.csv("Ab10trkin2_AllFeatures.csv")
#Set the start of the gene to 0 for easier plotting
Ab10trkin2_AllFeatures$start <- Ab10trkin2_AllFeatures$start-150358418
Ab10trkin2_AllFeatures$end <- Ab10trkin2_AllFeatures$end-150358418
#Set the y coordinate for plotting
Ab10trkin2_AllFeatures$y <- 1

#Separate genes and TE annotations for plotting
Ab10trkin2_AllFeatures_TE <- subset(Ab10trkin2_AllFeatures, repeattype == "dispersed_repeat")
Ab10trkin2_AllFeatures_Ex <- subset(Ab10trkin2_AllFeatures, repeattype != "dispersed_repeat")

#Load the All feautes file with both TE and gene annotations
K10L2trkin_AllFeatures <- read.csv("K10L2trkin_AllFeatures.csv")
#Set the start of the gene to 0 for easier plotting
K10L2trkin_AllFeatures$start <- K10L2trkin_AllFeatures$start-16280429
K10L2trkin_AllFeatures$end <- K10L2trkin_AllFeatures$end-16280429
#Set the y coordinate for plotting
K10L2trkin_AllFeatures$y <- 2

#Separate genes and TE annotations for plotting
K10L2trkin_AllFeatures_TE <- subset(K10L2trkin_AllFeatures, repeattype == "dispersed_repeat")
K10L2trkin_AllFeatures_Ex <- subset(K10L2trkin_AllFeatures, repeattype != "dispersed_repeat")

#Merge all 3 trkins exon files for plotting
trkin_all_EX <- rbind(Ab10trkin1_AllFeatures_Ex, Ab10trkin2_AllFeatures_Ex, K10L2trkin_AllFeatures_Ex)
#Separate out the annotation descriptions to isoate the part of the annotation (UTR or not) and rename them for better plotting. 
trkin_all_EX <- separate(trkin_all_EX, col = ID, into = c("ID", "Parent"), sep = ";")
trkin_all_EX$ID <- gsub("ID=", "", trkin_all_EX$ID)
trkin_all_EX <- separate(trkin_all_EX, col = ID, into = c("ID", "transcript", "Exon"), sep = "[.]")
trkin_all_EX$Exon <- gsub("exon", "", trkin_all_EX$Exon)
trkin_all_EX$Exon <- gsub("utr5p1", "5' UTR", trkin_all_EX$Exon)
trkin_all_EX$Exon <- gsub("utr3p1", "3' UTR", trkin_all_EX$Exon)
#Set the y coordinates for each exonic features label 
trkin_all_EX$y2 <- c(-0, -0.1, -0.2, 0, 0, -0.1, 0, -0.1, 0, -0.1, 0, -0.1, 0, -0.1, 0, -0.1, 0, -0.1, 0, -0.1, -0.2, 1, 0.9, 0.8, 1, 0.9, 1, 0.9, 1, 0.9, 1, 0.9, 1, 0.9, 1, 0.9, 1, 0.9, 1, 1, 0.9, 0.8, 2, 1.9, 1.8, 2, 1.9, 2, 1.9, 2, 1.9, 2, 1.9, 2, 1.9, 2, 1.9, 2, 1.9, 2, 2, 1.9, 1.8)

#Write out the processed file for plotting
write.csv(trkin_all_EX, "trkin_all_EX.csv", row.names = FALSE, quote = FALSE)

#Merge all the trkin TE files for plotting
trkin_all_TE <- rbind(Ab10trkin1_AllFeatures_TE, Ab10trkin2_AllFeatures_TE, K10L2trkin_AllFeatures_TE)

#This redefines the L1 and LINE unknown to just be LINE for cleaner plotting
trkin_all_TE$Superfamily <- gsub("LINE_unknown", "LINE", trkin_all_TE$Superfamily)
trkin_all_TE$Superfamily <- gsub("L1", "LINE", trkin_all_TE$Superfamily)
trkin_all_TE$Superfamily <- gsub("DTA", "hAT", trkin_all_TE$Superfamily)
trkin_all_TE$Superfamily <- gsub("LTR_unknown", "unknown LTR", trkin_all_TE$Superfamily)
trkin_all_TE$Superfamily <- gsub("DTH", "Harbinger", trkin_all_TE$Superfamily)

# This reorders the superfamily factor for better plotting
trkin_all_TE$Superfamily <- factor(trkin_all_TE$Superfamily, levels=c("Gypsy", "LINE", "Copia", "Harbinger",  "unknown LTR",  "hAT"))
summary(as.factor(trkin_all_TE$Superfamily))

#This removes any simple repeats for better plotting
trkin_all_TE <- subset(trkin_all_TE, Superfamily != "Simple Repeat")

#Read in the trkin specific TEs 
SPEC_TE <- read.csv("trkin_specific_TE.csv", header = TRUE)

#This creates the annotations
Ab10trkin1grob <- text_grob("Ab10 trkin1", face = "bold", color = "black", rot = 90)
Ab10trkin2grob <- text_grob("Ab10 trkin2", face = "bold", color = "black", rot = 90)
K10L2trkingrob <- text_grob("K10L2 trkin", face = "bold", color = "black", rot = 90)

#This cleaves only the first description, which has superfamily for plotting
trkin_all_TE$MotifSimp <- lapply(trkin_all_TE$Motif , function(x) {strsplit(x, "_")[[1]][1]})
trkin_all_TE$MotifSimp <- as.character(trkin_all_TE$MotifSimp)

#This generates the plot
pdf("trkins.pdf", height = 8, width = 6)
ggplot() +
  theme_classic() +
  geom_segment(data=trkin_all_TE, size = 10,  aes(x = start, y = y, xend = end, yend = y, colour = Superfamily)) +
  scale_colour_brewer(palette = "Set2") +
  geom_segment(data=trkin_all_EX, size = 7, aes(x = start, y = y+.2, xend = end, yend = y+.2), color = "black") +
  geom_segment(aes(x = 0, y = 0.2, xend = 111420, yend = 0.2), arrow = arrow(length = unit(0.5, "cm"))) +
  geom_segment(aes(x = 99521, y = 1.2, xend = 0, yend = 1.2), arrow = arrow(length = unit(0.5, "cm"))) +
  geom_segment(aes(x = 76901, y = 2.2, xend = 0, yend = 2.2), arrow = arrow(length = unit(0.5, "cm"))) +
  geom_segment(data=SPEC_TE, aes(x=x, y=y-0.01, xend=xend, yend=yend-0.01), size = 4, color = "navy") +
  geom_text(data=trkin_all_EX, aes(x=start, y= y2+0.5, label=Exon)) +
  ylab("") +
  xlab("") +
  xlim(-5000, 120000) +
  scale_x_continuous(breaks=seq(0,120000,10000), labels = ~ format(.x, scientific = FALSE)) +
  annotation_custom(Ab10trkin1grob, xmin=-5000, xmax=-5000,  ymin=0, ymax=0) +
  annotation_custom(Ab10trkin2grob, xmin=-5000, xmax=-5000,  ymin=1, ymax=1) +
  annotation_custom(K10L2trkingrob, xmin=-5000, xmax=-5000,  ymin=2, ymax=2) +
  theme(legend.position = "bottom", axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_text(angle = 45, vjust=0.5, size = 8)) +
  guides(color=guide_legend(ncol=3, nrow=2)) +
  coord_cartesian(clip="off")
dev.off()




