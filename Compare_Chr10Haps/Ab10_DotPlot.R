###This section is from online and parses the .delta file
library(dplyr)
library(magrittr)
library(GenomicRanges)
library(knitr)
library(ggplot2)
library(tidyr)
library(ggbio)
library(ggpubr)
library(RColorBrewer)
options(ucscChromosomeNames=FALSE)

#This loads the coords
coords <- read.delim("~/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/K10L2_DotPlot/Ab10_v_Ab10_nucmer.coords", header=FALSE)
#Drop the first two lines of header
coords <- coords[-c(1,2,3),]
#Define the column names
names(coords) <- c("s1", "e1","s2","e2", "LEN1","LEN2", "%IDY","LENR", "LENQ", "COVR", "COVQ", "TAGS", "TAGS2")

#set columns as numeric
coords$s1 <- as.numeric(coords$s1)
coords$s2 <- as.numeric(coords$s2)
coords$e1 <- as.numeric(coords$e1)
coords$e2 <- as.numeric(coords$e2)
coords$`%IDY`<- as.numeric(coords$`%IDY`)
coords$LEN1<- as.numeric(coords$LEN1)
coords$LEN2<- as.numeric(coords$LEN1)
coords$COVR<- as.numeric(coords$COVR)

#filter to a percent identify of 80 and a length of 2500
coords_filt <- subset(coords, `%IDY` >= 80 & LEN1 >= 3000 & COVR >= 0.010)

#Mummer has known issues handeling tandem repeats, I am removing them from the plotting
coords_filt <- subset(coords_filt, !(s1 >= 142472000-141115174 & s1 <= 146699300-141115174) & !(s2 >= 142472000-141115174 & s2 <= 146699300-141115174) & !(s1 >= 150656000-141115174 & s1 <= 153145000-141115174) & !(s2 >= 150656000-141115174 & s2 <= 153145000-141115174) & !(s1 >= 157485200-141115174 & s1 <= 159356550-141115174) & !(s2 >= 157485200-141115174 & s2 <= 159356550-141115174) & !(s1 >= 174433450-141115174 & s1 <= 182846100-141115174) & !(s2 >= 174433450-141115174 & s2 <= 182846100-141115174))

#Convert to MB for easier
coords_filt$s1_MB <- coords_filt$s1/1000000
coords_filt$s2_MB <- coords_filt$s2/1000000

#####################################
#This sections generates the base of the dot plot
#####################################

Rgrob <- text_grob("R1", face = "bold", color = "orchid4", size = 25)
Tr1grob <- text_grob("TR1", face = "bold", color = "deepskyblue", rot=60, size = 25)
Sharedgrob <- text_grob("Shared", face = "bold", color = "darkgoldenrod1", rot=60, size = 25)
Trkingrob <- text_grob("trkin", face = "bold", color = "blue", rot=60, size = 25)
Kindrgrob <- text_grob("kindr", face = "bold", color = "hotpink", size = 25)
K180grob <- text_grob("Knob180", face = "bold", color = "darkorange3", size = 25)


#These are descriptions of relevant regions
pdf("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/Paper/TRKIN_Paper/Ab10_Ab10_DotPlot.pdf", height = 8, width = 8)
dotplot <- ggplot() +
  geom_point(data=coords_filt, alpha = 0.1, aes(x=s1_MB, y=s2_MB, color=`%IDY` )) +
  scale_colour_viridis_c(direction=-1, breaks = c(89, 99), labels = c("85","100")) +
  scale_x_continuous(breaks=c(0,10,20,30,40,50), limits=c(-1.8, 54)) +
  scale_y_continuous(breaks=c(0,10,20,30,40,50), limits=c(-12, 54)) +
  theme_bw() +
  labs(x='Ab10 (MB)', y='Ab10 (MB)', color="Percent Identity") +
  theme(axis.text = element_text(size = 30), axis.title = element_text(size = 35, face = "bold"), legend.text = element_text(size = 25), legend.position = "bottom", legend.title = element_text(size = 25, face = "bold")) +
  #####################Ab10 Annotations
  #R1
  annotation_custom(Rgrob, xmin=0, xmax=0, ymin=-2800000, ymax=-2800000) +
  ####TR1 Knob1
  geom_rect(mapping=aes(xmin=142.472000-141.115174, xmax=146.699300-141.115174, ymin=-2.000000, ymax=-0.050000), alpha=0.5, fill="deepskyblue", color = NA ) +
  geom_rect(mapping=aes(xmin=150.656000-141.115174, xmax=153.145000-141.115174, ymin=-2.000000, ymax=-0.050000), alpha=0.5, fill="deepskyblue", color = NA ) +
  geom_rect(mapping=aes(xmin=157.485200-141.115174, xmax=159.356550-141.115174, ymin=-2.000000, ymax=-0.050000), alpha=0.5, fill="deepskyblue", color = NA ) +
  annotation_custom(Tr1grob, xmin=142.472000-141.115174, xmax=146.699300-141.115174, ymin=-7, ymax=-5) +
  annotation_custom(Tr1grob, xmin=150.656000-141.115174, xmax=153.145000-141.115174, ymin=-7, ymax=-5) +
  annotation_custom(Tr1grob, xmin=157.485200-141.115174+1, xmax=159.356550-141.115174+1, ymin=-7, ymax=-5) +
  #####Knob 180
  geom_rect(mapping=aes(xmin=174.433450-141.115174, xmax=182.846100-141.115174, ymin=-2.000000, ymax=-0.050000), alpha=0.5, fill="darkorange3", color = NA ) +
  annotation_custom(K180grob, xmin=174.433450-141.115174, xmax=182.846100-141.115174, ymin=-4, ymax=-4) +
  ######Shared Region
  geom_rect(mapping=aes(xmin=0, xmax=142.472000-141.115174, ymin=-2.000000, ymax=-0.050000), alpha=0.5, fill="darkgoldenrod1", color = NA ) +
  geom_rect(mapping=aes(xmin=152.050000-141.115174, xmax=156.350000-141.115174, ymin=-2.000000, ymax=-0.050000), alpha=0.5, fill="darkgoldenrod1", color = NA ) +
  geom_rect(mapping=aes(xmin=158.250000-141.115174, xmax=166.820000-141.115174, ymin=-2.000000, ymax=-0.050000), alpha=0.5, fill="darkgoldenrod1", color = NA ) +
  annotation_custom(Sharedgrob, xmin=0-2, xmax=142.472000-141.115174-2, ymin=-10, ymax=-7) +
  annotation_custom(Sharedgrob, xmin=152.050000-141.115174, xmax=156.350000-141.115174, ymin=-10, ymax=-7) +
  annotation_custom(Sharedgrob, xmin=154.050000-141.115174+2, xmax=167.820000-141.115174+2, ymin=-10, ymax=-7) +
  #######trkin
  geom_rect(mapping=aes(xmin=148.964528-141.115174, xmax=149.082763-141.115174, ymin=-2.000000, ymax=-0.050000), fill="blue", color = NA ) +
  annotation_custom(Trkingrob, xmin=148.964528-141.115174-1, xmax=149.082763-141.115174-1, ymin=-7, ymax=-5) +
  geom_rect(mapping=aes(xmin=150.358622-141.115174, xmax=150.457752-141.115174, ymin=-2.000000, ymax=-0.050000), fill="blue", color = NA ) +
  #######kindr
  geom_rect(mapping=aes(xmin=189.326066-141.115174, xmax=190.330226-141.115174, ymin=-2.000000, ymax=-0.050000), alpha=0.5, fill="hotpink", color = NA ) +
  annotation_custom(Kindrgrob, xmin=188.326066-141.115174+1.8, xmax=189.330226-141.115174+1.8, ymin=-4, ymax=-4) +
  coord_cartesian(clip="off")
  dotplot
dev.off()

nrpd2e2grob <- text_grob("9 Copies\nnrpd2/e2", face = "bold", color = "red", size = 20)

pdf("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/Paper/TRKIN_Paper/Ab10_Ab10_DotPlot_nrpd2e2.pdf", height=8, width=8)
Share <- dotplot +
  scale_x_continuous(breaks=c(42,43,44,45,46,47,48,49), limits=c(41.730926, 48.210892)) +
  scale_y_continuous(breaks=c(42,43,44,45,46,47,48,49), limits=c(41.730926, 48.210892)) +
  annotation_custom(nrpd2e2grob, xmin=186.734627-141.115174, xmax=187.669213-141.115174, ymin=186.734627-141.115174-1.1, ymax=187.669213-141.115174-1.1)
  
Share
dev.off()

