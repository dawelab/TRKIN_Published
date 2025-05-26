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
coords <- read.delim("~/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/K10L2_DotPlot/OGAb10_v_HiFiAb10_nucmer.coords", header=FALSE)
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
coords_filt <- subset(coords_filt, !(s1 >= 142324970-140935118 & s1 <= 146549396-140935118)
                      & !(s1 >= 150508711-140935118 & s1 <= 152995857-140935118)
                      & !(s1 >= 157286224-140935118 & s1 <= 159150445-140935118)
                      & !(s1 >= 174216963-140935118 & s1 <= 182945674-140935118)
                      
                      & !(s2 >= 142472000-141115174 & s2 <= 146699300-141115174)
                      & !(s2 >= 150656000-141115174 & s2 <= 153145000-141115174)
                      & !(s2 >= 157485200-141115174 & s2 <= 159356550-141115174)
                      & !(s2 >= 174433450-141115174 & s2 <= 182846100-141115174))

#Convert to MB for easier
coords_filt$s1_MB <- coords_filt$s1/1000000
coords_filt$s2_MB <- coords_filt$s2/1000000

#####################################
#This sections generates the base of the dot plot
#####################################

Tr1grob <- text_grob("TR1", face = "bold", color = "deepskyblue", rot=60, size = 25)
Sharedgrob <- text_grob("Shared", face = "bold", color = "darkgoldenrod1", rot=60, size = 25)
Trkingrob <- text_grob("trkin", face = "bold", color = "blue", rot=60, size = 25)
Kindrgrob <- text_grob("kindr", face = "bold", color = "hotpink", size = 25)
K180grob <- text_grob("Knob\n180", face = "bold", color = "darkorange3", size = 25)

Tr1groby <- text_grob("TR1", face = "bold", color = "deepskyblue", size = 25)
Sharedgroby <- text_grob("Shared", face = "bold", color = "darkgoldenrod1", size = 25)
Trkingroby <- text_grob("trkin", face = "bold", color = "blue", size = 25)
Kindrgroby <- text_grob("kindr", face = "bold", color = "hotpink", size = 25)
K180groby <- text_grob("Knob\n180", face = "bold", color = "darkorange3", size = 25)


#These are descriptions of relevant regions
pdf("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/Paper/TRKIN_Paper/OGAb10_HiFiAb10_DotPlot.pdf", height = 10, width = 10)
dotplot <- ggplot() +
  geom_point(data=coords_filt, alpha = 0.3, aes(x=s1_MB, y=s2_MB, color=`%IDY` )) +
  scale_colour_viridis_c(direction=-1, breaks = c(89, 99), labels = c("85","100")) +
  scale_x_continuous(breaks=c(0,10,20,30,40,50), limits=c(-13, 55)) +
  scale_y_continuous(breaks=c(0,10,20,30,40,50), limits=c(-13, 54)) +
  theme_bw() +
  labs(x='Ab10 v1 (MB)', y='Ab10 v2 (MB)', color="Percent Identity") +
  theme(axis.text = element_text(size = 30), axis.title = element_text(size = 35, face = "bold"), legend.text = element_text(size = 25), legend.position = "bottom", legend.title = element_text(size = 25, face = "bold")) +
  #####################OG Ab10 Annotations
  ####TR1 Knob1
  geom_rect(mapping=aes(xmin=142.324970-140.935118, xmax=146.549396-140.935118, ymin=-2.000000, ymax=-0.050000), alpha=0.5, fill="deepskyblue", color = NA ) +
  geom_rect(mapping=aes(xmin=150.508711-140.935118, xmax=152.995857-140.935118, ymin=-2.000000, ymax=-0.050000), alpha=0.5, fill="deepskyblue", color = NA ) +
  geom_rect(mapping=aes(xmin=157.286224-140.935118, xmax=159.150445-140.935118, ymin=-2.000000, ymax=-0.050000), alpha=0.5, fill="deepskyblue", color = NA ) +
  annotation_custom(Tr1grob, xmin=142.324970-140.935118, xmax=146.549396-140.935118, ymin=-7, ymax=-5) +
  annotation_custom(Tr1grob, xmin=150.508711-140.935118, xmax=152.995857-140.935118, ymin=-7, ymax=-5) +
  annotation_custom(Tr1grob, xmin=157.286224-140.935118+1, xmax=159.150445-140.935118+1, ymin=-7, ymax=-5) +
  #####Knob 180
  geom_rect(mapping=aes(xmin=174.216963-140.935118, xmax=182.945674-140.935118, ymin=-2.000000, ymax=-0.050000), alpha=0.5, fill="darkorange3", color = NA ) +
  annotation_custom(K180grob, xmin=174.216963-140.935118, xmax=182.945674-140.935118, ymin=-7, ymax=-7) +
  ######Shared Region
  geom_rect(mapping=aes(xmin=0, xmax=142.324970-140.935118, ymin=-2.000000, ymax=-0.050000), alpha=0.5, fill="darkgoldenrod1", color = NA ) +
  geom_rect(mapping=aes(xmin=152.050000-140.935118, xmax=156.880000-140.935118, ymin=-2.000000, ymax=-0.050000), alpha=0.5, fill="darkgoldenrod1", color = NA ) +
  geom_rect(mapping=aes(xmin=158.250000-140.935118, xmax=167.721000-140.935118, ymin=-2.000000, ymax=-0.050000), alpha=0.5, fill="darkgoldenrod1", color = NA ) +
  annotation_custom(Sharedgrob, xmin=0-2, xmax=142.324970-140.935118-2, ymin=-10, ymax=-7) +
  annotation_custom(Sharedgrob, xmin=152.050000-140.935118, xmax=156.880000-140.935118, ymin=-10, ymax=-7) +
  annotation_custom(Sharedgrob, xmin=158.250000-140.935118+1, xmax=167.721000-140.935118+1, ymin=-10, ymax=-7) +
  #######trkin
  geom_rect(mapping=aes(xmin=148.818350-140.935118, xmax=148.930352-140.935118, ymin=-2.000000, ymax=-0.050000), fill="blue", color = NA ) +
  annotation_custom(Trkingrob, xmin=148.818350-140.935118-1, xmax=150.308372-140.935118-1, ymin=-7, ymax=-5) +
  geom_rect(mapping=aes(xmin=150.208845-140.935118, xmax=150.308372-140.935118, ymin=-2.000000, ymax=-0.050000), fill="blue", color = NA ) +
  #######kindr
  geom_rect(mapping=aes(xmin=189.437000-140.935118, xmax=190.268773-140.935118, ymin=-2.000000, ymax=-0.050000), alpha=0.5, fill="hotpink", color = NA ) +
  annotation_custom(Kindrgrob, xmin=189.437000-140.935118+1.8, xmax=190.268773-140.935118+1.8, ymin=-5, ymax=-5) +
  
  
  
  
  
  #####################HiFi Ab10 Annotations

  ####TR1 Knob1
  geom_rect(mapping=aes(ymin=142.472000-141.115174, ymax=146.699300-141.115174, xmin=-2.000000, xmax=-0.050000), alpha=0.5, fill="deepskyblue", color = NA ) +
  geom_rect(mapping=aes(ymin=150.656000-141.115174, ymax=153.145000-141.115174, xmin=-2.000000, xmax=-0.050000), alpha=0.5, fill="deepskyblue", color = NA ) +
  geom_rect(mapping=aes(ymin=157.485200-141.115174, ymax=159.356550-141.115174, xmin=-2.000000, xmax=-0.050000), alpha=0.5, fill="deepskyblue", color = NA ) +
  annotation_custom(Tr1groby, ymin=142.472000-141.115174, ymax=146.699300-141.115174, xmin=-7, xmax=-5) +
  annotation_custom(Tr1groby, ymin=150.656000-141.115174, ymax=153.145000-141.115174, xmin=-7, xmax=-5) +
  annotation_custom(Tr1groby, ymin=157.485200-141.115174+1, ymax=159.356550-141.115174+1, xmin=-7, xmax=-5) +
  #####Knob 180
  geom_rect(mapping=aes(ymin=174.433450-141.115174, ymax=182.846100-141.115174, xmin=-2.000000, xmax=-0.050000), alpha=0.5, fill="darkorange3", color = NA ) +
  annotation_custom(K180groby, ymin=174.433450-141.115174, ymax=182.846100-141.115174, xmin=-7, xmax=-7) +
  ######Shared Region
  geom_rect(mapping=aes(ymin=0, ymax=142.472000-141.115174, xmin=-2.000000, xmax=-0.050000), alpha=0.5, fill="darkgoldenrod1", color = NA ) +
  geom_rect(mapping=aes(ymin=152.050000-141.115174, ymax=156.350000-141.115174, xmin=-2.000000, xmax=-0.050000), alpha=0.5, fill="darkgoldenrod1", color = NA ) +
  geom_rect(mapping=aes(ymin=158.250000-141.115174, ymax=166.820000-141.115174, xmin=-2.000000, xmax=-0.050000), alpha=0.5, fill="darkgoldenrod1", color = NA ) +
  annotation_custom(Sharedgroby, ymin=0-2, ymax=142.472000-141.115174-2, xmin=-10, xmax=-7) +
  annotation_custom(Sharedgroby, ymin=152.050000-141.115174+2, ymax=156.350000-141.115174+2, xmin=-10, xmax=-7) +
  annotation_custom(Sharedgroby, ymin=154.050000-141.115174+2, ymax=167.820000-141.115174+2, xmin=-10, xmax=-7) +
  #######trkin
  geom_rect(mapping=aes(ymin=148.964528-141.115174, ymax=149.082763-141.115174, xmin=-2.000000, xmax=-0.050000), fill="blue", color = NA ) +
  annotation_custom(Trkingroby, ymin=148.964528-141.115174-1, ymax=149.082763-141.115174-1, xmin=-7, xmax=-5) +
  geom_rect(mapping=aes(ymin=150.358622-141.115174, ymax=150.457752-141.115174, xmin=-2.000000, xmax=-0.050000), fill="blue", color = NA ) +
  #######kindr
  geom_rect(mapping=aes(ymin=189.326066-141.115174, ymax=190.330226-141.115174, xmin=-2.000000, xmax=-0.050000), alpha=0.5, fill="hotpink", color = NA ) +
  annotation_custom(Kindrgroby, ymin=188.326066-141.115174+1.8, ymax=189.330226-141.115174+1.8, xmin=-7, xmax=-7) +
  coord_cartesian(clip="off")
  dotplot
dev.off()


Kindrgrob <- text_grob("kindr", face = "bold", color = "red", size = 25)

pdf("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/Paper/TRKIN_Paper/OGAb10_HiFiAb10_DotPlot_DistalTip.pdf", height=10, width=10)
Distal <- dotplot +
  scale_x_continuous(breaks=c(42:54), limits=c(42,54.5)) +
  scale_y_continuous(breaks=c(42:54), limits=c(42, 54)) +
  annotation_custom(Kindrgrob, xmin=48, xmax=50, ymin=47, ymax=48)
Distal
dev.off()

