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
coords <- read.delim("K10L2_v_Ab10_nucmer.coords", header=FALSE)
#This drops the first row
coords <- coords[-c(1),]
names(coords) <- c("s1", "e1","s2","e2", "LEN1","LEN2", "%IDY","LENR", "LENQ", "COVR", "COVQ", "TAGS", "TAGS2")
coords$s1 <- as.numeric(coords$s1)
coords$s2 <- as.numeric(coords$s2)
coords$e1 <- as.numeric(coords$e1)
coords$e2 <- as.numeric(coords$e2)
coords$`%IDY`<- as.numeric(coords$`%IDY`)
coords$LEN1<- as.numeric(coords$LEN1)
coords$LEN2<- as.numeric(coords$LEN1)

coords_filt <- subset(coords, `%IDY` >= 80 & LEN1 >= 2500)

#####################################
#This section calculates the average percent identity across the trkin bearing region
#####################################
coords_filt2 <-  subset(coords_filt, s1 >= 11227636 & e1 <= 14057185 & s2 >= 5584126 & e2 <= 9540826)

#This plots the distribution of percent identity 
hist(coords_filt2$`%IDY`)
#This calculates the mean percent identity
mean(coords_filt2$`%IDY`)

#Mummer has known issues handeling tandem repeats, I am removing them from the plotting
coords_filt <- subset(coords_filt, !(s2 >= 142472000-141115174 & s2 <= 146699300-141115174)  & !(s2 >= 150656000-141115174 & s2 <= 153145000-141115174) & !(s2 >= 157485200-141115174 & s2 <= 159356550-141115174) & !(s2 >= 174433450-141115174 & s2 <= 182846100-141115174) & !(s1 >= 7429619-2730186 & s1 <= 13957822-2730186) & !(s1 >= 15726572-2730186 & s1 <= 15898240-2730186) & !(s1 >= 16787371-2730186 & s1 <= 25024178-2730186) & !(s1 >= 25498094-2730186 & s1 <= 25502283-2730186))

#####################################
#This sections generates the base of the dot plot
#####################################

Rgrob <- text_grob("R1", face = "bold", color = "orchid4")
Tr1grob <- text_grob("TR-1", face = "bold", color = "deepskyblue")
Sharedgrob <- text_grob("Shared", face = "bold", color = "darkgoldenrod1")
Trkingrob <- text_grob("trkin", face = "bold", color = "blue")
Kindrgrob <- text_grob("kindr", face = "bold", color = "darkgreen")

Rgroby <- text_grob("R1", face = "bold", color = "orchid4", rot = 90)
Tr1groby <- text_grob("TR-1", face = "bold", color = "deepskyblue", rot = 90)
Sharedgroby <- text_grob("Shared", face = "bold", color = "darkgoldenrod1", rot = 90)
Trkingroby <- text_grob("trkin", face = "bold", color = "blue", rot = 90)
Kindrgroby <- text_grob("kindr", face = "bold", color = "hotpink", rot = 90)
K180groby <- text_grob("Knob180", face = "bold", color = "darkorange3", rot = 90)

#Convert to MB for easier
coords_filt$s1_MB <- coords_filt$s1/1000000
coords_filt$s2_MB <- coords_filt$s2/1000000

#These are descriptions of relevant regions
pdf("K10L2_Ab10_DotPlot.pdf", height = 12, width = 14)
dotplot <- ggplot() +
  geom_point(data=coords_filt, alpha = 0.3, aes(x=s1_MB, y=s2_MB, color=`%IDY` )) +
  scale_colour_viridis_c(direction=-1, breaks = c(87.67, 99), labels = c("85","100")) +
  theme_bw() +
  labs(x='K10L2', y='Ab10', color="Percent Identity") +
  xlim(-5.000000, 40.000000-2.730186) +
  ylim(-3.000000, 55.000000) +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 24, face = "bold"), axis.text.x = element_text(angle=45, hjust=0.5, vjust=0.5), legend.text = element_text(size = 20), legend.position = "right", legend.title = element_text(size = 24, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5, size = 30, face = "bold")) +
  #R1
  annotation_custom(Rgrob, xmin=2.736860-2.730186, xmax=2.736860-2.730186, ymin=-1.000000, ymax=-1.000000) +
  #TR1 Knobs
  geom_rect(mapping=aes(xmin=7.429619-2.730186, xmax=13.957822-2.730186, ymin=-2.000000, ymax=-0.050000), alpha=0.5, fill="deepskyblue", color = NA ) +
  annotation_custom(Tr1grob, xmin=15.726572-2.730186, xmax=15.898240-2.730186,  ymin=-2.800000*1.4, ymax=-2.800000*1.4) +
  geom_rect(mapping=aes(xmin=16.787371-2.730186, xmax=25.024178-2.730186, ymin=-2.000000, ymax=-0.050000), alpha=0.5, fill="deepskyblue", color = NA ) +
  annotation_custom(Tr1grob, xmin=16.787371-2.730186, xmax=25.024178-2.730186,  ymin=-2.800000, ymax=-2.800000) +
  geom_rect(mapping=aes(xmin=25.498094-2.730186, xmax=25.502283-2.730186, ymin=-2.000000, ymax=-0.050000), alpha=0.5, fill="deepskyblue", color = NA ) +
  #annotation_custom(Tr1grob, xmin=25.498094-2.730186, xmax=25.502283-2.730186,  ymin=-2.800000, ymax=-2.800000) +
  #TRKIN
  geom_rect(mapping=aes(xmin=16.328801-2.730186, xmax=16.357145-2.730186, ymin=-2.000000, ymax=-0.050000), fill="blue", color = NA ) +
  annotation_custom(Trkingrob, xmin=16.328801-2.730186, xmax=16.357145-2730186,  ymin=-2.800000, ymax=-2.800000) +
  #Shared Region
  geom_rect(mapping=aes(xmin=2.729416-2.730186, xmax=7.429619-2.730186, ymin=-2.000000, ymax=-0.050000), alpha=0.5, fill="darkgoldenrod1", color = NA ) +
  annotation_custom(Sharedgrob, xmin=2.729416-2.730186, xmax=7.429619-2.730186, ymin=-2.800000, ymax=-2.800000) +
  geom_rect(mapping=aes(xmin=25.502283-2.730186, xmax=31.891546-2.730186, ymin=-2.000000, ymax=-0.050000), alpha=0.5, fill="darkgoldenrod1", color = NA ) +
  annotation_custom(Sharedgrob, xmin=25.502283-2.730186, xmax=31.891546-2.730186, ymin=-2.800000, ymax=-2.800000) +

  #####################Ab10 Annotations
  #R1
  annotation_custom(Rgrob, ymin=0, ymax=0, xmin=-1000000, xmax=-1000000) +
  ####TR1 Knob1
  geom_rect(mapping=aes(ymin=142.472000-141.115174, ymax=146.699300-141.115174, xmin=-1.000000, xmax=-0.050000), alpha=0.5, fill="deepskyblue", color = NA ) +
  geom_rect(mapping=aes(ymin=150.656000-141.115174, ymax=153.145000-141.115174, xmin=-1.000000, xmax=-0.050000), alpha=0.5, fill="deepskyblue", color = NA ) +
  geom_rect(mapping=aes(ymin=157.485200-141.115174, ymax=159.356550-141.115174, xmin=-1.000000, xmax=-0.050000), alpha=0.5, fill="deepskyblue", color = NA ) +
  geom_rect(mapping=aes(ymin=150.656000-141.115174, ymax=153.145000-141.115174, xmin=-1.000000, xmax=-0.050000), alpha=0.5, fill="deepskyblue", color = NA ) +
  annotation_custom(Tr1groby, ymin=142.472000-141.115174, ymax=146.699300-141.115174, xmin=-2.000000, xmax=-2.000000) +
  annotation_custom(Tr1groby, ymin=150.656000-141.115174, ymax=153.145000-141.115174, xmin=-2.000000, xmax=-2.000000) +
  annotation_custom(Tr1groby, ymin=157.485200-141.115174, ymax=159.356550-141.115174, xmin=-2.000000, xmax=-2.000000) +
  #####Knob 180
  geom_rect(mapping=aes(ymin=174.433450-141.115174, ymax=182.846100-141.115174, xmin=-1.000000, xmax=-0.050000), alpha=0.5, fill="darkorange3", color = NA ) +
  annotation_custom(K180groby, ymin=174.433450-141.115174, ymax=182.846100-141.115174, xmin=-2.000000, xmax=-2.000000) +
  ######Shared Region
  geom_rect(mapping=aes(ymin=0, ymax=142.472000-141.115174, xmin=-1.000000, xmax=-0.050000), alpha=0.5, fill="darkgoldenrod1", color = NA ) +
  geom_rect(mapping=aes(ymin=10.934826, ymax=15.234826, xmin=-1.000000, xmax=-0.050000), alpha=0.5, fill="darkgoldenrod1", color = NA ) +
  geom_rect(mapping=aes(ymin=17.134826, ymax=25.704826, xmin=-1.000000, xmax=-0.050000), alpha=0.5, fill="darkgoldenrod1", color = NA ) +
  annotation_custom(Sharedgroby, ymin=10.934826+0.500000, ymax=15.234826+0.500000,xmin=-2.000000, xmax=-2.000000) +
  annotation_custom(Sharedgroby, ymin=17.134826, ymax=25.704826,xmin=-2.000000, xmax=-2.000000) +
  #######trkin
  geom_rect(mapping=aes(ymin=148.964528-141.115174, ymax=149.082763-141.115174, xmin=-1.000000, xmax=0.050000), fill="blue", color = NA ) +
  annotation_custom(Trkingroby, ymin=148.964528-141.115174, ymax=150.457752-141.115174, xmin=-2.000000, xmax=-2.000000) +
  geom_rect(mapping=aes(ymin=150.358622-141.115174, ymax=150.457752-141.115174, xmin=-1.000000, xmax=-0.050000), fill="blue", color = NA ) +
  #annotation_custom(Trkingroby, ymin=150.358622-141.115174, ymax=150.457752-141.115174, xmin=-2.000000, xmax=-2.000000) +
  #######kindr
  geom_rect(mapping=aes(ymin=189.326066-141.115174, ymax=190.330226-141.115174, xmin=-1.000000, xmax=-0.050000), alpha=0.5, fill="hotpink", color = NA ) +
  annotation_custom(Kindrgroby, ymin=188.326066-141.115174, ymax=190.330226-141.115174, xmin=-2.000000, xmax=-2.000000) +
  ######Kin10
  #geom_rect(mapping=aes(ymin=190.307761-141.115174, ymax=191.107857-141.115174, xmin=-1.000000, xmax=-0.050000), alpha=0.5, fill="hotpink", color = NA ) +
  #annotation_custom(Kin10groby, ymin=191.305711-141.115174, ymax=195.107857-141.115174,  xmin=-1.500000, xmax=-1.500000) +
  coord_cartesian(clip="off")
dotplot
dev.off()


#This zooms in on the trkin specific region
Trkin1grob <- text_grob("Trkin1", face = "italic", color = "black", size = 8)
Trkin2grob <- text_grob("Trkin2", face = "italic", color = "black", size = 8)

pdf("K10L2_Ab10_DotPlot_InterKnob.pdf", height=3, width=3)
Interknob <- dotplot +
  xlim(11, 14.5) +
  ylim(5, 10) +
  #######trkin
  annotate("segment", x =16.357145-2.730186+0.500000, xend =16.328801-2.730186+0.150000, y = 149.023646-141.115174, yend = 149.023646-141.115174, size = 0.5, arrow = arrow(angle = 45, length = unit(0.25, "cm")), color = "red") +
  annotate("segment", x =16.328801-2.730186-0.500000, xend =16.357145-2.730186-0.150000, y = 150.408187-141.115174, yend = 150.408187-141.115174, size = 0.5, arrow = arrow(angle = 45, length = unit(0.25, "cm")), color = "red") +
  annotation_custom(Trkin1grob, ymin=148.964528-141.115174, ymax=149.082763-141.115174, xmin=14.2, xmax=14.6) +
  annotation_custom(Trkin2grob, ymin=150.358622-141.115174, ymax=150.457752-141.115174, xmin=12.5, xmax=13.0) +
  scale_colour_viridis_c(direction=-1, breaks = c(90, 95, 99.8), labels=c(90, 95, 100)) +
  theme(axis.text = element_text(size = 8), axis.text.x = element_text(angle=0, hjust=0.5, vjust=0.5), axis.title = element_text(size = 10, face = "plain"), legend.text = element_text(size = 8), legend.position = "bottom", legend.title = element_text(size = 10, face = "plain")) +
  theme(plot.title = element_text(hjust = 0.5, size = 40))
Interknob
dev.off()

