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
coords <- read.delim("~/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/K10L2_DotPlot/K10L2_v_Ab10_nucmer.coords", header=FALSE)
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
#This sections generates the base of the dot plot
#####################################

Rgrob <- text_grob("R1", face = "bold", color = "orchid4")
Tr1grob <- text_grob("TR1", face = "bold", color = "deepskyblue")
Sharedgrob <- text_grob("Shared", face = "bold", color = "darkgoldenrod1")
Trkingrob <- text_grob("trkin", face = "bold", color = "blue")
Kindrgrob <- text_grob("kindr", face = "bold", color = "darkgreen")

Rgroby <- text_grob("R1", face = "bold", color = "orchid4", rot = 90)
Tr1groby <- text_grob("TR1", face = "bold", color = "deepskyblue", rot = 90)
Sharedgroby <- text_grob("Shared", face = "bold", color = "darkgoldenrod1", rot = 90)
Trkingroby <- text_grob("trkin", face = "bold", color = "blue", rot = 90)
Kindrgroby <- text_grob("kindr", face = "bold", color = "hotpink", rot = 90)
K180groby <- text_grob("Knob180", face = "bold", color = "darkorange3", rot = 90)

#These are descriptions of relevant regions
pdf("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/Paper/TRKIN_Paper/K10L2_Ab10_DotPlot.pdf", height = 12, width = 14)
dotplot <- ggplot() +
  geom_point(data=coords_filt, alpha = 0.3, aes(x=s1, y=s2, color=`%IDY` )) +
  scale_colour_viridis_c(direction=-1, breaks = c(87.67, 99), labels = c("85","100")) +
  theme_bw() +
  labs(x='K10L2', y='Ab10', color="Percent Identity") +
  xlim(-5000000, 40000000-2730186) +
  ylim(-3000000, 55000000) +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 24, face = "bold"), axis.text.x = element_text(angle=45, hjust=0.5, vjust=0.5), legend.text = element_text(size = 20), legend.position = "right", legend.title = element_text(size = 24, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5, size = 30, face = "bold")) +
  #R1
  annotation_custom(Rgrob, xmin=2736860-2730186, xmax=2736860-2730186, ymin=-1000000, ymax=-1000000) +
  #TR1 Knobs
  geom_rect(mapping=aes(xmin=7429619-2730186, xmax=13957822-2730186, ymin=-2000000, ymax=-50000), alpha=0.5, fill="deepskyblue", color = NA ) +
  annotation_custom(Tr1grob, xmin=7429619-2730186, xmax=13957822-2730186,  ymin=-2800000, ymax=-2800000) +
  geom_rect(mapping=aes(xmin=15726572-2730186, xmax=15898240-2730186, ymin=-2000000, ymax=-50000), alpha=0.5, fill="deepskyblue", color = NA ) +
  annotation_custom(Tr1grob, xmin=15726572-2730186, xmax=15898240-2730186,  ymin=-2800000*1.4, ymax=-2800000*1.4) +
  geom_rect(mapping=aes(xmin=16787371-2730186, xmax=25024178-2730186, ymin=-2000000, ymax=-50000), alpha=0.5, fill="deepskyblue", color = NA ) +
  annotation_custom(Tr1grob, xmin=16787371-2730186, xmax=25024178-2730186,  ymin=-2800000, ymax=-2800000) +
  geom_rect(mapping=aes(xmin=25498094-2730186, xmax=25502283-2730186, ymin=-2000000, ymax=-50000), alpha=0.5, fill="deepskyblue", color = NA ) +
  #annotation_custom(Tr1grob, xmin=25498094-2730186, xmax=25502283-2730186,  ymin=-2800000, ymax=-2800000) +
  #TRKIN
  geom_rect(mapping=aes(xmin=16328801-2730186, xmax=16357145-2730186, ymin=-2000000, ymax=-50000), fill="blue", color = NA ) +
  annotation_custom(Trkingrob, xmin=16328801-2730186, xmax=16357145-2730186,  ymin=-2800000, ymax=-2800000) +
  #Shared Region
  geom_rect(mapping=aes(xmin=2729416-2730186, xmax=7429619-2730186, ymin=-2000000, ymax=-50000), alpha=0.5, fill="darkgoldenrod1", color = NA ) +
  annotation_custom(Sharedgrob, xmin=2729416-2730186, xmax=7429619-2730186, ymin=-2800000, ymax=-2800000) +
  geom_rect(mapping=aes(xmin=25502283-2730186, xmax=31891546-2730186, ymin=-2000000, ymax=-50000), alpha=0.5, fill="darkgoldenrod1", color = NA ) +
  annotation_custom(Sharedgrob, xmin=25502283-2730186, xmax=31891546-2730186, ymin=-2800000, ymax=-2800000) + 
  
  #####################Ab10 Annotations
  #R1
  annotation_custom(Rgrob, ymin=0, ymax=0, xmin=-1000000, xmax=-1000000) +
  ####TR1 Knob1
  geom_rect(mapping=aes(ymin=142472000-141115174, ymax=146699300-141115174, xmin=-1000000, xmax=-50000), alpha=0.5, fill="deepskyblue", color = NA ) +
  geom_rect(mapping=aes(ymin=150656000-141115174, ymax=153145000-141115174, xmin=-1000000, xmax=-50000), alpha=0.5, fill="deepskyblue", color = NA ) +
  geom_rect(mapping=aes(ymin=157485200-141115174, ymax=159356550-141115174, xmin=-1000000, xmax=-50000), alpha=0.5, fill="deepskyblue", color = NA ) +
  geom_rect(mapping=aes(ymin=150656000-141115174, ymax=153145000-141115174, xmin=-1000000, xmax=-50000), alpha=0.5, fill="deepskyblue", color = NA ) +
  annotation_custom(Tr1groby, ymin=142472000-141115174, ymax=146699300-141115174, xmin=-2000000, xmax=-2000000) +
  annotation_custom(Tr1groby, ymin=150656000-141115174, ymax=153145000-141115174, xmin=-2000000, xmax=-2000000) +
  annotation_custom(Tr1groby, ymin=157485200-141115174, ymax=159356550-141115174, xmin=-2000000, xmax=-2000000) +
  #####Knob 180
  geom_rect(mapping=aes(ymin=174433450-141115174, ymax=182846100-141115174, xmin=-1000000, xmax=-50000), alpha=0.5, fill="darkorange3", color = NA ) +
  annotation_custom(K180groby, ymin=174433450-141115174, ymax=182846100-141115174, xmin=-2000000, xmax=-2000000) +
  ######Shared Region
  geom_rect(mapping=aes(ymin=0, ymax=142472000-141115174, xmin=-1000000, xmax=-50000), alpha=0.5, fill="darkgoldenrod1", color = NA ) +
  geom_rect(mapping=aes(ymin=10934826, ymax=15234826, xmin=-1000000, xmax=-50000), alpha=0.5, fill="darkgoldenrod1", color = NA ) +
  geom_rect(mapping=aes(ymin=17134826, ymax=25704826, xmin=-1000000, xmax=-50000), alpha=0.5, fill="darkgoldenrod1", color = NA ) +
  annotation_custom(Sharedgroby, ymin=10934826+500000, ymax=15234826+500000,xmin=-2000000, xmax=-2000000) +
  annotation_custom(Sharedgroby, ymin=17134826, ymax=25704826,xmin=-2000000, xmax=-2000000) +
  #######trkin
  geom_rect(mapping=aes(ymin=148964528-141115174, ymax=149082763-141115174, xmin=-1000000, xmax=-50000), fill="blue", color = NA ) +
  annotation_custom(Trkingroby, ymin=148964528-141115174, ymax=150457752-141115174, xmin=-2000000, xmax=-2000000) +
  geom_rect(mapping=aes(ymin=150358622-141115174, ymax=150457752-141115174, xmin=-1000000, xmax=-50000), fill="blue", color = NA ) +
  #annotation_custom(Trkingroby, ymin=150358622-141115174, ymax=150457752-141115174, xmin=-2000000, xmax=-2000000) +
  #######kindr
  geom_rect(mapping=aes(ymin=189326066-141115174, ymax=190330226-141115174, xmin=-1000000, xmax=-50000), alpha=0.5, fill="hotpink", color = NA ) +
  annotation_custom(Kindrgroby, ymin=188326066-141115174, ymax=190330226-141115174, xmin=-2000000, xmax=-2000000) +
  ######Kin10 
  #geom_rect(mapping=aes(ymin=190307761-141115174, ymax=191107857-141115174, xmin=-1000000, xmax=-50000), alpha=0.5, fill="hotpink", color = NA ) +
  #annotation_custom(Kin10groby, ymin=191305711-141115174, ymax=195107857-141115174,  xmin=-1500000, xmax=-1500000) +
  coord_cartesian(clip="off")
dotplot
dev.off()


#This zooms in on the trkin specific region
Trkin1grob <- text_grob("trkin 1", face = "bold", color = "black", size = 25)
Trkin2grob <- text_grob("trkin 2", face = "bold", color = "black", size = 25)

trkin

pdf("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/Paper/TRKIN_Paper/K10L2_Ab10_DotPlot_InterKnob.pdf", height=11, width=8)
Interknob <- dotplot +
  xlim(10000000, 15000000) +
  ylim(145000000-141115174, 155000000-141115174) +
  #######trkin
  annotate("segment", x =16357145-2730186+500000, xend =16328801-2730186+150000, y = 149023646-141115174, yend = 149023646-141115174, size = 1, arrow = arrow(angle = 45), color = "red") +
  annotate("segment", x =16328801-2730186-500000, xend =16357145-2730186-150000, y = 150408187-141115174, yend = 150408187-141115174, size = 1, arrow = arrow(angle = 45), color = "red") +
  geom_rect(mapping=aes(ymin=148964528-141115174, ymax=149082763-141115174, xmin=16328801-2730186, xmax=16357145-2730186 ), fill="red", color = NA, alpha = 0.5) +
geom_rect(mapping=aes(ymin=150358622-141115174, ymax=150457752-141115174, xmin=16328801-2730186, xmax=16357145-2730186), fill="red", color = NA, alpha = 0.5) +
  annotation_custom(Trkin1grob, ymin=148964528-141115174, ymax=149082763-141115174, xmin=16357145-2730186+1100000, xmax=16357145-2730186+1100000) +
  annotation_custom(Trkin2grob, ymin=150358622-141115174, ymax=150457752-141115174, xmin=16357145-2730186-1100000, xmax=16357145-2730186-1100000) +
  scale_colour_viridis_c(direction=-1, breaks = c(89, 99), labels = c("85","100")) +
  theme(axis.text = element_text(size = 30), axis.title = element_text(size = 34, face = "bold"), axis.text.x = element_text(angle=75, hjust=0.5, vjust=0.5), legend.text = element_text(size = 30), legend.position = "bottom", legend.title = element_text(size = 34, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5, size = 40, face = "bold"))
Interknob
dev.off()

