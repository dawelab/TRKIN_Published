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
coords <- read.delim("~/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/K10L2_DotPlot/N10_v_K10L2_nucmer.coords", header=FALSE)
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

coords_filt <- subset(coords, `%IDY` >= 80 & LEN1 >= 3000 & COVR >= 0.010)

#Mummer has known issues handeling tandem repeats, I am removing them from the plotting
coords_filt <- subset(coords_filt, !(s2 >= 7429619-2730186 & s2 <= 13957822-2730186) & !(s2 >= 15726572-2730186 & s2 <= 15898240-2730186) & !(s2 >= 16787371-2730186 & s2 <= 25024178-2730186) & !(s2 >= 25498094-2730186 & s2 <= 25502283-2730186))

#Convert to MB for easier plotting
coords_filt$s1_MB <- coords_filt$s1/1000000
coords_filt$s2_MB <- coords_filt$s2/1000000

#####################################
#This sections generates the base of the dot plot
#####################################

Rgrob <- text_grob("R1", face = "bold", color = "orchid4", rot = 90, size=25)
Tr1grob <- text_grob("TR1", face = "bold", color = "deepskyblue", rot = 90, size=25)
Sharedgrob <- text_grob("Shared", face = "bold", color = "darkgoldenrod1", rot = 90, size=25)
Trkingrob <- text_grob("trkin", face = "bold", color = "blue", rot = 90, size=25)



#These are descriptions of relevant regions
pdf("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/Paper/TRKIN_Paper/N10_K10L2_DotPlot.pdf", height=8, width=8)
dotplot <- ggplot() +
  geom_point(data=coords_filt, alpha = 0.2, aes(x=s1_MB, y=s2_MB, color=`%IDY` )) +
  scale_colour_viridis_c(direction=-1, breaks = c(91, 99.90), labels = c("85","100")) +
  scale_x_continuous(breaks=c(0,10,20,30), limits=c(-1.2, 12)) +
  scale_y_continuous(breaks=c(0,10,20,30), limits=c(-0.3, 30)) +
  theme_bw() +
  labs(x='N10 (MB)', y='K10L2 (MB)', color="Percent Identity") +
  theme(axis.text = element_text(size = 30), axis.title = element_text(size = 35, face = "bold"), legend.text = element_text(size = 25), legend.position = "bottom", legend.title = element_text(size = 25, face = "bold")) +
  #R1
  #annotation_custom(Rgrob, ymin=2.736860-2.730186, ymax=2.736860-2.730186, xmin=-1, xmax=-1) +
  #TR1 Knobs
  geom_rect(mapping=aes(ymin=7.429619-2.730186, ymax=13.957822-2.730186, xmin=-0.500000, xmax=-0.000500), alpha=0.5, fill="deepskyblue", color = NA ) +
  annotation_custom(Tr1grob, ymin=7.429619-2.730186, ymax=13.957822-2.730186, xmin=-1.1, xmax=-1.1) +
  
  geom_rect(mapping=aes(ymin=15.726572-2.730186, ymax=15.898240-2.730186, xmin=-0.500000, xmax=-0.000500), alpha=0.5, fill="deepskyblue", color = NA ) +
  #annotation_custom(Tr1grob, ymin=15726572-2730186, ymax=15898240-2730186, xmin=-800000, xmax=-800000) +
  
  geom_rect(mapping=aes(ymin=16.787371-2.730186, ymax=25.024178-2.730186, xmin=-0.500000, xmax=-0.000500), alpha=0.5, fill="deepskyblue", color = NA ) +
  annotation_custom(Tr1grob, ymin=16.787371-2.730186, ymax=25.024178-2.730186, xmin=-1.1, xmax=-1.1) +
  
  geom_rect(mapping=aes(ymin=25.498094-2.730186, ymax=25.502283-2.730186, xmin=-0.500000, xmax=-0.000500), alpha=0.5, fill="deepskyblue", color = NA ) +
  #annotation_custom(Tr1grob, ymin=25.498094-2.730186, ymax=25.502283-2.730186,  xmin=-1.1, xmax=-1.1) +
  
  #TRKIN
  geom_rect(mapping=aes(ymin=16.328801-2.730186, ymax=16.357145-2.730186, xmin=-0.500000, xmax=-0.000500), fill="blue", color = NA ) +
  annotation_custom(Trkingrob, ymin=16.328801-2.730186, ymax=16.357145-2.730186, xmin=-1.1, xmax=-1.1) +
  
  #Shared Region
  geom_rect(mapping=aes(ymin=2.729416-2.730186, ymax=7.429619-2.730186, xmin=-0.500000, xmax=-0.000500), alpha=0.5, fill="darkgoldenrod1", color = NA ) +
  annotation_custom(Sharedgrob, ymin=2.729416-2.730186, ymax=7.429619-2.730186, xmin=-1.1, xmax=-1.1) +
  
  geom_rect(mapping=aes(ymin=25.502283-2.730186, ymax=31.891546-2.730186, xmin=-0.500000, xmax=-0.000500), alpha=0.5, fill="darkgoldenrod1", color = NA ) +
  annotation_custom(Sharedgrob, ymin=25.502283-2.730186, ymax=31.891546-2.730186, xmin=-1.1, xmax=-1.1) +
  coord_cartesian(clip="off")
dotplot
dev.off()

jmjcgrob <- text_grob("3 Copies\nJmjC Containing\n Protein", face = "bold", color = "red", size = 20)


pdf("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/Paper/TRKIN_Paper/N10_K10L2_DotPlot_SharedDup.pdf", height=8, width=8)
Share <- dotplot +
  scale_x_continuous(breaks=c(7.75, 7.85, 7.95, 8), limits=c(7.75, 8)) +
  scale_y_continuous(breaks=c(24.5, 25, 25.5, 26), limits=c(24.5, 26)) +
  annotation_custom(jmjcgrob, xmin=7.85, xmax=7.85, ymin=24.5, ymax=25)
Share
dev.off()

