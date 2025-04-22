library(readxl)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(gvlma)
library(EnvStats)

setwd("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/Paper/TRKIN_Published/Comp_Assay")

DATA_1 <- read_excel("~/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/Ab10_K10L2_Competition_Data_Reformat.xlsx")
DATA_1$Round <- 1

DATA_2 <- read_excel("~/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/Ab10_K10L2_Competition_Data_round2.xlsx")
DATA_2$Round <- 2
#This drops the Ab10+ K10L2- class because they all carried Cas9 which negated the results. 
DATA_2 <- subset(DATA_2, Genotype != "Ab10-I trkin + K10L2 trkin -")


DATA_3 <- read_excel("~/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/trkin_Competition_Assay_Summer_2024.xlsx")
DATA_3$Round <- 3

DATA_4 <- read_excel("~/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/Ab10_K10L2_CompAssay_Round4.xlsx")
DATA_4$Round <- 4

DATA <- rbind(DATA_1, DATA_2, DATA_3, DATA_4)

DATA <- subset(DATA, Total >= 50)

#Set the Genotype round and cas9 as a factor
DATA$Cas9 <- as.factor(DATA$Cas9)
DATA$Round <- as.factor(DATA$Round)
DATA$Genotype <- factor(DATA$Genotype, levels=c("True Positive", "Ab10-I trkin + K10L2 trkin +", "Ab10-I trkin - K10L2 trkin -", "Ab10-I trkin - K10L2 trkin +", "Ab10-I trkin + K10L2 trkin -", "Ab10-I trkin - N10", "Ab10-I trkin + N10"))

DATA$Drive <- as.numeric(DATA$Drive)
DATA$Total <- as.numeric(DATA$Total)


model <- pairwise.wilcox.test(DATA$Drive, DATA$Genotype)

SIG <- as.data.frame(model$p.value)

#This writes out the significance table
write.csv(SIG, file="KW_Test_Sig.csv")
#I went in and manually edited this file for ease of interpretation

#This converts percentage to proportion
DATA$Prop_Drive <- DATA$R/DATA$Total

#Define the names to be more readable
x.labels <- c("Ab10\n1(+) 2(+)\n K10L2\n (+)", "Ab10\n1(+) 2(-)\n K10L2\n (+)", "Ab10\n1(-) 2(-)\n K10L2\n(-)", "Ab10\n1(-) 2(-)\nK10L2\n(+)", "Ab10\n1(+) 2(-)\n K10L2\n(-)", "Ab10\n1(-) 2(-)\nN10", "Ab10\n1(+) 2(-)\nN10")

SIG <- data.frame(Genotype = c("True Positive", "Ab10-I trkin + K10L2 trkin +", "Ab10-I trkin - K10L2 trkin -", "Ab10-I trkin - K10L2 trkin +", "Ab10-I trkin + K10L2 trkin -", "Ab10-I trkin - N10", "Ab10-I trkin + N10"), y=c(rep(1,7)), Label=c("a", "b", "c", "a,b", "d", "c", "c,d"))

pdf("Effect_Of_trkin_on_Ab10K10L2_Seg_AllComp.pdf", height = 6, width=9.5)
a <- ggplot() +
  geom_jitter(data=DATA, aes(x =Genotype, y = Prop_Drive, color = Round, size= Total, shape = Cas9), alpha=0.7) +
  geom_boxplot(data=DATA, aes(x =Genotype, y = Prop_Drive), alpha = 0.5, outlier.shape = NA) +
  labs(x="trkin Genotype", y="Proportion Ab10", size = "Kernel\nNumber", color = "Season") +
  scale_x_discrete(labels = x.labels) +
  scale_y_continuous(breaks = c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
  scale_colour_viridis_d(option = "plasma", end = 0.7) +
  geom_text(data=SIG, aes(x=Genotype, y=y, label=Label), size = 7) +
  theme(axis.title = element_text(size = 20), axis.text.y = element_text(size = 18), axis.text.x = element_text(size = 18), legend.title = element_text(size = 20), legend.text = element_text(size = 18), legend.position = "right")
a
dev.off()

