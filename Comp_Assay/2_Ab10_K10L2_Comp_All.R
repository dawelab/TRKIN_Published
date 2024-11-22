library(readxl)
library(tidyverse)
library(ggplot2)
library(ggpubr)

setwd("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/Paper/TRKIN_Paper")

DATA_1 <- read_excel("~/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/Ab10_K10L2_Competition_Data_Reformat.xlsx")
DATA_1$Round <- 1

DATA_2 <- read_excel("~/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/Ab10_K10L2_Competition_Data_round2.xlsx")
DATA_2$Round <- 2
#This drops the Ab10+ K10L2- class because they all carried Cas9 which negated the results. 
DATA_2 <- subset(DATA_2, Genotype != "Ab10-I trkin + K10L2 trkin -")


DATA_3 <- read_excel("~/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/trkin_Competition_Assay_Summer_2024.xlsx")
DATA_3$Round <- 3

DATA <- rbind(DATA_1, DATA_2, DATA_3)

DATA <- subset(DATA, Total >= 50)

#Set the Genotype round and cas9 as a factor
DATA$Cas9 <- as.factor(DATA$Cas9)
DATA$Round <- as.factor(DATA$Round)
DATA$Genotype <- factor(DATA$Genotype, levels=c("Ab10-I trkin + K10L2 trkin +", "Ab10-I trkin - K10L2 trkin -", "Ab10-I trkin - K10L2 trkin +", "Ab10-I trkin + K10L2 trkin -", "Ab10-I trkin - N10", "Ab10-I trkin + N10"))

DATA$Drive <- as.numeric(DATA$Drive)
DATA$Total <- as.numeric(DATA$Total)

model <- aov(DATA$Drive ~ DATA$Cas9 + DATA$Round + DATA$Genotype, data = DATA)

summary(model)

SIG <- as.data.frame(TukeyHSD(model, conf.level=.95)[["DATA$Genotype"]])
SIG$sig <- ifelse(SIG$`p adj` <= 0.05 & SIG$`p adj` >= 0.01, "*", " ")
SIG$sig <- ifelse(SIG$`p adj` <= 0.01 & SIG$`p adj` >= 0.001, "**", SIG$sig)
SIG$sig <- ifelse(SIG$`p adj` <= 0.001 & SIG$`p adj` >= 0.0001, "***", SIG$sig)
SIG$sig <- ifelse(SIG$`p adj` <= 0.0001, "****", SIG$sig)

#This converts percentage to proportion
DATA$Prop_Drive <- DATA$R/DATA$Total

#Define the names to be more readable
x.labels <- c("Ab10\n1 + 2 -\n K10L2\n +", "Ab10\n1 - 2 -\n K10L2\n-", "Ab10\n1 - 2 -\nK10L2\n+", "Ab10\n1 + 2 -\n K10L2\n-", "Ab10\n1 - 2 -\nN10", "Ab10\n1 + 2 -\nN10")

pdf("Effect_Of_trkin_on_Ab10K10L2_Seg_AllComp.pdf", height = 6, width=9)
a <- ggplot(DATA, aes(x =Genotype, y = Prop_Drive)) +
  geom_jitter(aes(x =Genotype, y = Prop_Drive, color = Round, size= Total, shape = Cas9), alpha=0.7) +
  geom_boxplot(aes(x =Genotype, y = Prop_Drive), alpha = 0.5, outlier.shape = NA) +
  labs(x="trkin Genotype", y="Proportion Ab10-I", size = "Kernel\nNumber", color = "Season") +
  scale_x_discrete(labels = x.labels) +
  scale_y_continuous(breaks = c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
  scale_colour_viridis_d(option = "plasma", end = 0.7) +
  geom_signif(comparisons = list(c("Ab10-I trkin + K10L2 trkin +", "Ab10-I trkin - K10L2 trkin -")), annotation = c("***"), tip_length = 0.03, y_position = c(.90, .91)) +
  geom_signif(comparisons = list(c("Ab10-I trkin - K10L2 trkin -", "Ab10-I trkin - K10L2 trkin +")), annotation = c("****"), tip_length = 0.03, y_position = c(.90, .91)) +
  geom_signif(comparisons = list(c("Ab10-I trkin + K10L2 trkin -", "Ab10-I trkin - N10")), annotation = c("**"), tip_length = 0.03, y_position = c(.90, .91)) +
  geom_signif(comparisons = list(c("Ab10-I trkin - K10L2 trkin +", "Ab10-I trkin - N10")), annotation = c("****"), tip_length = 0.03, y_position = c(.95, .96)) +
  geom_signif(comparisons = list(c("Ab10-I trkin - N10", "Ab10-I trkin - K10L2 trkin -")), annotation = c("*"), tip_length = 0.03, y_position = c(1.00, 1.01)) +
  geom_signif(comparisons = list(c("Ab10-I trkin - K10L2 trkin +", "Ab10-I trkin + N10")), annotation = c("****"), tip_length = 0.03, y_position = c(1.05, 1.06)) +
  geom_signif(comparisons = list(c("Ab10-I trkin + K10L2 trkin +", "Ab10-I trkin - N10")), annotation = c("****"), tip_length = 0.03, y_position = c(1.10, 1.11)) +
  geom_signif(comparisons = list(c("Ab10-I trkin + K10L2 trkin +", "Ab10-I trkin + N10")), annotation = c("***"), tip_length = 0.03, y_position = c(1.15, 1.16)) +
  #geom_signif(comparisons = list(c("Ab10-I trkin + N10", "Ab10-I trkin - N10")), annotation = c("NS"), tip_length = 0.01, y_position = c(.90, .91)) +
  #geom_signif(comparisons = list(c("Ab10-I trkin - K10L2 trkin -", "Ab10-I trkin + K10L2 trkin -")), annotation = c("NS"), tip_length = 0.03, y_position = c(.90, .91)) +
  #geom_signif(comparisons = list(c("Ab10-I trkin + K10L2 trkin -", "Ab10-I trkin - K10L2 trkin +")), annotation = c("NS"), tip_length = 0.03, y_position = c(.90, .91)) +
  #geom_signif(comparisons = list(c("Ab10-I trkin + N10", "Ab10-I trkin - K10L2 trkin -")), annotation = c("NS"), tip_length = 0.03, y_position = c(1.00, 1.01)) +
  #geom_signif(comparisons = list(c("Ab10-I trkin + K10L2 trkin -", "Ab10-I trkin + N10")), annotation = c("NS"), tip_length = 0.03, y_position = c(1.15, 1.16)) +
  #geom_signif(comparisons = list(c("Ab10-I trkin + K10L2 trkin +", "Ab10-I trkin + K10L2 trkin -")), annotation = c("NS"), tip_length = 0.03, y_position = c(1.25, 1.26)) +
  #geom_signif(comparisons = list(c("Ab10-I trkin + K10L2 trkin +", "Ab10-I trkin - K10L2 trkin +")), annotation = c("NS"), tip_length = 0.03, y_position = c(1.30, 1.31)) +
  theme(axis.title = element_text(size = 20), axis.text.y = element_text(size = 18), axis.text.x = element_text(size = 18), legend.title = element_text(size = 20), legend.text = element_text(size = 18), legend.position = "right")
a
dev.off()



pdf("Effect_Of_trkin_on_Ab10K10L2_Seg_FocusedComp.pdf", height = 6, width=8)
b <- ggplot(DATA, aes(x =Genotype, y = Prop_Drive)) +
  geom_jitter(aes(x =Genotype, y = Prop_Drive, color = Round, size= Total), alpha=0.7) +
  geom_boxplot(aes(x =Genotype, y = Prop_Drive), alpha = 0.5, outlier.shape = NA) +
  labs(x="trkin Genotype", y="Proportion Ab10-I", size = "Kernel\nNumber", color = "Season") +
  scale_x_discrete(labels = x.labels) +
  scale_y_continuous(breaks = c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
  scale_colour_viridis_d(option = "plasma", end = 0.7) +
  geom_signif(comparisons = list(c("Ab10-I trkin + K10L2 trkin +", "Ab10-I trkin - K10L2 trkin -")), annotation = c("***"), tip_length = 0.03, y_position = c(.90, .91)) +
  #geom_signif(comparisons = list(c("Ab10-I trkin + K10L2 trkin +", "Ab10-I trkin - K10L2 trkin +")), annotation = c("NS"), tip_length = 0.03, y_position = c(.95, .96)) +
  #geom_signif(comparisons = list(c("Ab10-I trkin + K10L2 trkin +", "Ab10-I trkin + K10L2 trkin -")), annotation = c("NS"), tip_length = 0.03, y_position = c(.95, .96)) +
  geom_signif(comparisons = list(c("Ab10-I trkin + K10L2 trkin +", "Ab10-I trkin - N10")), annotation = c("****"), tip_length = 0.03, y_position = c(.95, .96)) +
  geom_signif(comparisons = list(c("Ab10-I trkin + K10L2 trkin +", "Ab10-I trkin + N10")), annotation = c("***"), tip_length = 0.03, y_position = c(1, 1.05)) +
  #geom_signif(comparisons = list(c("Ab10-I trkin + N10", "Ab10-I trkin - N10")), annotation = c("NS"), tip_length = 0.01, y_position = c(.90, .91)) +
  #geom_signif(comparisons = list(c("Ab10-I trkin - N10", "Ab10-I trkin - K10L2 trkin -")), annotation = c("*"), tip_length = 0.03, y_position = c(.90, .91)) +
  #geom_signif(comparisons = list(c("Ab10-I trkin - K10L2 trkin -", "Ab10-I trkin + K10L2 trkin -")), annotation = c("NS"), tip_length = 0.03, y_position = c(.90, .91)) +
  #geom_signif(comparisons = list(c("Ab10-I trkin + K10L2 trkin -", "Ab10-I trkin - K10L2 trkin +")), annotation = c("NS"), tip_length = 0.03, y_position = c(.90, .91)) +
#geom_signif(comparisons = list(c("Ab10-I trkin + K10L2 trkin -", "Ab10-I trkin - N10")), annotation = c("***"), tip_length = 0.03, y_position = c(.95, .96)) +
#geom_signif(comparisons = list(c("Ab10-I trkin + N10", "Ab10-I trkin - K10L2 trkin -")), annotation = c("NS"), tip_length = 0.03, y_position = c(1.00, 1.01)) +
#geom_signif(comparisons = list(c("Ab10-I trkin - K10L2 trkin -", "Ab10-I trkin - K10L2 trkin +")), annotation = c("****"), tip_length = 0.03, y_position = c(1.00, 1.01)) +
# geom_signif(comparisons = list(c("Ab10-I trkin - K10L2 trkin +", "Ab10-I trkin - N10")), annotation = c("****"), tip_length = 0.03, y_position = c(1.10, 1.11)) +
#geom_signif(comparisons = list(c("Ab10-I trkin + K10L2 trkin -", "Ab10-I trkin + N10")), annotation = c("NS"), tip_length = 0.03, y_position = c(1.15, 1.16)) +
# geom_signif(comparisons = list(c("Ab10-I trkin - K10L2 trkin +", "Ab10-I trkin + N10")), annotation = c("****"), tip_length = 0.03, y_position = c(1.20, 1.21)) +
  theme(axis.title = element_text(size = 20), axis.text.y = element_text(size = 18), axis.text.x = element_text(size = 18), legend.title = element_text(size = 20), legend.text = element_text(size = 18), legend.position = "bottom")

b
dev.off()

