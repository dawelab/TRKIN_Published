library(ggplot2)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(readxl)
library(gvlma)
library(gamlss)

#This loads in the data from the seed counts
DATA <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/Hawaii_Seed_Count.csv")
colnames(DATA) <- c("Number", "Individual", "R_(purple)", "R_Counter", "r_(yellow)", "r_Counter", "Total", "Sorter", "Sort_Date", "Date_Counted", "Person_Counting", "Drive", "Number_Defective_Kernels", "Prop_Defective", "Defective_Date", "Defective_Number_that_Day", "Notes")

#This loads in the field data
MAP <- read_excel("~/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/Hawaii_MappingFile_v2.xlsx")
colnames(MAP) <- c("Genotype", "Number", "MB_Number", "Event_Number", "Column_Number", "Row_Number", "Edge")

#This merges the seed count and field data
BOTH <- merge(DATA, MAP, by="Number")

#This converts several columns to the appropriate class
BOTH$Genotype <- as.factor(BOTH$Genotype)
BOTH$R_Counter <- as.factor(BOTH$R_Counter)
BOTH$r_Counter <- as.factor(BOTH$R_Counter)
BOTH$Sorter <- as.factor(BOTH$Sorter)
BOTH$Sort_Date <- as.factor(BOTH$Sort_Date)
BOTH$Date_Counted <- as.factor(BOTH$Date_Counted)
BOTH$Person_Counting <- as.factor(BOTH$Person_Counting)
BOTH$Defective_Date <- as.factor(BOTH$Defective_Date)


BOTH$`R_(purple)` <- as.numeric(BOTH$`R_(purple)`)
BOTH$`r_(yellow)` <- as.numeric(BOTH$`r_(yellow)`)
BOTH$Drive <- as.numeric(BOTH$Drive)
BOTH$Number_Defective_Kernels <- as.numeric(BOTH$Number_Defective_Kernels)
BOTH$Prop_Defective <- as.numeric(BOTH$Prop_Defective)
BOTH$Total <- as.numeric(BOTH$Total)
BOTH$Defective_Number_that_Day <- as.numeric(BOTH$Defective_Number_that_Day)

#This generates a histogram of total ear size
ggplot(data=BOTH, aes(x=Total)) +
  geom_histogram(bins=50)
ggsave("Hawaii_TotalHistogram.png")

ggplot(data=BOTH, aes(x=Total, y=Drive)) +
  geom_point()
ggsave("Hawaii_TotalvDrive.png")

#I am going to filter any ears below 50 total kernels, I don't feel comfortable scoring drive on anything less than that. There doesn't seem to be any consistent trend with very small ears having nosier drive. I am also dropping loose kernels because I don't think drive seen there is biologically relevant

BOTH_sub <- subset(BOTH, Total >= 50 & Individual != "loose")

#######################################################################################
#This tests if which corn counter was used has a significant effect on the number of kernels counted
#######################################################################################
COUNT <- data.frame(Count = c(BOTH_sub$`R_(purple)`, BOTH_sub$`r_(yellow)`), Counter = c(BOTH_sub$R_Counter, BOTH_sub$r_Counter))

COUNT$Count <- as.integer(COUNT$Count)
COUNT$Counter <- as.factor(COUNT$Counter)

#This isn't data isn't quite normally distributed and neither are the residuals, but I think that it is close enough that it is better to use the parametric statistic
hist(sqrt(COUNT$Count))
test <- lm(sqrt(Count) ~ Counter, data = COUNT)
gvlma(test)


testmodel <- aov(sqrt(Count) ~ Counter, data = COUNT)
summary(testmodel)

hist(residuals(testmodel))

#####################################
#Counter has no effect
#####################################


#######################################################################################
#This tests if the person sorting the corn affected the drive obtained
#######################################################################################

SORT <- BOTH_sub[,c("Drive", "Sorter", "Sort_Date")]

#drive is not quite normally distributed, but this transformation gets me as close as I can
sortmodel <- lm(sqrt(max(SORT$Drive+1) - SORT$Drive) ~ Sorter, Sort_Date, data = SORT)
gvlma(sortmodel)
hist(sqrt(max(SORT$Drive+1) - SORT$Drive))
hist(residuals(sortmodel))

sortmodel <- aov(sqrt(max(SORT$Drive+1) - SORT$Drive) ~ Sorter + Sort_Date, data = SORT)
summary(sortmodel)

#####################################
#Neither Sorter nor sort date have a significant effect 
#####################################



#######################################################################################
#Drive does not meet ANOVA assumtions but this transformation gets me as close as I can
#######################################################################################
#################This one gets me closest, but doesn't actually fit all assumptions 
#this plots the data 

#This drops samples with less than 50% drive which would be recombinats that do not carry Ab10. This drops only one sample
BOTH_sub_filt <- subset(BOTH_sub, Drive >= 0.50)

#This generates histograms for each samples drive 
ggplot(data = BOTH_sub_filt, aes(x=Drive)) +
  geom_histogram(bins=50)
ggsave("UnTranDrive_Hist.png")

ggplot(data = BOTH_sub_filt, aes(x=Drive)) +
  facet_wrap(BOTH_sub_filt$Genotype) +
  geom_histogram(bins=50)
ggsave("UnTranDrive_HistByGenotype.png")

#This runs the largest model to test assumptions
fit1=lm(Drive ~ Column_Number + Row_Number + Edge + Sorter + Sort_Date + Genotype, data = BOTH_sub_filt)
gvmodelfit1 <- gvlma(fit1)
gvmodelfit1
#This model does not meet the assumptions

#This generates a q:q plot for the residuals against the normal distribution. It's clearly not normal but it't pretty close, it looks like the problem is outliers 
qqnorm(residuals(fit1))

model_1 <- aov(Drive~ Column_Number + Row_Number + Edge + Sorter + Genotype, data = BOTH_sub_filt)
summary(model_1)
TukeyHSD(model_1, which="Genotype")
#With unfiltered, untransformed data TRKIN mutation 2 has significantly higher drive than both functional trkin 

#This square root transorms the data in an effort to help with the slight positive skew. It also drops the samples with <55% drive, which are likely N10 recombinants 
BOTH_sub_filt$Adj_Drive <- -(sqrt(max(BOTH_sub_filt$Drive+1) - BOTH_sub_filt$Drive))
BOTH_sub_filt <- subset(BOTH_sub_filt, Drive >= 0.55)

#This plots the data transformation
ggplot(data=BOTH_sub_filt, aes(x=Drive, y=Adj_Drive)) +
  facet_wrap(BOTH_sub_filt$Genotype) +
  geom_point()
ggsave("TranDrive_Point.png")

#This generates histograms for each sample after transformation 
ggplot(data = BOTH_sub_filt, aes(x=Adj_Drive)) +
  geom_histogram(bins=50)
ggsave("TranDrive_Hist.png")

ggplot(data = BOTH_sub_filt, aes(x=Adj_Drive)) +
  facet_wrap(BOTH_sub_filt$Genotype) +
  geom_histogram(bins=50)
ggsave("TranDrive_HistByGenotype.png")

#This runs the largest model to test assumptions
fit1=lm(Adj_Drive ~ Column_Number + Row_Number + Edge + Sorter + Sort_Date + Genotype, data = BOTH_sub_filt)
gvmodelfit1 <- gvlma(fit1)
gvmodelfit1
#This model does not meet the assumptions, but it is slightly closer to meeting the assumptions than the untransformed unfiltered data 

#This generates a q:q plot for the residuals against the normal distribution. It's clearly not normal but it't pretty close, it looks like the problem is outliers 
qqnorm(residuals(fit1))

#I am chooseing to go with this model which drops the sample with abberantly low drive and transmorms the data to reduce the skew
#This model includes all possible variable affecting drive
model_1 <- aov(Adj_Drive~ Column_Number + Row_Number + Edge + Sorter + Genotype, data = BOTH_sub_filt)
summary(model_1)
TukeyHSD(model_1, which="Genotype")

model_1 <- aov(Drive~ Column_Number + Row_Number + Edge + Sorter + Genotype, data = BOTH_sub_filt)
summary(model_1)
TukeyHSD(model_1, which="Genotype")
#With filtered and transformed data TRKIN mutation 2 has significantly higher drive than both functional trkin, but is not different from trkin mutation 1

Adj_Drive_Sig <- data.frame(x=c(1, 2, 3), y=c(1, 1, 1), Sig=c("a", "a", "b") )

pdf("Ab10_Drive_By_trkin_Genotype.pdf", height=10, width=4)
ggplot() +
  geom_jitter(data=BOTH_sub_filt, aes(x=Genotype, y=Drive), height = 0, color= "grey60", alpha=0.5) +
  geom_boxplot(data=BOTH_sub_filt, aes(x=Genotype, y=Drive), alpha=0.5, outlier.shape = NA) +
  geom_text(data=Adj_Drive_Sig, aes(x=x, y=y, label=Sig), size = 10) +
  scale_x_discrete(labels=c("TRKIN Mutation 2" = "1 -\n2 -", "TRKIN Mutation 1" = "1 -\n2 +", "Functional TRKIN" = "1 +\n2 -")) +
  xlab("Ab10-I trkin Genotype") +
  ylab("Proportion Ab10-I") +
  geom_vline(xintercept=1.5, color = "grey80", linetype="dashed", size = 0.5) +
  geom_vline(xintercept=2.5, color = "grey80", linetype="dashed", size = 0.5) +
  theme(plot.title = element_text(hjust=0.5, size = 30), plot.subtitle = element_text(hjust=0.5, size = 10), axis.title = element_text(size = 20), axis.text.y = element_text(size = 20), axis.text.x = element_text(size = 20), legend.title = element_text(size = 20), legend.text = element_text(size = 20), legend.position = "right") 
dev.off()

#####################################
#Summary: Drive is not normally distributed and does not meet all the assumptions of an anova, however the data appears nearly normal. Thus I chose to use parametric statistics for this test. In all models drive in trkin mutation 1 is not different from frunction trkin. Drive in trkin mutation 2 is consistently significantly higher from functional trkin, in some models it is significantly higher that trkin mutation 1. 
#####################################


#######################################################################################
#Total does meet ANOVA assumptions
#######################################################################################
#I am dropping samples with less than 55% drive in the event that they are recombinants not carrying Ab10 due to the possible fitness affects
ggplot(data = BOTH_sub_filt, aes(x=Total)) +
  geom_histogram(bins=30)
ggsave("Total_Hist.png")


ggplot(data = BOTH_sub_filt, aes(x=Total)) +
  facet_wrap(BOTH_sub_filt$Genotype) +
  geom_histogram(bins=30)
ggsave("Total_HistByGenotype.png")

BOTH_sub_filt <- subset(BOTH_sub, Drive >= 0.55)

fit2=lm(Total ~ Genotype, data = BOTH_sub_filt)
gvmodelfit2 <- gvlma(fit2)
gvmodelfit2

model_2 <- aov(Total ~ Column_Number + Row_Number + Edge + Sorter + Genotype, data = BOTH_sub_filt)
summary(model_2)
TukeyHSD(model_1, which="Sorter")

#####################################
#Only edge has a significant effect here
#####################################

pdf("KernelNumber_By_trkin_Genotype.pdf", height=10, width=4)
ggplot() +
  geom_jitter(data=BOTH_sub_filt, aes(x=Genotype, y=Total), height = 0, color= "grey60", alpha=0.5) +
  geom_boxplot(data=BOTH_sub_filt, aes(x=Genotype, y=Total), alpha=0.5, outlier.shape = NA) +
  scale_x_discrete(labels=c("TRKIN Mutation 2" = "1 -\n2 -", "TRKIN Mutation 1" = "1 -\n2 +", "Functional TRKIN" = "1 +\n2 -")) +
  xlab("Ab10-I trkin Genotype") +
  ylab("Kernel Number") +
  geom_vline(xintercept=1.5, color = "grey80", linetype="dashed", size = 0.5) +
  geom_vline(xintercept=2.5, color = "grey80", linetype="dashed", size = 0.5) +
  theme(plot.title = element_text(hjust=0.5, size = 30), plot.subtitle = element_text(hjust=0.5, size = 10), axis.title = element_text(size = 20), axis.text.y = element_text(size = 20), axis.text.x = element_text(size = 20), legend.title = element_text(size = 20), legend.text = element_text(size = 20), legend.position = "right") 
dev.off()

#####################################
#Summary: Trkin genotype has no affect on total kernel number 
#####################################

#This generates histograms for each samples proportion of defective kernels
ggplot(data = BOTH_sub_filt, aes(x=Prop_Defective)) +
  geom_histogram(bins=30)
ggsave("PropDefective_Hist.png")

ggplot(data = BOTH_sub_filt, aes(x=Prop_Defective)) +
  facet_wrap(BOTH_sub_filt$Genotype) +
  geom_histogram(bins=30)
ggsave("PropDefective_HistByGenotype.png")


#This is a highly zero inflated distribution and cannot be coerced to anything near normal. I explored looking into different statistics meant to handle the high 0 inflation, but they all seems more complicated than was necessary. I am just going to go with a Kruskal-Wallas
fit3=lm(Prop_Defective ~ Genotype, data = BOTH_sub)
gvmodelfit3 <- gvlma(fit3)
gvmodelfit3

model_2 <- kruskal.test(Prop_Defective ~ Genotype, data = BOTH_sub_filt)
model_2
#This reports that the model is significant 

pairwise.wilcox.test(BOTH_sub_filt$Prop_Defective, BOTH_sub_filt$Genotype)
#This reports that only mutation 1 is significantly different 

Abor_Sig <- data.frame(x=c(1, 2, 3), y=c(0.35, 0.35, 0.35), Sig=c("a", "b", "ab"))

pdf("PropDefective_By_trkin_Genotype.pdf", height=10, width=4)
ggplot() +
  geom_jitter(data=BOTH_sub, aes(x=Genotype, y=Prop_Defective), color= "grey60", height = 0, alpha=0.5) +
  geom_boxplot(data=BOTH_sub, aes(x=Genotype, y=Prop_Defective), alpha=0.5, alpha=0.5, outlier.shape = NA) +
  geom_text(data=Abor_Sig, aes(x=x, y=y, label=Sig), size = 10) +
  scale_x_discrete(labels=c("TRKIN Mutation 2" = "1 -\n2 -", "TRKIN Mutation 1" = "1 -\n2 +", "Functional TRKIN" = "1 +\n2 -")) +
  xlab("Ab10-I trkin Genotype") +
  ylab("Proportion Defective Kernels") +
  geom_vline(xintercept=1.5, color = "grey80", linetype="dashed", size = 0.5) +
  geom_vline(xintercept=2.5, color = "grey80", linetype="dashed", size = 0.5) +
  theme(plot.title = element_text(hjust=0.5, size = 30), plot.subtitle = element_text(hjust=0.5, size = 10), axis.title = element_text(size = 20), axis.text.y = element_text(size = 20), axis.text.x = element_text(size = 20), legend.title = element_text(size = 20), legend.text = element_text(size = 20), legend.position = "right") 
dev.off()
