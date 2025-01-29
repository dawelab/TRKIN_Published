library(readxl)
library(ggplot2)
library(gvlma)

setwd("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R\ Sessions/Paper/TRKIN_Published/Homozygous_Fitness")

DATA <- read_excel("~/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/Trkin_MappingPopulation_v2.xlsx")
DATA$`Height (cm)` <- as.numeric(DATA$`Height (cm)`)
DATA$`Height_Bin` <- as.factor(DATA$`Height_Bin`)
DATA$Pot <- as.factor(DATA$Pot)
DATA$`Position in row (1 is front)` <- as.numeric(DATA$`Position in row (1 is front)`)
DATA$`Row (1 is right, 2 is left)` <- as.factor(DATA$`Row (1 is right, 2 is left)`)
DATA$`Pollen (Half was on 9/10/23) (1 is first half, 2 is second half)`<- as.factor(DATA$`Pollen (Half was on 9/10/23) (1 is first half, 2 is second half)`)
DATA$`Silk (Half was on 9/14/24 1 is first half, 2 is second half)` <- as.factor(DATA$`Silk (Half was on 9/14/24 1 is first half, 2 is second half)`)

ggplot(DATA, aes(`Height (cm)`)) +
  geom_histogram(bins=50)


#This fails on link assumption suggesting that my response is not truly continuous, but only just barley I think it is legitiamate to use this
model <- lm(`Height (cm)` ~ Pot + `Position in row (1 is front)` + Genotype, data=DATA)
gvlma(model)
summary(model)
#Trkin has no effect on plant height

#This doesn't exactly meet assumptions but it is close
model2 <- lm(`total kernels` ~ Pot + `Position in row (1 is front)`+ `Pollen (Half was on 9/10/23) (1 is first half, 2 is second half)`+ DATA$`Silk (Half was on 9/14/24 1 is first half, 2 is second half)` + Genotype, data=DATA)
gvlma(model2)
summary(model2)
#Trkin has no effect on total kernel number, though plants that shed pollen in the second half had fewer kernels 

#This drops all samples with NAs 
DATA_Sub <- subset(DATA, DATA$`total kernels` != 0)

#This almost meets assumptions
model3 <- lm(`avg_kernel weight` ~ Pot + `Position in row (1 is front)`+ `Pollen (Half was on 9/10/23) (1 is first half, 2 is second half)`+ `Silk (Half was on 9/14/24 1 is first half, 2 is second half)` + Genotype, data=DATA_Sub)
gvlma(model3)
summary(model3)
#Trkin has no effect on kernel weight 

#Pollen was scored in halves, so I need to use a logistic regression 
model3 <- glm(`Pollen (Half was on 9/10/23) (1 is first half, 2 is second half)` ~ Pot + `Position in row (1 is front)`+  Genotype, data=DATA_Sub, family="binomial")
summary(model3)
#Nothing is significant

#Silk was scored in halves, so I need to use a logistic regression 
model4 <- glm(`Silk (Half was on 9/14/24 1 is first half, 2 is second half)` ~ Pot + `Position in row (1 is front)`+  Genotype, data=DATA_Sub, family="binomial")
summary(model4)
#Nothing is significant

pdf("trkinGenotype_v_Height.pdf", height=4, width=4)
ggplot(data=DATA, aes(x=Genotype, y=`Height (cm)`)) +
  geom_jitter(color="grey30", alpha=0.5) +
  geom_boxplot(alpha=0.5) +
  scale_x_discrete(labels=c("-2 Homozygous" = "1(-) 2(-)\n1(-) 2(-)", "Het for -2 and No edit" = "1(-) 2(-)\n1(+) 2(-)", "No Edit Homozygous"="1(+) 2(-)\n1(+) 2(-)")) +
  labs(y="Height (cm)", x="Ab10 trkin Genotype") +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 15))
dev.off()

facet_names <- c(
  `-2 Homozygous` = "1 - 2 -\n1 - 2 -",
  `Het for -2 and No edit` = "1 - 2 -\n1 + 2 -",
  `No Edit Homozygous` = "1 + 2 -\n1 + 2 -"
)

pdf("trkinGenotype_v_MaleFloweringTime.pdf", height=4, width=4)
ggplot(data=DATA, aes(x=`Pollen (Half was on 9/10/23) (1 is first half, 2 is second half)`)) +
  geom_bar(stat="count") +
  scale_y_continuous(breaks=c(1:12)) +
  scale_x_discrete(labels=c("1"="1", "2"="2", "No Tassel" = "NA")) +
  facet_wrap(Genotype ~ ., labeller = as_labeller(facet_names), ) +
  labs(x="Male Flowering Time") +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 15), strip.text.x = element_text(size = 18))
dev.off()

pdf("trkinGenotype_v_FemaleFloweringTime.pdf", height=4, width=4)
ggplot(data=DATA, aes(x=`Silk (Half was on 9/14/24 1 is first half, 2 is second half)`)) +
  geom_bar(stat="count") +
  scale_y_continuous(breaks=c(1:14)) +
  scale_x_discrete(labels=c("1"="1", "2"="2")) +
  facet_wrap(Genotype ~ ., labeller = as_labeller(facet_names), ) +
  labs(x="Female Flowering Time") +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 15), strip.text.x = element_text(size = 18))
dev.off()
  
pdf("trkinGenotype_v_TotalKernelCount.pdf", height=4, width=4)
ggplot(data=DATA, aes(x=Genotype, y=`total kernels`)) +
  geom_jitter(color="grey30", alpha=0.5) +
  geom_boxplot(alpha=0.5) +
  scale_x_discrete(labels=c("-2 Homozygous" = "1(-) 2(-)\n1(-) 2(-)", "Het for -2 and No edit" = "1(-) 2(-)\n1(+) 2(-)", "No Edit Homozygous"="1(+) 2(-)\n1(+) 2(-)")) +
  labs(y="Kernel Number", x="Ab10 trkin Genotype") +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 15))
dev.off()

pdf("trkinGenotype_v_AvgKernelWeight.pdf", height=4, width=4)
ggplot(data=DATA, aes(x=Genotype, y=as.numeric(`avg_kernel weight`))) +
  geom_jitter(color="grey30", alpha=0.5) +
  geom_boxplot(alpha=0.5) +
  scale_x_discrete(labels=c("-2 Homozygous" = "1(-) 2(-)\n1(-) 2(-)", "Het for -2 and No edit" = "1(-) 2(-)\n1(+) 2(-)", "No Edit Homozygous"="1(+) 2(-)\n1(+) 2(-)")) +
  labs(y="Average Kernel Weight", x="Ab10 trkin Genotype") +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 15))
dev.off()


Expected <- c(rep(0, 10), rep(1, 19), rep(2, 10))
Actual <- DATA$`Functional trkin`

EX <- as.data.frame(table(Expected))
colnames(EX) <- c("copies", "expected")

AC <-  as.data.frame(table(Actual))
colnames(AC) <- c("copies", "observed")

TABLE <- merge(EX, AC)
TABLE <- TABLE[,-c(1)]

chisq.test(TABLE)

pdf("trkinGenotype_BarPlot.pdf", height=4, width=4) 
ggplot(data=DATA, aes(x=Genotype)) +
  geom_bar(stat="count")+
  scale_x_discrete(labels=c("-2 Homozygous" = "1(-) 2(-)\n1(-) 2(-)", "Het for -2 and No edit" = "1(-) 2(-)\n1(+) 2(-)", "No Edit Homozygous"="1(+) 2(-)\n1(+) 2(-)")) +
  labs(x="Ab10 trkin Genotype") +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 15), strip.text.x = element_text(size = 18))
dev.off()

#Trkin does not affect Ab10 segregartion 
