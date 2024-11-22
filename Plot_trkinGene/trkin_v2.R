library(readxl)
library(ggplot2)
library(gvlma)

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

ggplot(data=DATA, aes(x=Genotype, y=`Height (cm)`)) +
  geom_boxplot() +
  geom_point() +
  ggtitle("Effect of trkin Genotype on Height") +
  labs(subtitle = "LM Height ~ Pot + Position + Genotype\n NS") +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5))
ggsave("trkinGenotype_v_Height.png")


ggplot(data=DATA, aes(x=`Pollen (Half was on 9/10/23) (1 is first half, 2 is second half)`)) +
  geom_bar(stat="count") +
  scale_y_continuous(breaks=c(1:12)) +
  facet_wrap(facet=~Genotype) +
  ggtitle("Effect of trkin Genotype on Male Flowering Time") +
  labs(subtitle = "GLM (family=binomial) Male Flowering Time ~ Pot + Position + Genotype\n NS") +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5)) +
  xlab("Male Flowering Time\n1=1st Half, 2=2nd Half")
ggsave("trkinGenotype_v_MaleFloweringTime.png")
  
ggplot(data=DATA, aes(x=`Silk (Half was on 9/14/24 1 is first half, 2 is second half)`)) +
  geom_bar(stat="count") +
  scale_y_continuous(breaks=c(1:14)) +
  facet_wrap(facet=~Genotype) +
  ggtitle("Effect of trkin Genotype on Female Flowering Time") +
  labs(subtitle = "GLM (family=binomial) Female Flowering Time ~ Pot + Position + Genotype\n NS") +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5)) +
  xlab("Female Flowering Time\n1=1st Half, 2=2nd Half")
ggsave("trkinGenotype_v_FemaleFloweringTime.png")

ggplot(data=DATA, aes(x=Genotype, y=`total kernels`)) +
  geom_boxplot() +
  geom_point() +
  ggtitle("Effect of trkin Genotype on Kernel Count") +
  labs(subtitle = "LM Total Kernel Count ~ Pot + Position + Male Flowering Time +\n Female Flowering Time + Genotype\n NS") +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5))
ggsave("trkinGenotype_v_TotalKernelCount.png")

ggplot(data=DATA, aes(x=Genotype, y=as.numeric(`avg_kernel weight`))) +
  geom_boxplot() +
  geom_point() +
  ggtitle("Effect of trkin Genotype on Avg Kernel Weight") +
  labs(subtitle = "LM Avg Kernel Weight ~ Pot + Position + Male Flowering Time +\n Female Flowering Time + Genotype\n NS") +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5))
ggsave("trkinGenotype_v_AvgKernelWeight.png")


Expected <- c(rep(0, 10), rep(1, 19), rep(2, 10))
Actual <- DATA$`Functional trkin`

EX <- as.data.frame(table(Expected))
colnames(EX) <- c("copies", "expected")

AC <-  as.data.frame(table(Actual))
colnames(AC) <- c("copies", "observed")

TABLE <- merge(EX, AC)
TABLE <- TABLE[,-c(1)]

chisq.test(TABLE)

ggplot(data=DATA, aes(x=Genotype)) +
  geom_bar(stat="count")+
  ggtitle("Ab10 Homozygous trkin Segregation") +
  labs(subtitle = "Chisquared=1.26, df=2, p=0.0.53") +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5))
ggsave("trkinGenotype_BarPlot.png")

#Trkin does not affect Ab10 segregartion 