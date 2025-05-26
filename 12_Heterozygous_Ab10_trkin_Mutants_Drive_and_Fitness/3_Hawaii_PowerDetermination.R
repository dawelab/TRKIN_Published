library(readxl)
library(tidyverse)

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

#Filters to only the Ab10 trkin1+ 
BOTH_SUB <- subset(BOTH, Genotype == "Functional TRKIN")

#I am now getting the mean and standard deviation for both 
Ab10_mean <- mean(BOTH_SUB$Prop_Defective)
Ab10_sd <- sd(BOTH_SUB$Prop_Defective)


#This determines the total number of plants per class I had
Count <- summary(as.factor(BOTH$Genotype))

#Drive need to simulate normal 
i=1
n=500
diff= c(0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.017, 0.018, 0.019, 0.02)
OUT <- data.frame("Iteration"=c(NA), "Functional vs Mut1" = c(NA), "Functional vs Mut2"= c(NA), "Mut1 vs Mut2"= c(NA))
OUT <- OUT[-c(1),]
FIN <- data.frame("Difference"=c(NA), "Power" = c(NA))
FIN <- FIN[-c(1),]

for(j in diff) {
  print(j)
  for(i in 1:n) {
    Ab10_sim <- rnorm(Count[1], mean = Ab10_mean, sd = Ab10_sd)
    trkinmut1_sim <- rnorm(Count[2], mean = (Ab10_mean + j), sd = Ab10_sd)
    trkinmut2_sim <- rnorm(Count[3], mean = (Ab10_mean + j), sd = Ab10_sd)
    Drive <- c(Ab10_sim, trkinmut1_sim, trkinmut2_sim)
    Genotype <- c(rep("trkin", Count[1]), rep("mut1", Count[2]), rep("mut2", Count[3]))
    SimDrive <- cbind(Drive, Genotype)
    SimDrive <- as.data.frame(SimDrive, row.names = NULL, optional = FALSE, make.names = TRUE, stringsAsFactors = FALSE)
    model_1 <- aov(Drive~ Genotype, data = SimDrive)
    summary(model_1)
    Tukey <- TukeyHSD(model_1, which="Genotype")
    Table <- as.data.frame(Tukey$Genotype)
    INDEX1 <- grep("trkin-mut1", rownames(Table))
    p1 <- Table[INDEX1, 4]
    Sig1 <- ifelse(p1 <= 0.05, "yes", "no")
    INDEX2 <- grep("trkin-mut2", rownames(Table))
    p2 <- Table[INDEX2, 4]
    Sig2 <- ifelse(p2 <= 0.05, "yes", "no")
    INDEX3 <- grep("mut2-mut1", rownames(Table))
    p3 <- Table[INDEX3, 4]
    Sig3 <- ifelse(p3 <= 0.05, "yes", "no")
    ROW <- cbind(i, Sig1, Sig2, Sig3)
    colnames(ROW) <- c("Iteration","Functional vs Mut1","Functional vs Mut2","Mut1 vs Mut2")
    OUT <<- rbind (OUT,ROW)
  }
  colnames(OUT) <- c("Iteration","Functional vs Mut1","Functional vs Mut2","Mut1 vs Mut2")
  TEST <- subset(OUT, `Functional vs Mut1` == "yes" &  `Functional vs Mut2` == "yes" & `Mut1 vs Mut2` == "no")
  Power <- nrow(TEST)/n
  ROW2 <- cbind(j, Power)
  colnames(ROW2) <- c("Difference", "Power")
  FIN <<- rbind(FIN, ROW2)
  OUT <- data.frame("Iteration"=c(NA), "Functional vs Mut1" = c(NA), "Functional vs Mut2"= c(NA), "Mut1 vs Mut2"= c(NA))
  OUT <- OUT[c(1),]
  OUT <- OUT[-c(1),]
}


FIN$Percent_Change <- FIN$Difference*100
FIN$Percent_Power <- FIN$Power*100

################Kernel Abortion Need to simulate 0 inflated data set
#I am now getting the mean and standard deviation for both 
#I am now getting the mean and standard deviation for both 
#This removes NAs
BOTH_SUBNA <- subset(BOTH_SUB, is.na(Prop_Defective) == FALSE)
Ab10_mean <- mean(BOTH_SUBNA$Prop_Defective)
Ab10_sd <- sd(BOTH_SUBNA$Prop_Defective)

i=1
n=500
diff= c(0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.017, 0.018, 0.019, 0.02)
OUT <- data.frame("Iteration"=c(NA), "Functional vs Mut1" = c(NA), "Functional vs Mut2"= c(NA), "Mut1 vs Mut2"= c(NA))
OUT <- OUT[-c(1),]
FIN <- data.frame("Difference"=c(NA), "Power" = c(NA))
FIN <- FIN[-c(1),]

for(j in diff) {
  print(j)
  for(i in 1:n) {
    Ab10_sim <- rnorm(Count[1], mean = Ab10_mean, sd = Ab10_sd)
    trkinmut1_sim <- rnorm(Count[2], mean = (Ab10_mean + j), sd = Ab10_sd)
    trkinmut2_sim <- rnorm(Count[3], mean = (Ab10_mean + j), sd = Ab10_sd)
    Defect <- c(Ab10_sim, trkinmut1_sim, trkinmut2_sim)
    Genotype <- c(rep("trkin", Count[1]), rep("mut1", Count[2]), rep("mut2", Count[3]))
    SimDefect <- cbind( Defect, Genotype)
    SimDefect <- as.data.frame(SimDefect, row.names = NULL, optional = FALSE, make.names = TRUE, stringsAsFactors = FALSE)
    SimDefect$Defect <- as.numeric(SimDefect$Defect)
 
    model_2 <- kruskal.test( Defect ~ Genotype, data = SimDefect)
    wilcox <- pairwise.wilcox.test(SimDefect$Defect, SimDefect$Genotype)
    Table <- as.data.frame(wilcox$p.value)
    INDEX1 <- grep("mut2", rownames(Table))
    INDEX2 <- grep("trkin", rownames(Table))
    INDEX3 <- grep("mut1", colnames(Table))
    INDEX4 <- grep("mut2", colnames(Table))
    p1 <- Table[INDEX2, INDEX3]
    Sig1 <- ifelse(p1 <= 0.05, "yes", "no")
    p2 <- Table[INDEX2, INDEX4]
    Sig2 <- ifelse(p2 <= 0.05, "yes", "no")
    p3 <- Table[INDEX1, INDEX3]
    Sig3 <- ifelse(p3 <= 0.05, "yes", "no")
    ROW <- cbind(i, Sig1, Sig2, Sig3)
    colnames(ROW) <- c("Iteration","Functional vs Mut1","Functional vs Mut2","Mut1 vs Mut2")
    OUT <<- rbind (OUT,ROW)
  }
  colnames(OUT) <- c("Iteration","Functional vs Mut1","Functional vs Mut2","Mut1 vs Mut2")
  TEST <- subset(OUT, `Functional vs Mut1` == "yes" &  `Functional vs Mut2` == "yes" & `Mut1 vs Mut2` == "no")
  Power <- nrow(TEST)/n
  ROW2 <- cbind(j, Power)
  colnames(ROW2) <- c("Difference", "Power")
  FIN <<- rbind(FIN, ROW2)
  OUT <- data.frame("Iteration"=c(NA), "Functional vs Mut1" = c(NA), "Functional vs Mut2"= c(NA), "Mut1 vs Mut2"= c(NA))
  OUT <- OUT[c(1),]
  OUT <- OUT[-c(1),]
}


FIN$Percent_Change <- FIN$Difference*100
FIN$Percent_Power <- FIN$Power*100




#################Total need to simulate normal distribution

#I am now getting the mean and standard deviation for both 
Ab10_mean <- mean(BOTH_SUB$Total)
Ab10_sd <- sd(BOTH_SUB$Total)


#Drive need to simulate normal 
i=1
n=500
diff= c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100)
OUT <- data.frame("Iteration"=c(NA), "Functional vs Mut1" = c(NA), "Functional vs Mut2"= c(NA), "Mut1 vs Mut2"= c(NA))
OUT <- OUT[-c(1),]
FIN <- data.frame("Difference"=c(NA), "Power" = c(NA))
FIN <- FIN[-c(1),]

for(j in diff) {
  print(j)
  for(i in 1:n) {
    Ab10_sim <- rnorm(Count[1], mean = Ab10_mean, sd = Ab10_sd)
    trkinmut1_sim <- rnorm(Count[2], mean = (Ab10_mean + j), sd = Ab10_sd)
    trkinmut2_sim <- rnorm(Count[3], mean = (Ab10_mean + j), sd = Ab10_sd)
    Total <- c(Ab10_sim, trkinmut1_sim, trkinmut2_sim)
    Genotype <- c(rep("trkin", Count[1]), rep("mut1", Count[2]), rep("mut2", Count[3]))
    SimTotal <- cbind(Total, Genotype)
    SimTotal <- as.data.frame(SimTotal, row.names = NULL, optional = FALSE, make.names = TRUE, stringsAsFactors = FALSE)
    model_1 <- aov(Total ~ Genotype, data = SimTotal)
    summary(model_1)
    Tukey <- TukeyHSD(model_1, which="Genotype")
    Table <- as.data.frame(Tukey$Genotype)
    INDEX1 <- grep("trkin-mut1", rownames(Table))
    p1 <- Table[INDEX1, 4]
    Sig1 <- ifelse(p1 <= 0.05, "yes", "no")
    INDEX2 <- grep("trkin-mut2", rownames(Table))
    p2 <- Table[INDEX2, 4]
    Sig2 <- ifelse(p2 <= 0.05, "yes", "no")
    INDEX3 <- grep("mut2-mut1", rownames(Table))
    p3 <- Table[INDEX3, 4]
    Sig3 <- ifelse(p3 <= 0.05, "yes", "no")
    ROW <- cbind(i, Sig1, Sig2, Sig3)
    colnames(ROW) <- c("Iteration","Functional vs Mut1","Functional vs Mut2","Mut1 vs Mut2")
    OUT <<- rbind (OUT,ROW)
  }
  colnames(OUT) <- c("Iteration","Functional vs Mut1","Functional vs Mut2","Mut1 vs Mut2")
  TEST <- subset(OUT, `Functional vs Mut1` == "yes" &  `Functional vs Mut2` == "yes" & `Mut1 vs Mut2` == "no")
  Power <- nrow(TEST)/n
  ROW2 <- cbind(j, Power)
  colnames(ROW2) <- c("Difference", "Power")
  FIN <<- rbind(FIN, ROW2)
  OUT <- data.frame("Iteration"=c(NA), "Functional vs Mut1" = c(NA), "Functional vs Mut2"= c(NA), "Mut1 vs Mut2"= c(NA))
  OUT <- OUT[c(1),]
  OUT <- OUT[-c(1),]
}


FIN$Percent_Change <- FIN$Difference*100
FIN$Percent_Power <- FIN$Power*100
