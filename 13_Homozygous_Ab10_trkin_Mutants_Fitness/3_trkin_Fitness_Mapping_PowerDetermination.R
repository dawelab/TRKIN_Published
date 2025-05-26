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

model <- lm(`Height (cm)` ~ Pot + `Position in row (1 is front)` + Genotype, data=DATA)

#This determines the total number of plants per class I had
Count <- summary(as.factor(DATA$Genotype))

#Calculate the mean and standard deviation
Ab10_mean <- mean(DATA$`Height (cm)`)
Ab10_sd <- sd(DATA$`Height (cm)`)

#Height need to simulate normal 
i=1
n=500
diff= c(1:60)
OUT <- data.frame("Iteration"=c(NA), "NoEdit vs HomMut" = c(NA), "NoEdit vs HetMut"= c(NA), "HomMut vs HetMut"= c(NA))
OUT <- OUT[-c(1),]
FIN <- data.frame("Difference"=c(NA), "Power" = c(NA))
FIN <- FIN[-c(1),]

for(j in diff) {
  print(j)
  for(i in 1:n) {
    HomMut <- rnorm(Count[1], mean = (Ab10_mean + j), sd = Ab10_sd)
    HetMut <- rnorm(Count[2], mean = (Ab10_mean + j/2), sd = Ab10_sd)
    NoEdit <- rnorm(Count[3], mean = Ab10_mean, sd = Ab10_sd)
    SimString <- c(HomMut, HetMut, NoEdit)
    Genotype <- c(rep("HomMut", Count[1]), rep("HetMut", Count[2]), rep("NoEdit", Count[3]))
    SimDF <- cbind(SimString, Genotype)
    SimDF <- as.data.frame(SimDF, row.names = NULL, optional = FALSE, make.names = TRUE, stringsAsFactors = FALSE)
    model_1 <- aov(SimString ~ Genotype, data = SimDF)
    Tukey <- TukeyHSD(model_1, which="Genotype")
    Table <- as.data.frame(Tukey$Genotype)
    INDEX1 <- grep("NoEdit-HomMut", rownames(Table))
    p1 <- Table[INDEX1, 4]
    Sig1 <- ifelse(p1 <= 0.05, "yes", "no")
    INDEX2 <- grep("NoEdit-HetMut", rownames(Table))
    p2 <- Table[INDEX2, 4]
    Sig2 <- ifelse(p2 <= 0.05, "yes", "no")
    INDEX3 <- grep("HomMut-HetMut", rownames(Table))
    p3 <- Table[INDEX3, 4]
    Sig3 <- ifelse(p3 <= 0.05, "yes", "no")
    ROW <- cbind(i, Sig1, Sig2, Sig3)
    colnames(ROW) <- c("Iteration","NoEdit-HomMut","NoEdit-HetMut","HomMut-HetMut")
    OUT <<- rbind (OUT,ROW)
  }
  TEST <- subset(OUT, `NoEdit-HomMut` == "yes")
  Power <- nrow(TEST)/n
  ROW2 <- cbind(j, Power)
  colnames(ROW2) <- c("Difference", "Power")
  FIN <<- rbind(FIN, ROW2)
  OUT <- data.frame("Iteration"=c(NA), "NoEdit-HomMut" = c(NA), "NoEdit-HetMut"= c(NA), "HomMut-HetMut"= c(NA))
  OUT <- OUT[c(1),]
  OUT <- OUT[-c(1),]
}


#Total Kernel Numbers

#Calculate the mean and standard deviation
Ab10_mean <- mean(DATA$`total kernels`)
Ab10_sd <- sd(DATA$`total kernels`)


i=1
n=500
diff= c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100, 115, 120, 125, 130)
OUT <- data.frame("Iteration"=c(NA), "NoEdit vs HomMut" = c(NA), "NoEdit vs HetMut"= c(NA), "HomMut vs HetMut"= c(NA))
OUT <- OUT[-c(1),]
FIN <- data.frame("Difference"=c(NA), "Power" = c(NA))
FIN <- FIN[-c(1),]

for(j in diff) {
  print(j)
  for(i in 1:n) {
    HomMut <- rnorm(Count[1]*100, mean = (Ab10_mean + j), sd = Ab10_sd)
    HetMut <- rnorm(Count[2]*100, mean = (Ab10_mean + j/2), sd = Ab10_sd)
    NoEdit <- rnorm(Count[3]*100, mean = Ab10_mean, sd = Ab10_sd)
    SimString <- c(HomMut, HetMut, NoEdit)
    Genotype <- c(rep("HomMut", Count[1]), rep("HetMut", Count[2]), rep("NoEdit", Count[3]))
    SimDF <- cbind(SimString, Genotype)
    SimDF <- as.data.frame(SimDF, row.names = NULL, optional = FALSE, make.names = TRUE, stringsAsFactors = FALSE)
    model_1 <- aov(SimString ~ Genotype, data = SimDF)
    Tukey <- TukeyHSD(model_1, which="Genotype")
    Table <- as.data.frame(Tukey$Genotype)
    INDEX1 <- grep("NoEdit-HomMut", rownames(Table))
    p1 <- Table[INDEX1, 4]
    Sig1 <- ifelse(p1 <= 0.05, "yes", "no")
    INDEX2 <- grep("NoEdit-HetMut", rownames(Table))
    p2 <- Table[INDEX2, 4]
    Sig2 <- ifelse(p2 <= 0.05, "yes", "no")
    INDEX3 <- grep("HomMut-HetMut", rownames(Table))
    p3 <- Table[INDEX3, 4]
    Sig3 <- ifelse(p3 <= 0.05, "yes", "no")
    ROW <- cbind(i, Sig1, Sig2, Sig3)
    colnames(ROW) <- c("Iteration","NoEdit-HomMut","NoEdit-HetMut","HomMut-HetMut")
    OUT <<- rbind (OUT,ROW)
  }
  TEST <- subset(OUT, `NoEdit-HomMut` == "yes")
  Power <- nrow(TEST)/n
  ROW2 <- cbind(j, Power)
  colnames(ROW2) <- c("Difference", "Power")
  FIN <<- rbind(FIN, ROW2)
  OUT <- data.frame("Iteration"=c(NA), "NoEdit-HomMut" = c(NA), "NoEdit-HetMut"= c(NA), "HomMut-HetMut"= c(NA))
  OUT <- OUT[c(1),]
  OUT <- OUT[-c(1),]
}

#Kernel Weight

#Calculate the mean and standard deviation
DATA$`avg_kernel weight` <- as.numeric(DATA$`avg_kernel weight`)
DATA_SUB <- subset(DATA, is.na(`avg_kernel weight`) == FALSE)

Ab10_mean <- mean(DATA_SUB$`avg_kernel weight`)
Ab10_sd <- sd(DATA_SUB$`avg_kernel weight`)

#Kernel Weight need to simulate normal 
i=1
n=500
diff= c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2)
OUT <- data.frame("Iteration"=c(NA), "NoEdit vs HomMut" = c(NA), "NoEdit vs HetMut"= c(NA), "HomMut vs HetMut"= c(NA))
OUT <- OUT[-c(1),]
FIN <- data.frame("Difference"=c(NA), "Power" = c(NA))
FIN <- FIN[-c(1),]

for(j in diff) {
  print(j)
  for(i in 1:n) {
    HomMut <- rnorm(Count[1], mean = (Ab10_mean + j), sd = Ab10_sd)
    HetMut <- rnorm(Count[2], mean = (Ab10_mean + j/2), sd = Ab10_sd)
    NoEdit <- rnorm(Count[3], mean = Ab10_mean, sd = Ab10_sd)
    SimString <- c(HomMut, HetMut, NoEdit)
    Genotype <- c(rep("HomMut", Count[1]), rep("HetMut", Count[2]), rep("NoEdit", Count[3]))
    SimDF <- cbind(SimString, Genotype)
    SimDF <- as.data.frame(SimDF, row.names = NULL, optional = FALSE, make.names = TRUE, stringsAsFactors = FALSE)
    model_1 <- aov(SimString ~ Genotype, data = SimDF)
    Tukey <- TukeyHSD(model_1, which="Genotype")
    Table <- as.data.frame(Tukey$Genotype)
    INDEX1 <- grep("NoEdit-HomMut", rownames(Table))
    p1 <- Table[INDEX1, 4]
    Sig1 <- ifelse(p1 <= 0.05, "yes", "no")
    INDEX2 <- grep("NoEdit-HetMut", rownames(Table))
    p2 <- Table[INDEX2, 4]
    Sig2 <- ifelse(p2 <= 0.05, "yes", "no")
    INDEX3 <- grep("HomMut-HetMut", rownames(Table))
    p3 <- Table[INDEX3, 4]
    Sig3 <- ifelse(p3 <= 0.05, "yes", "no")
    ROW <- cbind(i, Sig1, Sig2, Sig3)
    colnames(ROW) <- c("Iteration","NoEdit-HomMut","NoEdit-HetMut","HomMut-HetMut")
    OUT <<- rbind (OUT,ROW)
  }
  TEST <- subset(OUT, `NoEdit-HomMut` == "yes")
  Power <- nrow(TEST)/n
  ROW2 <- cbind(j, Power)
  colnames(ROW2) <- c("Difference", "Power")
  FIN <<- rbind(FIN, ROW2)
  OUT <- data.frame("Iteration"=c(NA), "NoEdit-HomMut" = c(NA), "NoEdit-HetMut"= c(NA), "HomMut-HetMut"= c(NA))
  OUT <- OUT[c(1),]
  OUT <- OUT[-c(1),]
}


#########################################

#Male Flowering Time

#Calculate the mean and standard deviation
DATA$`Pollen (Half was on 9/10/23) (1 is first half, 2 is second half)` <- as.numeric(DATA$`Pollen (Half was on 9/10/23) (1 is first half, 2 is second half)`)
DATA_SUB <- subset(DATA, is.na(`Pollen (Half was on 9/10/23) (1 is first half, 2 is second half)`) == FALSE & Genotype == "No Edit Homozygous" & `Pollen (Half was on 9/10/23) (1 is first half, 2 is second half)` != "No Tassel")
Ab10_probdf <- as.data.frame(summary(as.factor(DATA_SUB$`Pollen (Half was on 9/10/23) (1 is first half, 2 is second half)`)))
Ab10_prob <- Ab10_probdf[1,1]/(Ab10_probdf[1,1] + Ab10_probdf[2,1])


i=1
n=500
diff= c(-0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4)
OUT <- data.frame("Iteration"=c(NA), "NoEdit vs HomMut" = c(NA), "NoEdit vs HetMut"= c(NA), "HomMut vs HetMut"= c(NA))
OUT <- OUT[-c(1),]
FIN <- data.frame("Difference"=c(NA), "Power" = c(NA))
FIN <- FIN[-c(1),]

for(j in diff) {
  print(j)
  for(i in 1:n) {
    HomMut <- rbinom(Count[1], 1, prob = Ab10_prob+j)
    HetMut <- rbinom(Count[2], 1, prob = Ab10_prob+j/2)
    NoEdit <- rbinom(Count[3], 1, prob = Ab10_prob)
    SimString <- c(HomMut, HetMut, NoEdit)
    Genotype <- c(rep("HomMut", Count[1]), rep("HetMut", Count[2]), rep("NoEdit", Count[3]))
    SimDF <- cbind(SimString, Genotype)
    SimDF <- as.data.frame(SimDF, row.names = NULL, optional = FALSE, make.names = TRUE, stringsAsFactors = FALSE)
    SimDF$SimString <- as.factor(SimDF$SimString)
    SimDF$Genotype <- factor(x = SimDF$Genotype, c("NoEdit","HetMut","HomMut"))
    model3 <- glm(SimString ~ Genotype, data= SimDF, family="binomial")
    Summary <- summary(model3)
    Table <- Summary$coefficients
    INDEX1 <- grep("HetMut", rownames(Table))
    p1 <- Table[INDEX1, 4]
    Sig1 <- ifelse(p1 <= 0.05, "yes", "no")
    INDEX2 <- grep("HomMut", rownames(Table))
    p2 <- Table[INDEX2, 4]
    Sig2 <- ifelse(p2 <= 0.05, "yes", "no")
    ROW <- cbind(i, Sig1, Sig2)
    colnames(ROW) <- c("Iteration","HetMut","HomMut")
    OUT <<- rbind (OUT,ROW)
  }
  TEST <- subset(OUT, `HomMut` == "yes")
  Power <- nrow(TEST)/n
  ROW2 <- cbind(j, Power)
  colnames(ROW2) <- c("Difference", "Power")
  FIN <<- rbind(FIN, ROW2)
  OUT <- data.frame("Iteration"=c(NA), "NoEdit-HomMut" = c(NA), "NoEdit-HetMut"= c(NA), "HomMut-HetMut"= c(NA))
  OUT <- OUT[c(1),]
  OUT <- OUT[-c(1),]
}


#########################################

#Female Flowering Time

#Calculate the mean and standard deviation
DATA$`Silk (Half was on 9/14/24 1 is first half, 2 is second half)` <- as.numeric(DATA$`Silk (Half was on 9/14/24 1 is first half, 2 is second half)`)
DATA_SUB <- subset(DATA, is.na(`Silk (Half was on 9/14/24 1 is first half, 2 is second half)`) == FALSE & Genotype == "No Edit Homozygous")
Ab10_probdf <- as.data.frame(summary(as.factor(DATA_SUB$`Silk (Half was on 9/14/24 1 is first half, 2 is second half)`)))
Ab10_prob <- Ab10_probdf[2,1]/Ab10_probdf[1,1]


i=1
n=500
diff= c(-0.8,-0.7,-0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.1)
OUT <- data.frame("Iteration"=c(NA), "NoEdit vs HomMut" = c(NA), "NoEdit vs HetMut"= c(NA), "HomMut vs HetMut"= c(NA))
OUT <- OUT[-c(1),]
FIN <- data.frame("Difference"=c(NA), "Power" = c(NA))
FIN <- FIN[-c(1),]

for(j in diff) {
  print(j)
  for(i in 1:n) {
    HomMut <- rbinom(Count[1], 1, prob = Ab10_prob+j)
    HetMut <- rbinom(Count[2], 1, prob = Ab10_prob+j/2)
    NoEdit <- rbinom(Count[3], 1, prob = Ab10_prob)
    SimString <- c(HomMut, HetMut, NoEdit)
    Genotype <- c(rep("HomMut", Count[1]), rep("HetMut", Count[2]), rep("NoEdit", Count[3]))
    SimDF <- cbind(SimString, Genotype)
    SimDF <- as.data.frame(SimDF, row.names = NULL, optional = FALSE, make.names = TRUE, stringsAsFactors = FALSE)
    SimDF$SimString <- as.factor(SimDF$SimString)
    SimDF$Genotype <- factor(x = SimDF$Genotype, c("NoEdit","HetMut","HomMut"))
    model3 <- glm(SimString ~ Genotype, data= SimDF, family="binomial")
    Summary <- summary(model3)
    Table <- Summary$coefficients
    INDEX1 <- grep("HetMut", rownames(Table))
    p1 <- Table[INDEX1, 4]
    Sig1 <- ifelse(p1 <= 0.05, "yes", "no")
    INDEX2 <- grep("HomMut", rownames(Table))
    p2 <- Table[INDEX2, 4]
    Sig2 <- ifelse(p2 <= 0.05, "yes", "no")
    ROW <- cbind(i, Sig1, Sig2)
    colnames(ROW) <- c("Iteration","HetMut","HomMut")
    OUT <<- rbind (OUT,ROW)
  }
  TEST <- subset(OUT, `HomMut` == "yes")
  Power <- nrow(TEST)/n
  ROW2 <- cbind(j, Power)
  colnames(ROW2) <- c("Difference", "Power")
  FIN <<- rbind(FIN, ROW2)
  OUT <- data.frame("Iteration"=c(NA), "NoEdit-HomMut" = c(NA), "NoEdit-HetMut"= c(NA), "HomMut-HetMut"= c(NA))
  OUT <- OUT[c(1),]
  OUT <- OUT[-c(1),]
}

#I don't have enough data to detect any difference in flowering time or Kernel number. 

#########################################

#Between Ab10 competitiveness

#Calculating competition between Ab10s 

Expected <- c(rep(0, round(sum(Count)*.25)), rep(1, round(sum(Count)*.5)), rep(2, round(sum(Count)*.25)))
Sim <- c(rep(0, round(sum(Count)*.25)), rep(1, round(sum(Count)*.5)), rep(2, round(sum(Count)*.25)))



FIN <- data.frame("Difference"=c(NA), "Power" = c(NA))
FIN <- FIN[-c(1),]

  j="10%"
  Sim <- c(rep(0, round(n*0.16)), rep(1, round(n*0.48)), rep(2, round(n*0.36)))
  j="20%"
  Sim <- c(rep(0, round(n*0.09)), rep(1, round(n*0.42)), rep(2, round(n*0.49)))
  j="21%"
  Sim <- c(rep(0, round(n*0.0841)), rep(1, round(n*0.4118)), rep(2, round(n*0.5041)))
  j="23%"
  Sim <- c(rep(0, round(n*0.0729)), rep(1, round(n*0.3942)), rep(2, round(n*0.5329)))
  j="25%"
  Sim <- c(rep(0, round(n*0.0625)), rep(1, round(n*0.375)), rep(2, round(n*0.5625)))
  j="30%"
  Sim <- c(rep(0, round(n*0.01)), rep(1, round(n*0.18)), rep(2, round(n*0.81)))
  n=sum(Count)
  
  i=1
  n=149
  OUT <- data.frame("Iteration"=c(NA), "Different" = c(NA))
  OUT <- OUT[-c(1),]
  
  for(i in 1:5) {
    Expected <- c(rep(0, round(n*.25)), rep(1, round(n*.5)), rep(2, round(n*.25)))
    Sim <- c(rep(0, round(n*0.16)), rep(1, round(n*0.48)), rep(2, round(n*0.36)))
    EX <- as.data.frame(table(Expected))
    colnames(EX) <- c("copies", "expected")
    AC <-  as.data.frame(table(Sim))
    colnames(AC) <- c("copies", "observed")
    TABLE <- merge(EX, AC)
    TABLE <- TABLE[,-c(1)]
    test <- chisq.test(TABLE)
    p <- test$p.value
    Sig <- ifelse(p <= 0.05, "yes", "no")
    ROW <- cbind(i, Sig)
    colnames(ROW) <- c("Iteration", "Different")
    OUT <<- rbind (OUT,ROW)
  }
  
  
  TEST <- subset(OUT, Different == "yes")
  Power <- nrow(TEST)/n
  ROW2 <- cbind(j, Power)
  colnames(ROW2) <- c("Difference", "Power")
  FIN <<- rbind(FIN, ROW2)
  OUT <- data.frame("Iteration"=c(NA), "Different" = c(NA))
  OUT <- OUT[c(1),]
  OUT <- OUT[-c(1),]
