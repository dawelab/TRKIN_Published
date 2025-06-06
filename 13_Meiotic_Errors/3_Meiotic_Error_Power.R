
DATA <- read_excel("~/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/Micronuclei_Count_format2.xlsx")

DATA <- subset(DATA, Error_Type != "Minicell/Tetrad Cell")

DATA_1 <- subset(DATA, Error_Type == "% Tetrad Micronuclei" & Allele == "N10")
N10_mean <- mean(DATA_1$Perc_Error)
N10_sd <- sd(DATA_1$Perc_Error)
DATA_1 <- subset(DATA, Error_Type == "% Tetrad Micronuclei")
Count <- summary(as.factor(DATA_1$Allele))

#Looking at power
i=1
n=100
diff= c(1:50)
OUT <- data.frame("Iteration"=c(NA), "N10" = c(NA), "trkin_functional"= c(NA), "trkin_mut1"= c(NA), "trkin_mut2" = c(NA))
OUT <- OUT[-c(1),]
FIN <- data.frame("Difference"=c(NA), "Power" = c(NA))
FIN <- FIN[-c(1),]

for(j in diff) {
  print(j)
  for(i in 1:n) {
    N10 <- rnorm(Count[1], mean = (N10_mean), sd = N10_sd)
    func <- rnorm(Count[2], mean = (N10_mean + j), sd = N10_sd)
    mut1 <- rnorm(Count[3], mean = (N10_mean + j/2), sd = N10_sd)
    mut2 <- rnorm(Count[4], mean = (N10_mean + j/2), sd = N10_sd)
    SimString <- c(N10, func, mut1, mut2)
    Genotype <- c(rep("N10", Count[1]), rep("func", Count[2]), rep("mut1", Count[3]), rep("mut2", Count[4]))
    SimDF <- cbind(SimString, Genotype)
    SimDF <- as.data.frame(SimDF, row.names = NULL, optional = FALSE, make.names = TRUE, stringsAsFactors = FALSE)
    model_1 <- aov(SimString ~ Genotype, data = SimDF)
    Tukey <- TukeyHSD(model_1, which="Genotype")
    Table <- as.data.frame(Tukey$Genotype)
    INDEX1 <- grep("N10-func", rownames(Table))
    p1 <- Table[INDEX1, 4]
    Sig1 <- ifelse(p1 <= 0.05, "yes", "no")
    INDEX2 <- grep("mut2-func", rownames(Table))
    p2 <- Table[INDEX2, 4]
    Sig2 <- ifelse(p2 <= 0.05, "yes", "no")
    ROW <- cbind(i, Sig1, Sig2)
    colnames(ROW) <- c("Iteration","N10-func","func-mut2")
    OUT <<- rbind (OUT,ROW)
  }
  TEST <- subset(OUT, `N10-func` == "yes", `func-mut2` == "yes")
  Power <- nrow(TEST)/n
  ROW2 <- cbind(j, Power)
  colnames(ROW2) <- c("Difference", "Power")
  FIN <<- rbind(FIN, ROW2)
  OUT <- data.frame("Iteration"=c(NA), "N10" = c(NA), "trkin_functional"= c(NA), "trkin_mut1"= c(NA), "trkin_mut2" = c(NA))
  OUT <- OUT[c(1),]
  OUT <- OUT[-c(1),]
}


DATA_2 <- subset(DATA, Error_Type == "% Tetrad Microcyte" & Allele == "N10")
N10_mean <- mean(DATA_2$Perc_Error)
N10_sd <- sd(DATA_2$Perc_Error)
DATA_2 <- subset(DATA, Error_Type == "% Tetrad Micronuclei")
Count <- summary(as.factor(DATA_2$Allele))

#Looking at power
i=1
n=500
diff= c(0.000000000000000000001)
OUT <- data.frame("Iteration"=c(NA), "N10" = c(NA), "trkin_functional"= c(NA), "trkin_mut1"= c(NA), "trkin_mut2" = c(NA))
OUT <- OUT[-c(1),]
FIN <- data.frame("Difference"=c(NA), "Power" = c(NA))
FIN <- FIN[-c(1),]

for(j in diff) {
  print(j)
  for(i in 1:n) {
    N10 <- rnorm(Count[1], mean = (N10_mean), sd = N10_sd)
    func <- rnorm(Count[2], mean = (N10_mean + j), sd = N10_sd)
    mut1 <- rnorm(Count[3], mean = (N10_mean + j/2), sd = N10_sd)
    mut2 <- rnorm(Count[4], mean = (N10_mean + j/2), sd = N10_sd)
    SimString <- c(N10, func, mut1, mut2)
    Genotype <- c(rep("N10", Count[1]), rep("func", Count[2]), rep("mut1", Count[3]), rep("mut2", Count[4]))
    SimDF <- cbind(SimString, Genotype)
    SimDF <- as.data.frame(SimDF, row.names = NULL, optional = FALSE, make.names = TRUE, stringsAsFactors = FALSE)
    model_1 <- aov(SimString ~ Genotype, data = SimDF)
    Tukey <- TukeyHSD(model_1, which="Genotype")
    Table <- as.data.frame(Tukey$Genotype)
    INDEX1 <- grep("N10-func", rownames(Table))
    p1 <- Table[INDEX1, 4]
    Sig1 <- ifelse(p1 <= 0.05, "yes", "no")
    INDEX2 <- grep("mut2-func", rownames(Table))
    p2 <- Table[INDEX2, 4]
    Sig2 <- ifelse(p2 <= 0.05, "yes", "no")
    ROW <- cbind(i, Sig1, Sig2)
    colnames(ROW) <- c("Iteration","N10-func","func-mut2")
    OUT <<- rbind (OUT,ROW)
  }
  TEST <- subset(OUT, `N10-func` == "yes", `func-mut2` == "yes")
  Power <- nrow(TEST)/n
  ROW2 <- cbind(j, Power)
  colnames(ROW2) <- c("Difference", "Power")
  FIN <<- rbind(FIN, ROW2)
  OUT <- data.frame("Iteration"=c(NA), "N10" = c(NA), "trkin_functional"= c(NA), "trkin_mut1"= c(NA), "trkin_mut2" = c(NA))
  OUT <- OUT[c(1),]
  OUT <- OUT[-c(1),]
}

DATA_3 <- subset(DATA, Error_Type == "% Total Meiotic Errors" & Allele == "N10")
N10_mean <- mean(DATA_3$Perc_Error)
N10_sd <- sd(DATA_3$Perc_Error)
DATA_3 <- subset(DATA, Error_Type == "% Total Meiotic Errors")
Count <- summary(as.factor(DATA_3$Allele))

#Looking at power
i=1
n=100
diff= c(1:50)
OUT <- data.frame("Iteration"=c(NA), "N10" = c(NA), "trkin_functional"= c(NA), "trkin_mut1"= c(NA), "trkin_mut2" = c(NA))
OUT <- OUT[-c(1),]
FIN <- data.frame("Difference"=c(NA), "Power" = c(NA))
FIN <- FIN[-c(1),]

for(j in diff) {
  print(j)
  for(i in 1:n) {
    N10 <- rnorm(Count[1], mean = (N10_mean), sd = N10_sd)
    func <- rnorm(Count[2], mean = (N10_mean + j), sd = N10_sd)
    mut1 <- rnorm(Count[3], mean = (N10_mean + j/2), sd = N10_sd)
    mut2 <- rnorm(Count[4], mean = (N10_mean + j/2), sd = N10_sd)
    SimString <- c(N10, func, mut1, mut2)
    Genotype <- c(rep("N10", Count[1]), rep("func", Count[2]), rep("mut1", Count[3]), rep("mut2", Count[4]))
    SimDF <- cbind(SimString, Genotype)
    SimDF <- as.data.frame(SimDF, row.names = NULL, optional = FALSE, make.names = TRUE, stringsAsFactors = FALSE)
    model_1 <- aov(SimString ~ Genotype, data = SimDF)
    Tukey <- TukeyHSD(model_1, which="Genotype")
    Table <- as.data.frame(Tukey$Genotype)
    INDEX1 <- grep("N10-func", rownames(Table))
    p1 <- Table[INDEX1, 4]
    Sig1 <- ifelse(p1 <= 0.05, "yes", "no")
    INDEX2 <- grep("mut2-func", rownames(Table))
    p2 <- Table[INDEX2, 4]
    Sig2 <- ifelse(p2 <= 0.05, "yes", "no")
    ROW <- cbind(i, Sig1, Sig2)
    colnames(ROW) <- c("Iteration","N10-func","func-mut2")
    OUT <<- rbind (OUT,ROW)
  }
  TEST <- subset(OUT, `N10-func` == "yes", `func-mut2` == "yes")
  Power <- nrow(TEST)/n
  ROW2 <- cbind(j, Power)
  colnames(ROW2) <- c("Difference", "Power")
  FIN <<- rbind(FIN, ROW2)
  OUT <- data.frame("Iteration"=c(NA), "N10" = c(NA), "trkin_functional"= c(NA), "trkin_mut1"= c(NA), "trkin_mut2" = c(NA))
  OUT <- OUT[c(1),]
  OUT <- OUT[-c(1),]
}

DATA_3 <- subset(DATA, Error_Type == "% Dyad Micronuclei" & Allele == "N10")
N10_mean <- mean(DATA_3$Perc_Error)
N10_sd <- sd(DATA_3$Perc_Error)
DATA_3 <- subset(DATA, Error_Type == "% Dyad Micronuclei")
Count <- summary(as.factor(DATA_3$Allele))

#Looking at power
i=1
n=100
diff= c(1:50)
OUT <- data.frame("Iteration"=c(NA), "N10" = c(NA), "trkin_functional"= c(NA), "trkin_mut1"= c(NA), "trkin_mut2" = c(NA))
OUT <- OUT[-c(1),]
FIN <- data.frame("Difference"=c(NA), "Power" = c(NA))
FIN <- FIN[-c(1),]

for(j in diff) {
  print(j)
  for(i in 1:n) {
    N10 <- rnorm(Count[1], mean = (N10_mean), sd = N10_sd)
    func <- rnorm(Count[2], mean = (N10_mean + j), sd = N10_sd)
    mut1 <- rnorm(Count[3], mean = (N10_mean + j/2), sd = N10_sd)
    mut2 <- rnorm(Count[4], mean = (N10_mean + j/2), sd = N10_sd)
    SimString <- c(N10, func, mut1, mut2)
    Genotype <- c(rep("N10", Count[1]), rep("func", Count[2]), rep("mut1", Count[3]), rep("mut2", Count[4]))
    SimDF <- cbind(SimString, Genotype)
    SimDF <- as.data.frame(SimDF, row.names = NULL, optional = FALSE, make.names = TRUE, stringsAsFactors = FALSE)
    model_1 <- aov(SimString ~ Genotype, data = SimDF)
    Tukey <- TukeyHSD(model_1, which="Genotype")
    Table <- as.data.frame(Tukey$Genotype)
    INDEX1 <- grep("N10-func", rownames(Table))
    p1 <- Table[INDEX1, 4]
    Sig1 <- ifelse(p1 <= 0.05, "yes", "no")
    INDEX2 <- grep("mut2-func", rownames(Table))
    p2 <- Table[INDEX2, 4]
    Sig2 <- ifelse(p2 <= 0.05, "yes", "no")
    ROW <- cbind(i, Sig1, Sig2)
    colnames(ROW) <- c("Iteration","N10-func","func-mut2")
    OUT <<- rbind (OUT,ROW)
  }
  TEST <- subset(OUT, `N10-func` == "yes", `func-mut2` == "yes")
  Power <- nrow(TEST)/n
  ROW2 <- cbind(j, Power)
  colnames(ROW2) <- c("Difference", "Power")
  FIN <<- rbind(FIN, ROW2)
  OUT <- data.frame("Iteration"=c(NA), "N10" = c(NA), "trkin_functional"= c(NA), "trkin_mut1"= c(NA), "trkin_mut2" = c(NA))
  OUT <- OUT[c(1),]
  OUT <- OUT[-c(1),]
}
