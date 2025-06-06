GFP_Seg <- read_excel("GFP_Seg.xlsx")

model1 <- aov(GFP_Seg ~ Genotype, data=GFP_Seg)

summary(model1)

TukeyHSD(model1)


#Filters to only the Ab10 trkin1+ 
GFP_Seg_sub <- subset(GFP_Seg, Genotype == "Ab10-I trkin +")

#I am now getting the mean and standard deviation for both 
Ab10_mean <- mean(GFP_Seg_sub$GFP_Seg)
Ab10_sd <- sd(GFP_Seg_sub$GFP_Seg)

GFP_Seg_sub <- subset(GFP_Seg, Genotype == "K10L2 trkin +")
K10L2_mean <- mean(GFP_Seg_sub$GFP_Seg)
K10L2_sd <- sd(GFP_Seg_sub$GFP_Seg)

GFP_Seg_sub <- subset(GFP_Seg, Genotype == "N10")
N10_mean <- mean(GFP_Seg_sub$GFP_Seg)
N10_sd <- sd(GFP_Seg_sub$GFP_Seg)

#This determines the total number of plants per class I had
Count <- summary(as.factor(GFP_Seg$Genotype))

#Drive need to simulate normal 
i=1
n=500
diff= c(1:20)
OUT <- data.frame("Iteration"=c(NA), "Ab10_mut-Ab10_func" = c(NA), "K10L2_mut-K10L2_func"= c(NA))
OUT <- OUT[-c(1),]
FIN <- data.frame("Difference"=c(NA), "Power" = c(NA))
FIN <- FIN[-c(1),]

for(j in diff) {
  print(j)
  for(i in 1:n) {
    Ab10_mut <- rnorm(Count[1], mean = Ab10_mean-j, sd = Ab10_sd)
    Ab10_func <- rnorm(Count[2], mean = Ab10_mean, sd = Ab10_sd)
    K10L2_mut <- rnorm(Count[3], mean = K10L2_mean-j, sd = K10L2_sd)
    K10L2_func <- rnorm(Count[4], mean = K10L2_mean, sd = K10L2_sd)
    N10 <- rnorm(Count[5], mean = N10_mean, sd = N10_sd)
    Drive <- c(Ab10_mut, Ab10_func, K10L2_mut, K10L2_func, N10)
    Genotype <- c(rep("Ab10_mut", Count[1]), rep("Ab10_func", Count[2]), rep("K10L2_mut", Count[3]), rep("K10L2_func", Count[4]), rep("N10", Count[5]))
    SimDrive <- cbind(Drive, Genotype)
    SimDrive <- as.data.frame(SimDrive, row.names = NULL, optional = FALSE, make.names = TRUE, stringsAsFactors = FALSE)
    model_1 <- aov(Drive~ Genotype, data = SimDrive)
    summary(model_1)
    Tukey <- TukeyHSD(model_1, which="Genotype")
    Table <- as.data.frame(Tukey$Genotype)
    INDEX1 <- grep("Ab10_mut-Ab10_func", rownames(Table))
    p1 <- Table[INDEX1, 4]
    Sig1 <- ifelse(p1 <= 0.05, "yes", "no")
    INDEX2 <- grep("K10L2_mut-K10L2_func", rownames(Table))
    p2 <- Table[INDEX2, 4]
    Sig2 <- ifelse(p2 <= 0.05, "yes", "no")
    ROW <- cbind(i, Sig1, Sig2)
    colnames(ROW) <- c("Iteration","Ab10_mut-Ab10_func","K10L2_mut-K10L2_func")
    OUT <<- rbind (OUT,ROW)
  }
  TEST <- subset(OUT, `Ab10_mut-Ab10_func` == "yes" &  `K10L2_mut-K10L2_func` == "yes")
  Power <- nrow(TEST)/n
  ROW2 <- cbind(j, Power)
  colnames(ROW2) <- c("Difference", "Power")
  FIN <<- rbind(FIN, ROW2)
  OUT <- data.frame("Iteration"=c(NA), "Ab10_mut-Ab10_func" = c(NA), "K10L2_mut-K10L2_func"= c(NA))  
  OUT <- OUT[c(1),]
  OUT <- OUT[-c(1),]
}


FIN$Percent_Change <- FIN$Difference*100
FIN$Percent_Power <- FIN$Power*100
