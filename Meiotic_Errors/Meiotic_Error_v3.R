library(ggplot2)
library("readxl")
library("reshape2")
library(gvlma)

DATA <- read_excel("~/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/Micronuclei_Count_format2.xlsx")

DATA <- subset(DATA, Error_Type != "Minicell/Tetrad Cell")

colors <- c("blue4", "darkseagreen4", "darkred", "darkgoldenrod", "darkorchid")

# Doing the math
DATA_1 <- subset(DATA, Error_Type == "% Tetrad Micronuclei")
model_1 <- lm(Perc_Error ~ Allele, data = DATA_1)
gvlma(model_1)
model_1 <- aov(Perc_Error ~ Allele, data = DATA_1)
summary(model_1)

DATA_2 <- subset(DATA, Error_Type == "% Tetrad Microcyte")
model_2 <- lm(Perc_Error ~ Allele, data = DATA_2)
gvlma(model_2)
model_2 <- aov(Perc_Error ~ Allele, data = DATA_2)
summary(model_2)
#This doesn't meet all the assumptions, but it's close

DATA_3 <- subset(DATA, Error_Type == "% Total Meiotic Errors")
model_3 <- lm(Perc_Error ~ Allele, data = DATA_3)
gvlma(model_3)
model_3 <- aov(Perc_Error ~ Allele, data = DATA_3)
summary(model_3)

DATA_4 <- subset(DATA, Error_Type == "% Dyad Micronuclei")
model_4 <- lm(Perc_Error ~ Allele, data = DATA_4)
gvlma(model_4)
model_4 <- aov(Perc_Error ~ Allele, data = DATA_4)
summary(model_4)
#This doesn't meet all the assumptions, but it's close

pdf("TRKIN_Mut_Meiotic_Errors_FacetWrap.pdf", height = 6, width = 6)
p <- ggplot(data = DATA, aes(x=Allele)) +
  geom_boxplot(aes(y=Perc_Error), alpha=0.2, outlier.shape = NA) +
  geom_jitter(aes(y=Perc_Error, size=Total_Cells), alpha=0.5, color = "navy", height=0 ) +
  geom_vline(xintercept = 1.5, linetype="dotted", color="darkgrey") +
  geom_vline(xintercept = 2.5, linetype="dotted", color="darkgrey") +
  geom_vline(xintercept = 3.5, linetype="dotted", color="darkgrey") +
  scale_x_discrete(labels=c("None (N10)", "Functional","Mutant 1","Mutant 2")) +
  scale_color_manual(values = colors) +
  guides(color = FALSE) +
  scale_x_discrete(labels=c("trkin_mut2" = "Ab10\n1 -\n2 -", "trkin_mut1" = "Ab10\n1 -\n2 +", "trkin_functional" = "Ab10\n1 +\n2 -", "N10" = "N10")) +
  labs(size="Total Cell Number", x = "trkin Allele", y = "Percent Meiotic Error") +
  theme(plot.title = element_text(hjust=0.5, size = 20), axis.title = element_text(size = 18), axis.text = element_text(size = 15), legend.title = element_text(size = 18), legend.text = element_text(size = 15), plot.subtitle = element_text(hjust=0.5, size = 18), legend.position="bottom", strip.text.x = element_text(size = 18)) +
  facet_wrap(~Error_Type)
p  
dev.off()


pdf("TRKIN_Mut_Meiotic_Errors_FacetWrap_zoom.pdf", height = 6, width = 6)
p <- ggplot(data = DATA, aes(x=Allele)) +
  geom_boxplot(aes(y=Perc_Error), alpha=0.2, outlier.shape = NA) +
  geom_jitter(aes(y=Perc_Error, size=Total_Cells), alpha=0.5, color = "navy", height=0 ) +
  geom_vline(xintercept = 1.5, linetype="dotted", color="darkgrey") +
  geom_vline(xintercept = 2.5, linetype="dotted", color="darkgrey") +
  geom_vline(xintercept = 3.5, linetype="dotted", color="darkgrey") +
  scale_x_discrete(labels=c("None (N10)", "Functional","Mutant 1","Mutant 2")) +
  scale_color_manual(values = colors) +
  guides(color = FALSE) +
  scale_x_discrete(labels=c("trkin_mut2" = "Ab10\n1 -\n2 -", "trkin_mut1" = "Ab10\n1 -\n2 +", "trkin_functional" = "Ab10\n1 +\n2 -", "N10" = "N10")) +
  ylim(c(0,10)) +
  labs(size="Total Cell Number", x = "trkin Allele", y = "Percent Meiotic Error", color="Meiotic Error Type") +
  theme(plot.title = element_text(hjust=0.5, size = 20), axis.title = element_text(size = 18), axis.text = element_text(size = 15), legend.title = element_text(size = 18), legend.text = element_text(size = 15), plot.subtitle = element_text(hjust=0.5, size = 18), legend.position="bottom", strip.text.x = element_text(size = 18)) +
  facet_wrap(~Error_Type)
p  
dev.off()
