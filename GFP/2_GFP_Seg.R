library(readxl)
library(ggplot2)

GFP_Seg <- read_excel("~/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/GFP_Seg.xlsx")

model1 <- aov(GFP_Seg ~ Genotype, data=GFP_Seg)

summary(model1)

TukeyHSD(model1)

SIG <- data.frame(Genotype=c("Ab10-I trkin -", "Ab10-I trkin +", "K10L2 trkin -", "K10L2 trkin +", "N10"), Location=c(90, 90, 60, 60, 60), Sig= c("a", "a", "b", "b", "b"))

pdf("GFP_Seg.pdf", height = 4, width = 6) 
a <- ggplot() +
  geom_jitter(data=GFP_Seg, aes(x=Genotype, y=GFP_Seg, size=Total), alpha=0.5, height=0, color="grey30") +
  geom_boxplot(data=GFP_Seg, aes(x=Genotype, y=GFP_Seg), alpha=0.5, outlier.shape=NA) +
  geom_text(data=SIG, aes(x=Genotype, y=Location, label=Sig), size = 10) +
  ylab("Proportion 4L Mixed Knob") +
  xlab("trkin Genotype") +
  guides(color="none") +
  scale_x_discrete(labels=c("Ab10-I trkin -" = "Ab10\n1 -\n2 -", "Ab10-I trkin +" = "Ab10\n1 +\n2 -", "K10L2 trkin -" = "K10L2\n -", "K10L2 trkin +"="K10L2\n +",  "N10" = "N10")) +
  theme(plot.title = element_text(hjust=0.5, size = 20), plot.subtitle = element_text(hjust=0.5, size = 15), axis.title = element_text(size = 18), axis.text.y = element_text(size = 15), axis.text.x = element_text(size = 8, angle = 0), legend.title = element_text(size = 18), legend.text = element_text(size = 15), legend.position = "right") 
a
dev.off()



