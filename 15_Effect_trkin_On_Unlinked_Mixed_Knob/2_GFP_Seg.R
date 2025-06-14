library(readxl)
library(ggplot2)

setwd("")

GFP_Seg <- read_excel("GFP_Seg.xlsx")

model1 <- aov(GFP_Seg ~ Genotype, data=GFP_Seg)

summary(model1)

TukeyHSD(model1)

SIG <- data.frame(Genotype=c("Ab10-I trkin -", "Ab10-I trkin +", "K10L2 trkin -", "K10L2 trkin +", "N10"), Location=c(90, 90, 60, 60, 60), Sig= c("a", "a", "b", "b", "b"))

pdf("GFP_Seg.pdf", height = 4, width = 4) 
a <- ggplot() +
  geom_jitter(data=GFP_Seg, aes(x=Genotype, y=GFP_Seg, size=Total), alpha=0.5, height=0, color="grey30") +
  geom_boxplot(data=GFP_Seg, aes(x=Genotype, y=GFP_Seg), alpha=0.5, outlier.shape=NA) +
  geom_text(data=SIG, aes(x=Genotype, y=Location, label=Sig), size = 3) +
  ylab("Proportion 4L Mixed Knob") +
  xlab("trkin Genotype") +
  guides(color="none") +
  scale_x_discrete(labels=c("Ab10-I trkin -" = "Ab10\n1(-)\n2(-)", "Ab10-I trkin +" = "Ab10\n1(+)\n2(-)", "K10L2 trkin -" = "K10L2\n (-)", "K10L2 trkin +"="K10L2\n (+)",  "N10" = "N10")) +
  theme(plot.title = element_text(hjust=0.5, size = 10), plot.subtitle = element_text(hjust=0.5, size = 10), axis.title = element_text(size = 10), axis.text.y = element_text(size = 10), axis.text.x = element_text(size = 10, angle = 0), legend.title = element_text(size = 10), legend.text = element_text(size = 10), legend.position = "right") 
a
dev.off()



