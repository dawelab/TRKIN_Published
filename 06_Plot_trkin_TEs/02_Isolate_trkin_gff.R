setwd("")

#Load both annotation files
K10L2 <- read.delim("B73_Ab10_HiFi_v2.gene.sorted.gff3", header=FALSE, comment.char="#")
colnames(K10L2) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")


Ab10 <- read.delim("CI66_K10L2_v1.gene.sorted.gff3", header=FALSE, comment.char="#")
colnames(Ab10) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

#Isolate the K10L2 trkin annotaions
K10L2_trkinreg <- subset(K10L2, seqname == "K10L2" & start >= 7429619 & end <= 25502283)
K10L2_trkinreg$seqname <- "K10L2"
#Write out the K10L2 trkin annotations
write.table(K10L2_trkinreg, "CI66_K10L2_v1.trkinreg.gff3", sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)

#Isolate the Ab10 trkin annotations
Ab10_trkinreg <- subset(Ab10, seqname == "chr10" & start >= 142472000 & end <= 153145000)
#Write out the Ab1o trkin region
write.table(Ab10_trkinreg, "B73_Ab10_HiFi_v2.trkinreg.gff3", sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)
