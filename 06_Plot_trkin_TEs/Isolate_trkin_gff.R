K10L2 <- read.delim("/Volumes/Transcend/pasa_db_K10L2.gene_structures_post_PASA_updates.3815155.gff3", header=FALSE, comment.char="#")
colnames(K10L2) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")


Ab10 <- read.delim("/Volumes/Transcend/pasa_db_Ab10.gene_structures_post_PASA_updates.280839.gff3", header=FALSE, comment.char="#")
colnames(Ab10) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")


K10L2_trkinreg <- subset(K10L2, seqname == "K10L2" & start >= 7429619 & end <= 25502283)
K10L2_trkinreg$seqname <- "K10L2"
write.table(K10L2_trkinreg, "/Volumes/Transcend/pasa_db_K10L2.gene_structures_post_PASA_updates.3815155.trkinreg.gff3", sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)

Ab10_trkinreg <- subset(Ab10, seqname == "chr10" & start >= 142472000 & end <= 153145000)
write.table(Ab10_trkinreg, "/Volumes/Transcend/pasa_db_Ab10.gene_structures_post_PASA_updates.280839.trkinreg.gff3", sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)
