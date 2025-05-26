Ab10_New <- read.delim("/Volumes/Transcend/Ab10_HiFi_v2_corrected.trkinreg.K10L2liftoff.newonly.gff3", header=FALSE, comment.char="#")
colnames(Ab10_New) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

Ab10_New <- subset(Ab10_New, seqname == "chr10" & start >= 142472000 & end <= 153145000)

#this renames the genes so that they don't overlab with existing annotations and I don't need to renumber all the genes
Ab10_New$attribute <- gsub("=g", "=gK", Ab10_New$attribute)

write.table(Ab10_New, "/Volumes/Transcend/Ab10_HiFi_v2_corrected.trkinreg.K10L2liftoff.newonly.fix.gff3", sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)


K10L2_New <- read.delim("/Volumes/Transcend/CI66_K10L2_v1.trkinreg.Ab10liftoff.newonly.gff3", header=FALSE, comment.char="#")
colnames(K10L2_New) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

K10L2_New <- subset(K10L2_New, seqname == "K10L2" & start >= 7429619 & end <= 25502283)

#this renames the genes so that they don't overlab with existing annotations and I don't need to renumber all the genes
K10L2_New$attribute <- gsub("=g", "=gA", K10L2_New$attribute)

write.table(K10L2_New, "/Volumes/Transcend/CI66_K10L2_v1.trkinreg.Ab10liftoff.newonly.fix.gff3", sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)
