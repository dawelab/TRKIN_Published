
setwd(")

#Load in the new Ab10 annotations
Ab10_New <- read.delim("B73_Ab10_HiFi_v2.trkinreg.K10L2liftoff.newonly.gff3", header=FALSE, comment.char="#")
colnames(Ab10_New) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

#Isolate newly liftedover Ab10 annotations in the trkin region
Ab10_New <- subset(Ab10_New, seqname == "chr10" & start >= 142472000 & end <= 153145000)

#Rename the genes so that they don't overlab with existing annotations and are marked as liftedoff
Ab10_New$attribute <- gsub("=g", "=gK", Ab10_New$attribute)

#Write out the updated file 
write.table(Ab10_New, "B73_Ab10_HiFi_v2.trkinreg.K10L2liftoff.newonly.fix.gff3", sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)

#Load in the new K10L2 annotations
K10L2_New <- read.delim("CI66_K10L2_v1.trkinreg.Ab10liftoff.newonly.gff3", header=FALSE, comment.char="#")
colnames(K10L2_New) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

#Isolate newly liftedover K10L2 annotations in the trkin region
K10L2_New <- subset(K10L2_New, seqname == "K10L2" & start >= 7429619 & end <= 25502283)

#Rename the genes so that they don't overlab with existing annotations and are marked as liftedoff
K10L2_New$attribute <- gsub("=g", "=gA", K10L2_New$attribute)

#Write out the updated file 
write.table(K10L2_New, "CI66_K10L2_v1.trkinreg.Ab10liftoff.newonly.fix.gff3", sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)
