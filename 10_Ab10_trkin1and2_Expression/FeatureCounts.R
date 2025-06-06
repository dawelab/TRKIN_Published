library(Rsubread)

setwd("")

files <- c("Ab10IMMR_1_v_B73Ab10.s.bam", "Ab10IMMR_2_v_B73Ab10.s.bam", "Ab10IMMR_3_v_B73Ab10.s.bam", "K10L2_1_v_B73Ab10.s.bam", "K10L2_2_v_B73Ab10.s.bam", "K10L2_3_v_B73Ab10.s.bam", "K10L2_4_v_B73Ab10.s.bam", "K10L2_5_v_B73Ab10.s.bam", "N10_1_Ab10Sib_v_B73Ab10.s.bam", "N10_1_v_B73Ab10.s.bam", "N10_2_Ab10Sib_v_B73Ab10.s.bam", "N10_2_v_B73Ab10.s.bam", "N10_3_Ab10Sib_v_B73Ab10.s.bam", "N10_3_v_B73Ab10.s.bam")

COUNTS <- featureCounts(files, annot.ext="/path/to/Liftoff/B73_Ab10_HiFi_v2.gene.gtf", isGTFAnnotationFile=TRUE, isPairedEnd=TRUE, minMQS=20, countMultiMappingReads=FALSE, largestOverlap=TRUE, nthreads=12)

write.table(COUNTS$counts, file = "FeatureCounts_counts.txt", row.names = TRUE)
write.table(COUNTS$annotation, file = "FeatureCounts_annotation.txt", row.names = FALSE)
write.table(COUNTS$target, file = "FeatureCounts_target.txt", row.names = FALSE)
write.table(COUNTS$stat, file = "FeatureCounts_stat.txt", row.names = FALSE)
