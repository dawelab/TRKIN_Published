library(Rsubread)

setwd("")

files <- c("B73Ab10_V11_middle_MN02041.s.bam", "B73Ab10_V18_tassel_MN02061.s.bam", "B73Ab10_R1_anther_MN02082.s.bam", "B73Ab10_16DAP_endosperm_MN02092.s.bam", "B73Ab10_V11_tip_MN02052.s.bam", "B73Ab10_16DAP_embryo_MN02101.s.bam", "B73Ab10_8DAS_root_MN02012.s.bam", "B73Ab10_16DAP_endosperm_MN02091.s.bam", "B73Ab10_16DAP_embryo_MN02102.s.bam", "B73Ab10_R1_anther_MN02081.s.bam", "B73Ab10_V11_tip_MN02051.s.bam", "B73Ab10_V18_ear_MN02072.s.bam", "B73Ab10_8DAS_root_MN02011.s.bam", "B73Ab10_V11_base_MN02032.s.bam", "B73Ab10_V18_ear_MN02071.s.bam", "B73Ab10_8DAS_shoot_MN02022.s.bam", "B73Ab10_V11_base_MN02031.s.bam", "B73Ab10_8DAS_shoot_MN02021.s.bam", "B73Ab10_V11_middle_MN02042.s.bam", "B73Ab10_V18_tassel_MN02062.s.bam")

COUNTS <- featureCounts(files, annot.ext="/path/to/Liftoff/B73_Ab10_HiFi_v2.gene.gtf", isGTFAnnotationFile=TRUE, isPairedEnd=TRUE, minMQS=20, countMultiMappingReads=FALSE, largestOverlap=TRUE, nthreads=12)

write.table(COUNTS$counts, file = "FeatureCounts_counts.txt", row.names = TRUE)
write.table(COUNTS$annotation, file = "FeatureCounts_annotation.txt", row.names = FALSE)
write.table(COUNTS$target, file = "FeatureCounts_target.txt", row.names = FALSE)
write.table(COUNTS$stat, file = "FeatureCounts_stat.txt", row.names = FALSE)
