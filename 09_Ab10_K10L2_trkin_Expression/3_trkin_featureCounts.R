library(Rsubread)

setwd("/scratch/mjb51923/TRKIN_CRISPR/out/RNA_Aln/")

files <- c("Ab10IMMR_1_v_B73Ab10.s.bam", "Ab10IMMR_2_v_B73Ab10.s.bam", "Ab10IMMR_3_v_B73Ab10.s.bam", "K10L2_1_v_B73Ab10.s.bam", "K10L2_2_v_B73Ab10.s.bam", "K10L2_3_v_B73Ab10.s.bam", "K10L2_4_v_B73Ab10.s.bam", "K10L2_5_v_B73Ab10.s.bam", "N10_1_Ab10Sib_v_B73Ab10.s.bam", "N10_1_v_B73Ab10.s.bam", "N10_2_Ab10Sib_v_B73Ab10.s.bam", "N10_2_v_B73Ab10.s.bam", "N10_3_Ab10Sib_v_B73Ab10.s.bam", "N10_3_v_B73Ab10.s.bam")

featureCounts(files, annot.ext="/scratch/mjb51923/annotations/Zm-B73_AB10-REFERENCE-NAM-1.0_Zm00043a.1.gtf", isGTFAnnotationFile=TRUE, isPairedEnd=TRUE, minMQS=20, countMultiMappingReads=FALSE, largestOverlap=TRUE, nthreads=12)
