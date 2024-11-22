library(tidyverse)
library(reshape2)
library(splitstackshape)
library(karyoploteR)
library(readxl)
library(ggplot2)
library(pafr)
library(Rsamtools)

#Set the working directory
setwd("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/Paper/TRKIN_Paper")

#Load in the GFF files and the EnTap files for K10L2 and Ab10
K10L2_GFF <- read.delim("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/OrthoFinder/CI66_K10L2.genes.edit.gff3", header = FALSE)
colnames(K10L2_GFF) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
K10L2_entap  <- read.delim("~/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/OrthoFinder/K10L2_entap_results.tsv")

Ab10_GFF <- read.delim("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/OrthoFinder/HiFiAb10.genes.edit.gff3", header = FALSE)
#This drops an unnecessary column at the top
Ab10_GFF <- Ab10_GFF[-c(1),]
colnames(Ab10_GFF) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
Ab10_entap <- read.delim("~/University_of_Georgia/Dawe_Lab_Documents/Trkin_CRISPR/R Sessions/OrthoFinder/Ab10_entap_results.tsv")

#This reformats the GFF file to be compatible with the entap file
K10L2_GFF_GENE <- subset(K10L2_GFF, feature == "gene")
K10L2_GFF_GENE_temp1 <- separate(K10L2_GFF_GENE, col=attribute, into= c("ID", "biotype", "logic_name"), sep= ';')
K10L2_GFF_GENE_temp1$biotype <- gsub("Name=", "", K10L2_GFF_GENE_temp1$biotype)
K10L2_GFF_GENE <- K10L2_GFF_GENE_temp1[,c("seqname", "start", "end", "biotype")]

#This reformats the GFF file to be compatible with the entap file
Ab10_GFF_GENE <- subset(Ab10_GFF, feature == "gene")
Ab10_GFF_GENE_temp1 <- separate(Ab10_GFF_GENE, col=attribute, into= c("ID", "biotype", "logic_name"), sep= ';')
Ab10_GFF_GENE_temp1$biotype <- gsub("Name=", "", Ab10_GFF_GENE_temp1$biotype)
Ab10_GFF_GENE <- Ab10_GFF_GENE_temp1[,c("seqname", "start", "end", "biotype")]

#This merges the GFF files and the entap output
Ab10 <- merge(Ab10_GFF_GENE, Ab10_entap, by.x="biotype", by.y="Query.Sequence")
K10L2 <- merge(K10L2_GFF_GENE, K10L2_entap, by.x="biotype", by.y="Query.Sequence")

#This isolates the SGE specific regions
Ab10_SPEC <- subset(Ab10, seqname == "chr10" & ((start >= 142472000 & end <= 153145000) | (start >= 167721000)))
colnames(Ab10_SPEC) <- c("gene", "chr", "start", "end", "Frame", "Subject.Sequence", "Percent.Identical", "Alignment.Length", "Mismatches", "Gap.Openings", "Query.Start", "Query.End", "Subject.Start", "Subject.End", "E.Value", "Coverage", "Description", "Species", "Taxonomic.Lineage", "Origin.Database", "Contaminant", "Informative", "UniProt.Database.Cross.Reference", "UniProt.Additional.Information", "UniProt.KEGG.Terms", "UniProt.GO.Biological", "UniProt.GO.Cellular", "UniProt.GO.Molecular", "EggNOG.Seed.Ortholog", "EggNOG.Seed.E.Value", "EggNOG.Seed.Score", "EggNOG.Predicted.Gene", "EggNOG.Tax.Scope", "EggNOG.Tax.Scope.Max", "EggNOG.Member.OGs", "EggNOG.Description", "EggNOG.KEGG.Terms", "EggNOG.GO.Biological", "EggNOG.GO.Cellular", "EggNOG.GO.Molecular", "EggNOG.Protein.Domains", "X")

K10L2_SPEC <- subset(K10L2, seqname == "K10L2" & start >= 7429619 & end <= 25502283)
colnames(K10L2_SPEC) <- c("gene", "chr", "start", "end", "Frame", "Subject.Sequence", "Percent.Identical", "Alignment.Length", "Mismatches", "Gap.Openings", "Query.Start", "Query.End", "Subject.Start", "Subject.End", "E.Value", "Coverage", "Description", "Species", "Taxonomic.Lineage", "Origin.Database", "Contaminant", "Informative", "UniProt.Database.Cross.Reference", "UniProt.Additional.Information", "UniProt.KEGG.Terms", "UniProt.GO.Biological", "UniProt.GO.Cellular", "UniProt.GO.Molecular", "EggNOG.Seed.Ortholog", "EggNOG.Seed.E.Value", "EggNOG.Seed.Score", "EggNOG.Predicted.Gene", "EggNOG.Tax.Scope", "EggNOG.Tax.Scope.Max", "EggNOG.Member.OGs", "EggNOG.Description", "EggNOG.KEGG.Terms", "EggNOG.GO.Biological", "EggNOG.GO.Cellular", "EggNOG.GO.Molecular", "EggNOG.Protein.Domains", "X")


#This writes out the file
write.csv(Ab10_SPEC, "Ab10SpecificGenes_FunctionalAnnotation.csv", row.names = FALSE, quote = FALSE)
write.csv(K10L2_SPEC, "K10L2SpecificGenes_FunctionalAnnotation.csv", row.names = FALSE, quote = FALSE)

