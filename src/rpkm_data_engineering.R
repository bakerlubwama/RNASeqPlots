# Getting Gene Length -----------------------------------------------------
library(GenomicFeatures)
library(GenomeInfoDb)
library(rtracklayer)
library(tidyverse)

gtf_txdb <- makeTxDbFromGFF("data/raw/Mus_musculus.GRCm39.103.chr.gtf")
gene_list <- genes(gtf_txdb)
gene_list <- as.data.frame(gene_list, header=TRUE)
head(gene_list)

#Loading counts data
counts <- read.csv("data/processed/counts_timepoints_summed.csv", header=TRUE)
head(counts)
counts <- counts[,-18]
head(counts)

#Add gene lengths to count data
counts_with_length <- merge(counts, gene_list, by.x="genes", by.y="gene_id")
head(counts_with_length)
counts_with_length <- as.data.frame(counts_with_length, header=TRUE)
counts_with_length <- counts_with_length[,-c(18,19,20,22)]
head(counts_with_length)

#save data
write.csv(counts_with_length, "data/processed/counts_with_length.csv")

#load normalised rpkm data
rpkm_norm <- read.csv("data/processed/rpkm_normalised.csv", header=TRUE, row.names = 1)
head(rpkm_norm)

#rpkm boxplot
pdf("plots/rpkm_boxplot.pdf")
boxplot(log2(rpkm_norm+1), xlab="Samples and Replicates (to the right)", ylab="Log Abundance",main="RPKM normalised counts")
dev.off()
