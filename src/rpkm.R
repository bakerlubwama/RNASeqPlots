# Getting Gene Length -----------------------------------------------------
library(GenomicFeatures)
library(GenomeInfoDb)
library(rtracklayer)
library(tidyverse)

gtf_txdb <- makeTxDbFromGFF("data/raw/Mus_musculus.GRCm39.103.chr.gtf")
gene_list <- genes(gtf_txdb)
gene_list <- as.data.frame(gene_list. header=TRUE)
gene_list <- as_tibble(gene_list)
head(gene_list)

#Loading counts data
counts <- read.csv("data/processed/counts_timepoints_summed.csv", header=TRUE, row.names=1)
gene_ids = rownames(counts)
head(gene_ids)

#Subsetting Gene list to get desired lengths
filter(gene_list, gene_id==gene_ids)
