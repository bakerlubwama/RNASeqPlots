#--------------------Scaling Methods--------------------
#TMM normalisation (Using edgeR) -----------------------------------------
library(edgeR)
#load data
x <- read.csv("data/processed/counts_timepoints_summed.csv", row.names="genes")
#delete last column (total gene counts)
x <- x[,-17]
head(x)
#comparing groups
group <- factor(c("0hr", "0hr", "1hr", "1hr", "6hr", "6hr", "12hr", "12hr",
                   "24hr", "24hr", "36hr", "36hr", "48hr", "48hr", "72hr", "72hr"))
y <- DGEList(counts=x, group=group)
#normalise for library size by calculating scaling factor using TMM (default method)
y <- calcNormFactors
#normalistion factors for each library 
y$samples
# count per million read (normalised count)
norm_counts <- cpm(y)
head(norm_counts)


# DESeq2 normalisation  ---------------------------------------------------
library(DESeq2)
#load data
x1 <- read.csv("data/processed/counts_timepoints_summed2.csv", row.names="genes")
cond <- read.csv("data/processed/metadata_no_ids.csv", row.names="Timepoint")
x1
cond$Replicate <- factor(cond$Replicate)
rownames(cond)
# keep only required columsn present in the sample information table
x1 <- x1[, rownames(cond)]
head(x1)
#get dds
dds <- DESeqDataSetFromMatrix(countData = x1, colData = cond, design = ~ Replicate)
dds <- estimateSizeFactors(dds)
#DESeq2 normalisation counts
y2 = counts(dds, normalized = TRUE)
head(y2)

#get size factors
sizeFactors(dds)




