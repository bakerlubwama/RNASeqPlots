library(GGally)
# Simulating Data Set - Data Engineering ----------------------------------
counts <- read.csv("data/raw/counts_raw.csv", header=TRUE, row.names=1)
gene_ids <- rownames(counts)
head(gene_ids)
gene_ids <- as.data.frame(gene_ids, header=TRUE)
head(gene_ids)
simulated_data <- gene_ids
head(simulated_data)


#load & merge normalised data
#rpm
rpm <- read.csv("data/processed/rpm_normalised.csv",header=TRUE)
rpm <- rpm[,1:2]
head(rpm)

add_to_sim <- function(sim_data, add, colname)
  {
  sim_data <- merge(sim_data, add, by.x="gene_ids", by.y="genes")
  head(sim_data)
  colnames(sim_data)[ncol(sim_data)] <- (colname)
  return(sim_data)
}

simulated_data <- add_to_sim(sim_data=simulated_data, add=rpm, colname="rpm")
head(simulated_data)

#rpkm
rpkm <- read.csv("data/processed/rpkm_normalised.csv",header=TRUE)
head(rpkm)
rpkm <- rpkm[,1:2]
head(rpkm)

simulated_data <- add_to_sim(sim_data=simulated_data, add=rpkm, colname="rpkm")
head(simulated_data)

#edgeR(tmm)
tmm <- read.csv("data/processed/tmm_normalised.csv", header=TRUE)
head(tmm)
tmm <- tmm[,1:2]
colnames(tmm)[1] <- "genes"
head(tmm)

simulated_data <- add_to_sim(sim_data = simulated_data, add=tmm, colname = "tmm")
head(simulated_data)

#DESeq
deseq <- read.csv("data/processed/deseq_normalised.csv", header=TRUE)
head(deseq)
deseq <- deseq[,1:2]
colnames(deseq)[1] <- "genes"
head(deseq)

simulated_data <- add_to_sim(sim_data = simulated_data, add=deseq, colname = "deseq")
head(simulated_data)

#Upper Quartile
uq <- read.csv("data/processed/upper_quartile_normalised.csv", header=TRUE)
head(uq)
uq <- uq[,1:2]
colnames(uq)[1] <- "genes"
head(uq)

simulated_data <- add_to_sim(sim_data = simulated_data, add=uq, colname = "uq")
head(simulated_data)

#Quantile 
quantile <- read.csv("data/processed/quantile_normalised_counts.csv", header=TRUE)
head(quantile)
quantile <- quantile[,1:2]
colnames(quantile)[1] <- "genes"
head(quantile)

simulated_data <- add_to_sim(sim_data = simulated_data, add=quantile, colname = "quantile")
head(simulated_data)



#Unnormalised
unnorm <- read.csv("data/raw/counts_raw.csv", header=TRUE)
head(unnorm)
unorm <- unnorm[,1:2]
colnames(unnorm)[1] <- "genes"
head(unorm)

simulated_data <- add_to_sim(sim_data = simulated_data, add=unorm, colname = "Unnormalised")
head(simulated_data)
write.csv(simulated_data, "data/processed/simulated_data.csv")

#ggpairs
simulated_data <- read.csv("data/processed/simulated_data.csv", header=TRUE, row.names=1)
simulated_data <- simulated_data[,-1]
head(simulated_data)
pdf('plots/ggpairs_simulated_data.pdf')
ggpairs(simulated_data)
dev.off()
