library(preprocessCore)
library(stats)
library(NOISeq)
data <- read.csv("data/processed/counts_timepoints_summed.csv", header=TRUE, row.names=1)
data <- data[,-17]
data_matrix <- as.matrix(data)
head(data_matrix)

# Quantile Normalisation --------------------------------------------------
data_norm <- normalize.quantiles(data_matrix)
head(data_norm)
dimnames(data_norm) <- list(rownames(data), colnames(data))
head(data_norm)
write.csv(data_norm,"data/processed/quantile_normalised_counts.csv")
pdf("plots/quantile_normalised.pdf")
boxplot(log2(data_norm+1), xlab="Samples and Replicates (to the right)", ylab="Log Abundance",main="Quantile normalised counts")
dev.off()


# Upper Quantile Normalisation --------------------------------------------
upper_quartile_normalised <- uqua(data_matrix)
head(upper_quartile_normalised)
write.csv(upper_quartile_normalised,"data/processed/upper_quartile_normalised.csv")

pdf("plots/upper_quartile_normalised_boxplot.pdf")
boxplot(log2(upper_quartile_normalised+1), xlab="Samples and Replicates (to the right)", ylab="Log Abundance",main="Upper quartile normalised counts")
dev.off()

#PCAs for Normalised Data 
# Initiate libraries
library(ggplot2)
library(ggpubr)

#"plots/pcas_quantile_normalised.pdf"
#Generating the PCA plots
#Generating the PCA plots
inc <- function(x)
{
  eval.parent(substitute(x <- x + 100))
}
myplots <- list()
r = 100

par(mar = rep(2, 4))

for (i in seq(1,10))
{
  pca <- prcomp(t(d[1:r,2:17]), scale=TRUE)
  plot(pca$x[,1], pca$x[,2])
  pca.data <- data.frame(Sample=rownames(pca$x), 		X=pca$x[,1], Y=pca$x[,2])
  pca.data
  
  # Variability of each principal component: pca.var
  pca.var <-pca$sdev^2
  # Variance explained by each principal component: pve
  pve <- pca.var / sum(pca.var)
  
  p1 <- eval(substitute(ggplot(data=pca.data, aes(x=X, y=Y, label=Sample)) + geom_text(aes(colour = factor(Sample, 
                                                                                                           levels=c("hr0","hr0.1","hr1","hr1.1","hr6","hr6.1","hr12","hr12.1","hr24","hr24.1","hr36","hr36.1","hr48","hr48.1","hr72","hr72.1"),
                                                                                                           labels=c("Replicate 1","Replicate 2","Replicate 1","Replicate 2","Replicate 1","Replicate 2","Replicate 1","Replicate 2","Replicate 1","Replicate 2","Replicate 1","Replicate 2","Replicate 1","Replicate 2","Replicate 1","Replicate 2")
  ))) + xlab(paste("PC1 (", round(pve[1],2), ") - Top " ,r, " genes", sep="")) + ylab(paste("PC2 (", round(pve[2],2), ") - Top ",r, " genes", sep="")) + theme_bw() + theme(legend.position="bottom") + stat_conf_ellipse(aes(colour = factor(Sample, 
                                                                                                                                                                                                                                            levels=c("hr0","hr0.1","hr1","hr1.1","hr6","hr6.1","hr12","hr12.1","hr24","hr24.1","hr36","hr36.1","hr48","hr48.1","hr72","hr72.1"),
                                                                                                                                                                                                                                            labels=c("Replicate 1","Replicate 2","Replicate 1","Replicate 2","Replicate 1","Replicate 2","Replicate 1","Replicate 2","Replicate 1","Replicate 2","Replicate 1","Replicate 2","Replicate 1","Replicate 2","Replicate 1","Replicate 2")
  )))))
  myplots[[i]] <- p1
  inc(r)
}
#Save to a pdf

pdf("plots/pcas_quantile_normalised.pdf")

myplots

dev.off()
