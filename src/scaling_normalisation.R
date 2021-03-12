# Visualising RPM Normalised Data -----------------------------------------
rpm_norm <- read.csv("data/processed/rpm_normalised.csv",header=TRUE,row.names=1)
head(rpm_norm)
pdf("plots/rpm_normalised.pdf")
boxplot(log2(rpm_norm+1), xlab="Samples and Replicates (to the right)", ylab="Log Abundance",main="RPM normalised counts")
dev.off()


library(EDASeq)
library(goseq)
data <- read.csv('data/processed/counts_timepoints_summed.csv',header=TRUE)
head(data)
data <- data[,-18]
head(data)
ensemble_list <- data[,1]
head(ensemble_list)
lengths <- getlength(ensemble_list,'mm9','ensGene')
#21stdec 2009-havgardens-26may2014
#may2014