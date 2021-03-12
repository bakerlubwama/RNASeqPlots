#load libraries
library(ggplot2)
library(GGally)

#load data
my.data <- read.csv("data/processed/counts_timepoints_summed2.csv", row.names="genes", header=TRUE)
#delete last column (total gene counts)
my.data <- my.data[,-17]
head(my.data)

#summary stats
pdf("plots/summary_stats_table.pdf")
summary(data)

#Using 1hr and 12hr samples only for ggpairs.
ggpairs_data <- my.data
ggpairs_data <- my.data[,c("X0hr","X0hr.1","X12hr","X12hr.1")]
head(ggpairs_data)
#ggpairs
pdf("plots/ggpairs.pdf")
ggpairs(ggpairs_data)
dev.off()

#boxplot
pdf("plots/boxplots_log.pdf")
boxplot(log2(my.data+1))
dev.off()
