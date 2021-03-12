library(tidyr)
library(ggplot2)
library(preprocessCore)
library(ggrepel)
library(gridExtra)
library(edgeR)
library(DESeq2)
library(NMF)
library(GGally)
library(gprofiler2)

# Reading & formatting the data 
df = read.csv("data/processed/rpkm_normalised.csv", sep=',', header=TRUE,row.names=1)
head(df)
#df = df[,-17]
df_log <- log(df+1)
head(df_log)

#Increment function
#inc <- function(x)
#{
#  eval.parent(substitute(x <- x + 2))
#}

#myplots <- list()
#r = 1
#par(mar = rep(2, 4))
#Save to pdf
#pdf("plots/density_plots_quantile_normalised.pdf")
#for (i in seq(1,8))
#{
#r1 = colnames(df_log[r])
#r2 = colnames(df_log[r+1])
#data <- data.frame(r1=df_log[,r],r2 = df_log[,r+1])
#plot(density(data[,1]), bty="n", main=paste("Quantile Normalised - Time Point: ",r1, sep=""), col=r, xlab="log(Gene Counts)") + lines(density(data[,2]), lty=3)
#myplots[[i]] <- p1
#inc(r)
#}
###dev.off()

#Density Plot Function
qc.plots<-function(cts,title){
  cts.tidy = pivot_longer(cts, cols=colnames(cts), names_to='sample', values_to='expression')
  # we now remove 0 counts to create more understandable visualizations
  cts.tidy$expression[cts.tidy$expression == 0] = NA
  cts.tidy$log.expression=log2(cts.tidy$expression)
  
  a<-ggplot(drop_na(cts.tidy), aes(x=log.expression, color=sample)) +
    geom_density(alpha=0.3)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
    theme(plot.title = element_text(size=10)) + labs(title=title)
  a
}


pdf("plots/density_plots_rpkm_normalised.pdf")
qc.plots(df_log, "RPKM Normalised")
dev.off()
