#Initiating library
library("ggplot2")

#Uploading data 
#X = read.csv("counts_raw.csv")





# Loading data ------------------------------------------------------------

#RPM Data
rpm_norm <- read.csv("data/processed/rpm_normalised.csv",header=TRUE,row.names=1)
rpm_norm_2 <- rpm_norm[,c("hr0","hr0.1","hr12","hr12.1")]
head(rpm_norm_2)

#TMM Data
tmm_norm <- read.csv("data/processed/tmm_normalised.csv",header=TRUE,row.names=1)
head(tmm_norm)
tmm_norm <- tmm_norm[,c("hr0","hr0.1","hr12","hr12.1")]
head(tmm_norm)

#DESeq Data
deseq_norm <- read.csv("data/processed/deseq_normalised.csv",header=TRUE,row.names=1)
head(deseq_norm)
deseq_norm <- deseq_norm[,c("X0hr","X0hr.1","X12hr","X12hr.1")]
head(deseq_norm)

#Quantile Data 
quantile_norm <- read.csv("data/processed/quantile_normalised_counts.csv", header=TRUE, row.names=1)
head(quantile_norm)
quantile_norm <- quantile_norm[,c("hr0","hr0.1","hr12","hr12.1")]
head(quantile_norm)

#Upper Quartile Data
uq_norm <- read.csv("data/processed/upper_quartile_normalised.csv", header=TRUE, row.names=1)
head(uq_norm)
uq_norm <- uq_norm[,c("hr0","hr0.1","hr12","hr12.1")]
head(uq_norm)

#RPKM data
rpkm_norm <- read.csv("data/processed/rpkm_normalised.csv", header=TRUE, row.names=1)
head(rpkm_norm)
rpkm_norm <- rpkm_norm[,c("hr0","hr0.1","hr12","hr12.1")]
head(rpkm_norm)
# Functions and variables -------------------------------------------------

#incrementation function
inc <- function(x)
{
  eval.parent(substitute(x <- x + 12))
}

#MA Plot function 
ma_plot <- function(X,start_col,end_col,title,output_name)
{
  t=0
  myplots <- list()
  for (i in seq(start_col, end_col, 2))
  {
    A = log2((X[,i] + X[,i+1]+2)/2)
    M = log2((X[,i] +1)/(X[,i+1]+1))
    #a.categ = floor(A/2)
    #a.categ = as.factor(a.categ)
    madata = data.frame(M,A)
    p1 <- eval(substitute(ggplot(data = madata, aes(x = A, y = M)) + geom_point() + geom_hline(yintercept=0, color = "red") + geom_hline(yintercept=0.5, color = "blue") + geom_hline(yintercept=-0.5, color = "blue") + geom_hline(yintercept=1, color = "green") + geom_hline(yintercept=-1, color = "green") + labs(title=paste(title,t,"hr"))))
    myplots[[i]] <- p1
    inc(t)
  }
  pdf(paste("plots/",output_name))
  myplots
} 


# Generating plots --------------------------------------------------------

#MA plot

ma_plot(X=rpkm_norm,start_col=1,end_col=4,title="RPKM Normalised (MA Plot)",output_name ="rpkm_ma.pdf")
dev.off()



