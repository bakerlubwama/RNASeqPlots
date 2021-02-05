#Initiating library
library("ggplot2")

#Uploading data 
X = read.csv("counts_raw.csv")
myplots <- list()

#Generating boxplots
for (i in seq(2, 17, 2))
{
  A = log2((X[,i] + X[,i+1]+2)/2)
  M = log2((X[,i] +1)/(X[,i+1]+1))
  a.categ = floor(A/2)
  a.categ = as.factor(a.categ)
  madata = data.frame(M,A)
  p1 <- eval(substitute(ggplot(data = madata, aes(x = a.categ, y = M)) + geom_boxplot() + geom_hline(yintercept=0, color = "red") + geom_hline(yintercept=0.5, color = "blue") + geom_hline(yintercept=-0.5, color = "blue") + geom_hline(yintercept=1, color = "green") + geom_hline(yintercept=-1, color = "green")))
  myplots[[i]] <- p1
}
pdf("boxplot_floor_over2.pdf")
myplots
dev.off()