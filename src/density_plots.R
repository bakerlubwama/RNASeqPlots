# Reading & formatting the data 
df = read.csv("counts_only.csv", sep=',', header=TRUE)
head(df)
df_log <- log(df+1)
head(df_log)

#Increment function
inc <- function(x)
{
  eval.parent(substitute(x <- x + 2))
}

#myplots <- list()
r = 1
#par(mar = rep(2, 4))
#Save to pdf
pdf("density_plots.pdf")
for (i in seq(1,8))
{
r1 = colnames(df_log[r])
r2 = colnames(df_log[r+1])
data <- data.frame(r1=df_log[,r],r2 = df_log[,r+1])
plot(density(data[,1]), bty="n", main=paste("Time Point: ",r1, sep=""), col=r, xlab="log(Gene Counts)") + lines(density(data[,2]), lty=3)
#myplots[[i]] <- p1
inc(r)
}
dev.off()

