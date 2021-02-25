# Read data 
d = read.csv("data/raw/counts_raw.csv", header=TRUE)
d = data.frame(d, row.names=1)
head(d)
#load library
library(optimbase)
library(IntClust)
#Define Jaccard Similarity Index Function
jaccard_index <- function(i, j, k)
{
  si = data.frame(d[1:k,i], row.names=rownames(d)[1:k])
  sj = data.frame(d[1:k,j], row.names=rownames(d)[1:k])
  colnames(si) <- c("counts")
  colnames(sj) <- c("counts")
  sorted_i = si[order(si$counts, decreasing=TRUE), ,drop=FALSE]
  sorted_j = sj[order(si$counts, decreasing=TRUE), ,drop=FALSE]
  si_genes = rownames(sorted_i)[1:k]
  sj_genes = rownames(sorted_j)[1:k]
  intersect = intersect(si_genes, sj_genes)
  union = union(si_genes, sj_genes)
  JSI = length(intersect)/length(union)
  return(JSI)
}

# Define increment function
inc <- function(x)
{
  eval.parent(substitute(x <- x + 100))
}
length(d)
#Initiate K
k = 100
#Save pdf
pdf("data/jaccard_heatmaps.pdf")
for (i in seq(1,10))
{
  jsi_matrix <- ones(nx=length(d), ny=length(d))
  for (i in seq(1, 15))
  {
    for (j in seq(i+1, 16))
    {
      JSI <- jaccard_index(i, j, k)
      jsi_matrix[i,j] <- JSI
      jsi_matrix[j,i] <- JSI
    }
  }
  jsi_matrix
  heatmap = heatmap(jsi_matrix, labRow= c("0hr","0hr.1","1hr","0hr.1","6hr","6hr.1","12hr","12hr.1","24hr","24hr.1","36hr","36hr.1","48hr","48hr.1","72hr","72hr.1"),
                    labCol= c("0hr","0hr.1","1hr","0hr.1","6hr","6hr.1","12hr","12hr.1","24hr","24hr.1","36hr","36hr.1","48hr","48hr.1","72hr","72hr.1"),
                    main=paste("Top", k, "Genes"), Rowv=NA, Colv=NA)
  inc(k)
}
dev.off()