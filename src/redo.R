#load libraries
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

#load data
unnormalised <- read.csv("data/processed2/unnormalised.csv",header=TRUE, row.names=1)
head(unnormalised)

rpm <- read.csv("data/processed2/rpm.csv",header=TRUE, row.names=1)
head(rpm)
rpm = rpm[,-17]
head(rpm)

rpkm <- read.csv("data/processed2/rpkm.csv",header=TRUE, row.names=1)
head(rpkm)
rpkm = rpkm[,-17]
head(rpkm)

tmm <- read.csv("data/processed2/tmm.csv",header=TRUE, row.names=1)
head(tmm)
tmm = tmm[,-17]
head(tmm)

deseq <- read.csv("data/processed2/deseq.csv",header=TRUE, row.names=1)
head(deseq)
deseq = deseq[,-17]

uq <- read.csv("data/processed2/upper_quartile.csv",header=TRUE, row.names=1)
head(uq)
uq = uq[,-17]
head(uq)

quantile <- read.csv("data/processed2/quantile.csv",header=TRUE, row.names=1)
head(quantile)
meta_data <- read.csv("data/processed2/meta.csv",header=TRUE)

#Filtering
unnormalised.no.zeros = unnormalised[rowSums(unnormalised) > 0,]
expression.threshold <- 20
unnormalised.filtered = unnormalised.no.zeros
unnormalised.filtered[unnormalised.filtered<expression.threshold]<-expression.threshold
unnormalised.filtered = unnormalised.filtered[rowSums(unnormalised.filtered)>expression.threshold*ncol(unnormalised.filtered),]
write.csv(unnormalised.filtered, "data/processed2/unnormalised_filtered.csv")

#Load filtered data
rpm.filtered = read.csv("data/processed2/rpm_filtered.csv", header=TRUE, row.names=1)

rpkm.filtered = read.csv("data/processed2/rpkm_filtered.csv", header=TRUE, row.names=1)

tmm.filtered = read.csv("data/processed2/tmm_filtered.csv", header=TRUE, row.names=1)

deseq.filtered = read.csv("data/processed2/deseq_filtered.csv", header=TRUE, row.names=1)

uq.filtered = read.csv("data/processed2/uq_filtered.csv", header=TRUE, row.names=1)

quantile.filtered = read.csv("data/processed2/quantile_filtered.csv", header=TRUE, row.names=1)

# Quality Control ---------------------------------------------------------

#Visualising distribution of abundance per sample 
qc.plots<-function(cts,title){
  cts.tidy = pivot_longer(cts, cols=colnames(cts), names_to='sample', values_to='expression')
  # we now remove 0 counts to create more understandable visualizations
  cts.tidy$expression[cts.tidy$expression == 0] = NA
  cts.tidy$log.expression=log2(cts.tidy$expression)
  a<-ggplot(drop_na(cts.tidy), aes(x=log.expression, color=sample)) +
    geom_density(alpha=0.3)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    theme(plot.title = element_text(size=10))
  b<-ggplot(drop_na(cts.tidy), aes(x=sample, y=log.expression,color=sample)) +
    geom_violin(alpha=0.3) +
    theme(legend.position = "none")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    theme(plot.title = element_text(size=10))+
    theme(axis.text.x=element_blank())
  c<-ggplot(drop_na(cts.tidy), aes(x=sample, y=log.expression,color=sample)) +
    geom_boxplot(alpha=0.3) +
    theme(legend.position = "none")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    theme(plot.title = element_text(size=10))+
    theme(axis.text.x=element_blank())
  grid.arrange(a,arrangeGrob(b,c),nrow=1,top=title,widths=2:1)
}

#unfiltered qc 
pdf("plots/QC/unfiltered_unnormalised.pdf")
qc.plots(unnormalised, 'Unfiltered and unnormalised')
dev.off()

pdf("plots/QC/unfiltered_rpm.pdf")
qc.plots(rpm, 'Unfiltered and RPM normalised')
dev.off()

pdf("plots/QC/unfiltered_rpkm.pdf")
qc.plots(rpkm, 'Unfiltered and RPKM normalised')
dev.off()

pdf("plots/QC/unfiltered_tmm.pdf")
qc.plots(tmm, 'Unfiltered and TMM (edgeR) normalised')
dev.off()

pdf("plots/QC/unfiltered_deseq.pdf")
qc.plots(deseq, 'Unfiltered and DESeq normalised')
dev.off()

pdf("plots/QC/unfiltered_uq.pdf")
qc.plots(uq, 'Unfiltered and Upper Quartile normalised')
dev.off()

pdf("plots/QC/unfiltered_quantile.pdf")
qc.plots(quantile, 'Unfiltered and Quantile normalised')
dev.off()

#filtered qc 
df("plots/QC/filtered_unnormalised.pdf")
qc.plots(unnormalised.filtered, 'filtered and unnormalised')
dev.off()

pdf("plots/QC/filtered_rpm.pdf")
qc.plots(rpm.filtered, 'filtered and RPM normalised')
dev.off()

pdf("plots/QC/filtered_rpkm.pdf")
qc.plots(rpkm.filtered, 'filtered and RPKM normalised')
dev.off()

pdf("plots/QC/filtered_tmm.pdf")
qc.plots(tmm.filtered, 'filtered and TMM (edgeR) normalised')
dev.off()

pdf("plots/QC/filtered_deseq.pdf")
qc.plots(deseq.filtered, 'filtered and DESeq normalised')
dev.off()

pdf("plots/QC/filtered_uq.pdf")
qc.plots(uq.filtered, 'filtered and Upper Quartile normalised')
dev.off()

pdf("plots/QC/filtered_quantile.pdf")
qc.plots(quantile.filtered, 'filtered and Quantile normalised')
dev.off()

# Alt Quantile ------------------------------------------------------------

#alt_quantile = data.frame(normalize.quantiles(as.matrix(unnormalised)),row.names=rownames(unnormalised))
#colnames(alt_quantile)=colnames(unnormalised)

#pdf("plots/QC/unfiltered_quantile_alt.pdf")
#qc.plots(alt_quantile, 'Unfiltered and Quantile normalised')
#dev.off()


# GGPAIRS -----------------------------------------------------------------

#unfiltered
pdf("plots/ggpairs/unfiltered_unnormalised.pdf")
ggpairs(unnormalised[,c('SRR7624365','SRR7624366','SRR7624371','SRR7624372')], title="Unfiltered and Unnormalised")
dev.off()

pdf("plots/ggpairs/unfiltered_rpm.pdf")
ggpairs(rpm[,c('SRR7624365','SRR7624366','SRR7624371','SRR7624372')], title="Unfiltered and RPM Normalised")
dev.off()

pdf("plots/ggpairs/unfiltered_rpkm.pdf")
ggpairs(rpkm[,c('SRR7624365','SRR7624366','SRR7624371','SRR7624372')], title="Unfiltered and RPKM Normalised")
dev.off()

pdf("plots/ggpairs/unfiltered_tmm.pdf")
ggpairs(tmm[,c('SRR7624365','SRR7624366','SRR7624371','SRR7624372')], title="Unfiltered and TMM (edgeR) Normalised")
dev.off()

pdf("plots/ggpairs/unfiltered_deseq.pdf")
ggpairs(deseq[,c('SRR7624365','SRR7624366','SRR7624371','SRR7624372')], title="Unfiltered and DESeq Normalised")
dev.off()

pdf("plots/ggpairs/unfiltered_uq.pdf")
ggpairs(uq[,c('SRR7624365','SRR7624366','SRR7624371','SRR7624372')], title="Unfiltered and Upper Quartile Normalised")
dev.off()

pdf("plots/ggpairs/unfiltered_quantile.pdf")
ggpairs(quantile[,c('SRR7624365','SRR7624366','SRR7624371','SRR7624372')], title="Unfiltered and Quantile Normalised")
dev.off()



#filtered
pdf("plots/ggpairs/filtered_unnormalised.pdf")
ggpairs(unnormalised.filtered[,c('SRR7624365','SRR7624366','SRR7624371','SRR7624372')], title="filtered and Unnormalised")
dev.off()

pdf("plots/ggpairs/filtered_rpm.pdf")
ggpairs(rpm.filtered[,c('SRR7624365','SRR7624366','SRR7624371','SRR7624372')], title="filtered and RPM Normalised")
dev.off()

pdf("plots/ggpairs/filtered_rpkm.pdf")
ggpairs(rpkm.filtered[,c('SRR7624365','SRR7624366','SRR7624371','SRR7624372')], title="filtered and RPKM Normalised")
dev.off()

pdf("plots/ggpairs/filtered_tmm.pdf")
ggpairs(tmm.filtered[,c('SRR7624365','SRR7624366','SRR7624371','SRR7624372')], title="filtered and TMM (edgeR) Normalised")
dev.off()

pdf("plots/ggpairs/filtered_deseq.pdf")
ggpairs(deseq.filtered[,c('SRR7624365','SRR7624366','SRR7624371','SRR7624372')], title="filtered and DESeq Normalised")
dev.off()

pdf("plots/ggpairs/filtered_uq.pdf")
ggpairs(uq.filtered[,c('SRR7624365','SRR7624366','SRR7624371','SRR7624372')], title="filtered and Upper Quartile Normalised")
dev.off()

pdf("plots/ggpairs/filtered_quantile.pdf")
ggpairs(quantile.filtered[,c('SRR7624365','SRR7624366','SRR7624371','SRR7624372')], title="filtered and Quantile Normalised")
dev.off()

#Comparing replicates
pdf("plots/replicates/unnormalised_filtered.pdf")
grid.arrange(ggplot(unnormalised.filtered, aes(x=SRR7624365, y=SRR7624366)) +
               geom_point() +
               xlab('0h rep 1 abundance') +
               ylab('0h rep 2 abundance'),
             ggplot(unnormalised.filtered, aes(x=log2(SRR7624365), y=log2(SRR7624366))) +
               geom_point() +
               xlab('0h rep 1 log2(abundance)') +
               ylab('0h rep 2 log2(abundance'),nrow=1,top='(Unnormalised) 0h, rep vs rep')

grid.arrange(ggplot(unnormalised.filtered, aes(x=SRR7624371, y=SRR7624372)) +
               geom_point() +
               xlab('12h rep 1 abundance') +
               ylab('12h rep 2 abundance'),
             ggplot(unnormalised.filtered, aes(x=log2(SRR7624371), y=log2(SRR7624372))) +
               geom_point() +
               xlab('12h rep 1 log2(abundance)') +
               ylab('12h rep 2 log2(abundance'),nrow=1,top='(Unnormalised) 12h, rep vs rep')
dev.off()


# MA plots ----------------------------------------------------------------

ma.plot = function(counts, meta, i, j, lower.lim=NA, upper.lim=NA, log.transformed=FALSE,title ){
  main.title = paste0(meta$id[i], ' v ', meta$id[j])
  main.title = title
  sub.title = paste(meta$type[i], meta$rep[i], 'v',
                    meta$type[j], meta$rep[j])
  # if already logtransformed we don't need to do it again
  if (log.transformed == TRUE){
    l1 = counts[,i]
    l2 = counts[,j]
  } else {
    # mask away the zeros
    zero.mask = !(counts[,i] == 0 | counts[,j] == 0)
    l1 = log2(counts[zero.mask, i])
    l2 = log2(counts[zero.mask, j])
  }
  m = l1 - l2
  a = 0.5 * (l1 + l2)
  data = data.frame(A = a, M = m)
  p = ggplot(data=data, aes(x=A, y=M, color='red', fill='red')) +
    geom_point(alpha=0.05)+
    theme(legend.position = "none") +
    geom_hline(yintercept=0.5,colour='gray60')+
    geom_hline(yintercept=-0.5,colour='gray60')+
    geom_hline(yintercept=1,colour='gray60')+
    geom_hline(yintercept=-1,colour='gray60')
  a.binned = cut(a, 20)
  data.binned = data.frame(A = a.binned, M = m)
  q = ggplot(data=data.binned) +
    geom_boxplot(aes(A, M)) +
    theme(axis.text.x=element_text(angle=90))+
    theme(legend.position = "none") +
    geom_hline(yintercept=0.5,colour='gray60')+
    geom_hline(yintercept=-0.5,colour='gray60')+
    geom_hline(yintercept=1,colour='gray60')+
    geom_hline(yintercept=-1,colour='gray60')
  # add ylim only if one of upper.lim, lower.lim is non-NA
  if (!is.na(lower.lim) | !is.na(upper.lim)){
    p = p + ylim(lower.lim, upper.lim)
    q = q + ylim(lower.lim, upper.lim)
  }
  grid.arrange(p, q, ncol = 2, top=paste0(main.title, '\n', sub.title))
}

#Filtered
pdf("plots/maplots/unnormalised_filtered.pdf")
for (i in 1:2){
  ma.plot(unnormalised.filtered[,c('SRR7624365','SRR7624366','SRR7624371','SRR7624372')],
          meta, (2*i)-1, 2*i, lower.lim= -5, upper.lim=5, title="Unnormalised (filtered)")
}
dev.off()

pdf("plots/maplots/rpm_filtered.pdf")
for (i in 1:2){
  ma.plot(rpm.filtered[,c('SRR7624365','SRR7624366','SRR7624371','SRR7624372')],
          meta, (2*i)-1, 2*i, lower.lim= -5, upper.lim=5, title="RPM Normalised (filtered)")
}
dev.off()

pdf("plots/maplots/rpkm_filtered.pdf")
for (i in 1:2){
  ma.plot(rpkm.filtered[,c('SRR7624365','SRR7624366','SRR7624371','SRR7624372')],
          meta, (2*i)-1, 2*i, lower.lim= -5, upper.lim=5, title="RPKM Normalised (filtered)")
}
dev.off()

pdf("plots/maplots/tmm_filtered.pdf")
for (i in 1:2){
  ma.plot(tmm.filtered[,c('SRR7624365','SRR7624366','SRR7624371','SRR7624372')],
          meta, (2*i)-1, 2*i, lower.lim= -5, upper.lim=5, title="TMM (edgeR) Normalised (filtered)")
}
dev.off()

pdf("plots/maplots/deseq_filtered.pdf")
for (i in 1:2){
  ma.plot(deseq.filtered[,c('SRR7624365','SRR7624366','SRR7624371','SRR7624372')],
          meta, (2*i)-1, 2*i, lower.lim= -5, upper.lim=5, title="DESeq Normalised (filtered)")
}
dev.off()

pdf("plots/maplots/uq_filtered.pdf")
for (i in 1:2){
  ma.plot(uq.filtered[,c('SRR7624365','SRR7624366','SRR7624371','SRR7624372')],
          meta, (2*i)-1, 2*i, lower.lim= -5, upper.lim=5, title="Upper Quartile Normalised (filtered)")
}
dev.off()

pdf("plots/maplots/quantile_filtered.pdf")
for (i in 1:2){
  ma.plot(quantile.filtered[,c('SRR7624365','SRR7624366','SRR7624371','SRR7624372')],
          meta, (2*i)-1, 2*i, lower.lim= -5, upper.lim=5, title="Quantile Normalised (filtered)")
}
dev.off()


# PCA  --------------------------------------------------------------------


# JSI ---------------------------------------------------------------------

jaccard.index = function(a, b){
  if ((length(a) == 0) & (length(b) == 0)){
    return(1)
  } else{
    u = length(union(a,b))
    i = length(intersect(a,b))
    return(i/u)
  }
}
leaf.labels=paste(meta$barcode,meta$type,meta$patient,meta$rep,sep='_')
jaccard.heatmap = function(counts, n.abundant, labels, heading){
  # colnames_counts <- c(t(inner(meta$type, meta$n_replicate, paste, sep = "_")))
  colnames_counts=paste0(meta$patient,meta$type,meta$rep)
  labels=labels[order(colnames_counts)]
  # colnames(reorderedcounts)=colnames_counts
  counts=counts[,order(colnames_counts)]
  n.samples = ncol(counts)
  hm = matrix(nrow=n.samples, ncol=n.samples)
  hm[] = 0
  for (i in 1:n.samples){
    for (j in 1:i){
      11
      i.gene.indices = order(counts[,i], decreasing=TRUE)[1:n.abundant]
      j.gene.indices = order(counts[,j], decreasing=TRUE)[1:n.abundant]
      hm[i, j] = jaccard.index(i.gene.indices, j.gene.indices)
      hm[j, i] = hm[i, j]
    }
  }
  title = paste0('Jaccard index of ', n.abundant, ' most abundant genes',heading)
  aheatmap(hm,
           color='Greys',
           Rowv = NA,
           Colv = NA,
           labRow=labels,
           labCol=labels,
           main=title,
           breaks=c(0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1),
           treeheight=0)
}

pdf('plots/jsi/unnormalised_filtered.pdf')
n.abundances = c(2000, 1000, 500, 200, 100, 50)
for (n in n.abundances){
  jaccard.heatmap(unnormalised.filtered[,c('SRR7624365','SRR7624366','SRR7624371','SRR7624372')], n, leaf.labels, heading=" (Unnormalised)")
}
dev.off()

pdf('plots/jsi/rpm_filtered.pdf')
n.abundances = c(2000, 1000, 500, 200, 100, 50)
for (n in n.abundances){
  jaccard.heatmap(rpm.filtered[,c('SRR7624365','SRR7624366','SRR7624371','SRR7624372')], n, leaf.labels, heading=" (RPM Normalised)")
}
dev.off()

pdf('plots/jsi/rpkm_filtered.pdf')
n.abundances = c(2000, 1000, 500, 200, 100, 50)
for (n in n.abundances){
  jaccard.heatmap(rpkm.filtered[,c('SRR7624365','SRR7624366','SRR7624371','SRR7624372')], n, leaf.labels, heading=" (RPKM Normalised)")
}
dev.off()

pdf('plots/jsi/tmm_filtered.pdf')
n.abundances = c(2000, 1000, 500, 200, 100, 50)
for (n in n.abundances){
  jaccard.heatmap(tmm.filtered[,c('SRR7624365','SRR7624366','SRR7624371','SRR7624372')], n, leaf.labels, heading=" (TMM [edgeR] Normalised)")
}
dev.off()

pdf('plots/jsi/deseq_filtered.pdf')
n.abundances = c(2000, 1000, 500, 200, 100, 50)
for (n in n.abundances){
  jaccard.heatmap(deseq.filtered[,c('SRR7624365','SRR7624366','SRR7624371','SRR7624372')], n, leaf.labels, heading=" (DESeq Normalised)")
}
dev.off()

pdf('plots/jsi/uq_filtered.pdf')
n.abundances = c(2000, 1000, 500, 200, 100, 50)
for (n in n.abundances){
  jaccard.heatmap(uq.filtered[,c('SRR7624365','SRR7624366','SRR7624371','SRR7624372')], n, leaf.labels, heading=" (Upper Quartile Normalised)")
}
dev.off()

pdf('plots/jsi/quantile_filtered.pdf')
n.abundances = c(2000, 1000, 500, 200, 100, 50)
for (n in n.abundances){
  jaccard.heatmap(quantile.filtered[,c('SRR7624365','SRR7624366','SRR7624371','SRR7624372')], n, leaf.labels, heading=" (Quantile Normalised)")
}
dev.off()


# Differential Expression -------------------------------------------------

#edgeR - RPM

edg = DGEList(quantile.filtered[,c('SRR7624365','SRR7624366','SRR7624371','SRR7624372')],group=c(1,1,2,2))
edg = estimateDisp(edg)

exact =exactTest(edg,dispersion=0.2)
exact$table$adjustedp = p.adjust(exact$table$PValue,method='BH')
print(length(rownames(exact$table[abs(exact$table$logFC) > 0.5 & exact$table$adjustedp<0.05,])))

edger_genes = rownames(exact$table[abs(exact$table$logFC) > 0.5 & exact$table$adjustedp<0.05,])

#DESeq2 - RPM

ds = DESeqDataSetFromMatrix(countData = round(quantile.filtered[,c('SRR7624365','SRR7624366','SRR7624371','SRR7624372')]),
                            colData = meta, design = ~ type)
ds = ds[rowSums(counts(ds)) > 0,]
ds$sizeFactor=rep(1:1,ncol(quantile.filtered[,c('SRR7624365','SRR7624366','SRR7624371','SRR7624372')]))
ds = DESeq(ds)

cts.ma = log2(counts(ds, normalized=TRUE) + 1)
cts.ds2.norm = as.data.frame(cts.ma)
ds.results = na.omit(results(ds, contrast=c('type','12h','0h'),
                             lfcThreshold=1.5,
                             altHypothesis='greaterAbs',
                             independentFiltering=TRUE ,
                             alpha=0.05))
deseq_genes = rownames(ds.results[abs(ds.results$log2FoldChange)>0.5 & ds.results$padj<0.05,])
print(length(deseq_genes))

#Enrichment

gprofiler_results = gprofiler2::gost(intersect(edger_genes,deseq_genes),
                                     organism='mmusculus',
                                     custom_bg = rownames(quantile.filtered[,c('SRR7624365','SRR7624366','SRR7624371','SRR7624372')]),
                                     sources=c('GO:BP','GO:MF','GO:CC','KEGG','REAC','TF','MIRNA'),
                                     correction_method='fdr')

pdf("plots/DE/quantile_filtered.pdf")
gostplot(gprofiler_results, capped = TRUE, interactive = FALSE)
dev.off()





