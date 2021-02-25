#Load the data
data = read.csv("data/raw/counts_raw.csv", header=TRUE)
data = as.matrix(data[,2:17])
#Load library
library(noisyr)

#Similarity calculation across samples
dist_ms <- calculate_distance_matrices_counts(data)

#Noise Quantification
thre_noise <- calculate_threshold_noise(expression.matrix=dist_ms$exp,
                                        abn.matrix=dist_ms$abn)

#Visualisation function
plot_distance_abundance(abn.matrix = dist_ms$abn, dist.matrix = dist_ms$dist)

#Noise removal
remove_noise_method(expression.matrix=dist_ms$exp, abn.matrix=dist_ms$abn,
                    dist.matrix=dist_ms$dist)
