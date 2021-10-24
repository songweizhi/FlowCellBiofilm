# Multidimensional Scaling
# https://www.statmethods.net/advstats/mds.html

# MDS flow cell biofilm
setwd("/Users/songweizhi/Dropbox/Research/Flow_cell/deepSNV_output_summary")

## Import data
snv_summary = read.delim('deepSNV_output_summary_existence.txt', row.names=1)

snv_summary
d <- dist(snv_summary)
d
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
fit
# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
     main="Metric	MDS",	type="n")
text(x, y, labels = row.names(mydata), cex=.7)


library(MASS)
d <- dist(snv_summary) # euclidean distances between the rows
fit <- isoMDS(d, k=2) # k is the number of dim
fit # view results

# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
     main="Nonmetric MDS", type="n")
text(x, y, labels = row.names(mydata), cex=.7)
