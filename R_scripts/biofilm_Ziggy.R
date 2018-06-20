library(vegan)


setwd("/Users/songweizhi/Desktop/aaa")

Dat = read.delim('deepSNV_output_summary_210_existence_cdc_D9.txt', row.names=1)
Dat
Fac1 = read.delim('deepSNV_output_summary_210_factor_D9.txt', row.names=1)
Fac1

tDat = t(Dat)
tDat

Dat1.bc = vegdist(tDat, method = "bray")
#Dat1.bc = vegdist(tDat, method = "euclidean")


Dat1.hclus = hclust(Dat1.bc, method = "average") # Cluster the distances
par(pty = 's') # square plot
plot(Dat1.hclus)

tDat1.mds = metaMDS(Dat1.bc, autotransform = F, trace = F, trymax=50)
tDat1.mds # stress = 0.130458










