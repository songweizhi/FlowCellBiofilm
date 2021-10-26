library(vegan)


setwd("/Users/songweizhi/Desktop/000")

Dat = read.delim('deepSNV_output_summary_210_existence_cdc_D9.txt', row.names=1)
Fac1 = read.delim('deepSNV_output_summary_210_factor_D9.txt', row.names=1)


snv_summary_file = '/Users/songweizhi/Desktop/000/SNV_QC_ncd_even_flk_depth_210_matrix_no_plasmid.txt'
snv_factor_file = '/Users/songweizhi/Desktop/000/Stats_210_factor.txt'

Dat = read.delim(snv_summary_file, row.names=1)
Fac1 = read.delim(snv_factor_file, row.names=1)


tDat = t(Dat)
tDat

Dat1.bc = vegdist(tDat, method = "bray")
#Dat1.bc = vegdist(tDat, method = "euclidean")


Dat1.hclus = hclust(Dat1.bc, method = "average") # Cluster the distances
par(pty = 's') # square plot
plot(Dat1.hclus)

tDat1.mds = metaMDS(Dat1.bc, autotransform = F, trace = F, trymax=50)
tDat1.mds # stress = 0.130458










