library(vegan) # for metaMDS
library(ecodist) # for nmds,  https://www.rdocumentation.org/packages/ecodist/versions/2.0.1/topics/nmds

# set working directory
#setwd("/Users/songweizhi/Desktop/222/OneStep_Stats")
setwd("/Users/songweizhi/FC_Biofilm/OneStep_Stats")
setwd("/Users/songweizhi/Desktop/Biofilm2021/OneStep_MinBoth_10_MinEach_1_StrandBias_10_DepthDiff_30")
setwd("/Users/songweizhi/Desktop/000")

vegdist_method = 'euclidean' # euclidean or jaccard or bray, manhattan
# Dissimilarity index, partial match to "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao" or "mahalanobis".
# Jaccard is useful for presence/absence analyses (composition only)

# get input file name
strain = '210' # 210 or D2
snv_summary_file = paste('SNV_QC_ncd_even_flk_depth_', strain, '_matrix_no_plasmid.txt', sep = '')
#snv_summary_file = paste('SNV_QC_ncd_even_flk_depth_', strain, '_matrix.txt', sep = '')
snv_summary_file = paste('SNV_QC_ncd_even_flk_depth_', strain, '_matrix_gene_level_no_plasmid.txt', sep = '')
snv_summary_file = paste('SNV_QC_ncd_even_flk_depth_', strain, '_matrix_gene_level.txt', sep = '')
snv_summary_file = paste('SNV_QC_ncd_even_flk_depth_', strain, '_matrix_no_plasmid.txt', sep = '')

snv_factor_file = paste('Stats_', strain, '_factor.txt', sep = '')

# plot out
plot_out = paste(snv_summary_file, '.', vegdist_method, '.pdf', sep = '')

# Import data
snv_summary = read.delim(snv_summary_file, row.names=1)
snv_factor = read.delim(snv_factor_file, row.names=1)

# transpose dataframe as vegan wants variables (SNVs) as columns
snv_summary_t = t(snv_summary)

# calculate relative abundance
#snv_summary_t_norm = snv_summary_t / rowSums(snv_summary_t) * 100

## Calculate dissimilarity matrix, type ?vegdist for other options (e.g. Jaccard is useful for presence/absence analyses (composition only)).
snv_summary_t_ja = vegdist(snv_summary_t, method = vegdist_method)

# plot hclust
# snv_summary_t_ja_hclus = hclust(snv_summary_t_ja, method = "average") # Cluster the distances
# par(pty = 's') # square plot
# plot(snv_summary_t_ja_hclus)
# plot nMDS

# non-metric multi-dimensional scaling (nMDS), which is an ordination method.
# A low stress value (<0.2) and a linear pattern in a stressplot suggest a good mapping of your data.
#snv_summary_t_ja_mds = metaMDS(snv_summary_t_ja, autotransform = F, trace = F, trymax=50)
snv_summary_t_ja_mds = metaMDS(comm = snv_summary_t_ja)
#snv_summary_t_ja_mds = nmds(snv_summary_t_ja)

# plot stress
# stressplot(snv_summary_t_ja_mds)

####################################### Get MDS plot #######################################

# get plot_shape_list, Mono210: square; MonoD2: triangle; Coculture: circle
plot_shape_list  = c()
for (eachA in snv_factor$Species){
  if (eachA == 'Mono210'){
    current_shape = 15
  } else if (eachA == 'MonoD2'){
    current_shape = 17
  } else if (eachA == 'Coculture'){
    current_shape = 19
  }
  plot_shape_list = c(plot_shape_list, current_shape)
}

# plot
pdf(plot_out, width=5, height=5, pointsize=9)
par(mar=c(5,4,2,7))
plot(snv_summary_t_ja_mds$points, col=snv_factor$Label, pch=plot_shape_list)
text(snv_summary_t_ja_mds$points, labels=snv_factor$Time, cex= 0.7, pos=3)

# get legend_shape_list
legend_shape_list = c()
legend_shape_list_210 = c(19, 19, 19, 15, 15, 15)
legend_shape_list_D2 = c(19, 19, 19, 17, 17, 17)
if (strain == '210'){
  legend_shape_list = legend_shape_list_210
}else if (strain == 'D2'){
  legend_shape_list = legend_shape_list_D2
}

# add legend#
#legend("topleft", inset=c(1,0), xpd=TRUE, bty="n", legend = levels(snv_factor$Label),pch=legend_shape_list, col=as.numeric(as.factor(levels(snv_factor$Label))))

invisible(dev.off())

############################## Perform PERMANOVA analysis ##############################

# PERMANOVA analysis (http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/adonis.html)
adonis2(snv_summary_t_ja ~ snv_factor$Species, snv_summary, method = vegdist_method)
adonis2(snv_summary_t_ja ~ snv_factor$Replicate, snv_summary, method = vegdist_method)
#adonis2(snv_summary_t_ja ~ snv_factor$Species*snv_factor$Replicate, snv_summary, method = vegdist_method)

# clear variables/objects
rm(list=ls())
