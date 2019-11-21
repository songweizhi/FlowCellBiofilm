library(ecodist) # for nmds,  https://www.rdocumentation.org/packages/ecodist/versions/2.0.1/topics/nmds
# set working directory
setwd("/Users/songweizhi/Dropbox/Research/Flow_cell_datasets/deepSNV_output_summary/days")
#strain = 'mono210' # mono210 or monoD2 or Coculture
#data_type = 'existence' # frequency or existence
#data_type = 'frequency' # frequency or existence
vegdist_method = 'euclidean' # euclidean or jaccard or bray
# get input file name
#snv_summary_file = paste('deepSNV_output_summary_', strain, '_', data_type, '.txt', sep = '')
#snv_summary_file
#snv_factor_file = paste('deepSNV_output_summary_', strain, '_factor.txt', sep = '')
#snv_factor_file
snv_summary = read.delim('deepSNV_output_summary_D2_existence_cdc_D9.txt', row.names=1)
snv_summary
snv_factor = read.delim('deepSNV_output_summary_D2_factor_D9.txt', row.names=1)
#snv_factor
# Import data
#snv_summary = read.delim(snv_summary_file, row.names=1)
#snv_factor = read.delim(snv_factor_file, row.names=1)
# transpose dataframe as vegan wants variables (SNVs) as columns
snv_summary_t = t(snv_summary)
#snv_summary_t
# calculate relative abundance
# snv_summary_t_norm = snv_summary_t / rowSums(snv_summary_t) * 100
## Calculate dissimilarity matrix, type ?vegdist for other options (e.g. Jaccard is useful for presence/absence analyses (composition only)).
snv_summary_t_ja = vegdist(snv_summary_t, method = vegdist_method)
snv_summary_t_ja
#class(snv_summary_t_ja)
#?vegdist
# PERMANOVA analysis (http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/adonis.html)
adonis2(snv_summary_t_ja ~ snv_factor$Species, snv_summary, method = vegdist_method)
#rm(list=ls())
# PERMANOVA analysis (http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/adonis.html)
snv_summary
library(vegan) # for metaMDS
library(ecodist) # for nmds,  https://www.rdocumentation.org/packages/ecodist/versions/2.0.1/topics/nmds
# set working directory
setwd("/Users/songweizhi/Dropbox/Research/Flow_cell/deepSNV_output_summary/days")
#strain = 'mono210' # mono210 or monoD2 or Coculture
#data_type = 'existence' # frequency or existence
#data_type = 'frequency' # frequency or existence
vegdist_method = 'euclidean' # euclidean or jaccard or bray
# get input file name
#snv_summary_file = paste('deepSNV_output_summary_', strain, '_', data_type, '.txt', sep = '')
#snv_summary_file
#snv_factor_file = paste('deepSNV_output_summary_', strain, '_factor.txt', sep = '')
#snv_factor_file
snv_summary = read.delim('deepSNV_output_summary_D2_existence_D9_b.txt', row.names=1)
#snv_summary
snv_factor = read.delim('deepSNV_output_summary_D2_factor_D9_b.txt', row.names=1)
#snv_factor
# Import data
#snv_summary = read.delim(snv_summary_file, row.names=1)
#snv_factor = read.delim(snv_factor_file, row.names=1)
# transpose dataframe as vegan wants variables (SNVs) as columns
snv_summary_t = t(snv_summary)
#snv_summary_t
# calculate relative abundance
# snv_summary_t_norm = snv_summary_t / rowSums(snv_summary_t) * 100
## Calculate dissimilarity matrix, type ?vegdist for other options (e.g. Jaccard is useful for presence/absence analyses (composition only)).
snv_summary_t_ja = vegdist(snv_summary_t, method = vegdist_method)
#snv_summary_t_ja
#class(snv_summary_t_ja)
#?vegdist
# PERMANOVA analysis (http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/adonis.html)
snv_summary
snv_summary_t_ja
adonis2(snv_summary_t_ja ~ snv_factor$Species, snv_summary, method = vegdist_method)
#rm(list=ls())
# PERMANOVA analysis (http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/adonis.html)
snv_summary
snv_summary_t_ja
adonis2(snv_summary_t_ja ~ snv_factor$Species, snv_summary, method = vegdist_method)
?adonis2
adonis2(snv_summary_t_ja ~ snv_factor$Species, snv_summary, method = vegdist_method, permutations = 99)
snv_summary_t_ja
adonis2(snv_summary_t_ja ~ snv_factor$Species, snv_summary, method = vegdist_method, permutations = 99)
library(vegan) # for metaMDS
library(ecodist) # for nmds,  https://www.rdocumentation.org/packages/ecodist/versions/2.0.1/topics/nmds
# set working directory
setwd("/Users/songweizhi/Dropbox/Research/Flow_cell/deepSNV_output_summary/days")
#strain = 'mono210' # mono210 or monoD2 or Coculture
#data_type = 'existence' # frequency or existence
#data_type = 'frequency' # frequency or existence
vegdist_method = 'euclidean' # euclidean or jaccard or bray
# get input file name
#snv_summary_file = paste('deepSNV_output_summary_', strain, '_', data_type, '.txt', sep = '')
#snv_summary_file
#snv_factor_file = paste('deepSNV_output_summary_', strain, '_factor.txt', sep = '')
#snv_factor_file
snv_summary = read.delim('deepSNV_output_summary_D2_existence_D9_b.txt', row.names=1)
#snv_summary
snv_factor = read.delim('deepSNV_output_summary_D2_factor_D9_b.txt', row.names=1)
#snv_factor
# Import data
#snv_summary = read.delim(snv_summary_file, row.names=1)
#snv_factor = read.delim(snv_factor_file, row.names=1)
# transpose dataframe as vegan wants variables (SNVs) as columns
snv_summary_t = t(snv_summary)
#snv_summary_t
# calculate relative abundance
# snv_summary_t_norm = snv_summary_t / rowSums(snv_summary_t) * 100
## Calculate dissimilarity matrix, type ?vegdist for other options (e.g. Jaccard is useful for presence/absence analyses (composition only)).
snv_summary_t_ja = vegdist(snv_summary_t, method = vegdist_method)
#snv_summary_t_ja
#class(snv_summary_t_ja)
#?vegdist
# PERMANOVA analysis (http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/adonis.html)
snv_summary
snv_summary_t_ja
adonis2(snv_summary_t_ja ~ snv_factor$Species, snv_summary, method = vegdist_method)
?adonis2
#rm(list=ls())
library(vegan) # for metaMDS
library(ecodist) # for nmds,  https://www.rdocumentation.org/packages/ecodist/versions/2.0.1/topics/nmds
# set working directory
setwd("/Users/songweizhi/Dropbox/Research/Flow_cell_datasets/deepSNV_output_summary/days")
#strain = 'mono210' # mono210 or monoD2 or Coculture
#data_type = 'existence' # frequency or existence
#data_type = 'frequency' # frequency or existence
vegdist_method = 'euclidean' # euclidean or jaccard or bray
# get input file name
#snv_summary_file = paste('deepSNV_output_summary_', strain, '_', data_type, '.txt', sep = '')
#snv_summary_file
#snv_factor_file = paste('deepSNV_output_summary_', strain, '_factor.txt', sep = '')
#snv_factor_file
snv_summary = read.delim('deepSNV_output_summary_D2_existence_D9_b.txt', row.names=1)
#snv_summary
snv_factor = read.delim('deepSNV_output_summary_D2_factor_D9_b.txt', row.names=1)
#snv_factor
# Import data
#snv_summary = read.delim(snv_summary_file, row.names=1)
#snv_factor = read.delim(snv_factor_file, row.names=1)
# transpose dataframe as vegan wants variables (SNVs) as columns
snv_summary_t = t(snv_summary)
#snv_summary_t
# calculate relative abundance
# snv_summary_t_norm = snv_summary_t / rowSums(snv_summary_t) * 100
## Calculate dissimilarity matrix, type ?vegdist for other options (e.g. Jaccard is useful for presence/absence analyses (composition only)).
snv_summary_t_ja = vegdist(snv_summary_t, method = vegdist_method)
#snv_summary_t_ja
#class(snv_summary_t_ja)
#?vegdist
# PERMANOVA analysis (http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/adonis.html)
snv_summary
snv_summary_t_ja
adonis2(snv_summary_t_ja ~ snv_factor$Species, snv_summary, method = vegdist_method)
?adonis2
#rm(list=ls())
library(vegan) # for metaMDS
library(ecodist) # for nmds,  https://www.rdocumentation.org/packages/ecodist/versions/2.0.1/topics/nmds
# set working directory
setwd("/Users/songweizhi/Dropbox/Research/Flow_cell/deepSNV_output_summary/days")
#strain = 'mono210' # mono210 or monoD2 or Coculture
#data_type = 'existence' # frequency or existence
#data_type = 'frequency' # frequency or existence
vegdist_method = 'euclidean' # euclidean or jaccard or bray
# get input file name
#snv_summary_file = paste('deepSNV_output_summary_', strain, '_', data_type, '.txt', sep = '')
#snv_summary_file
#snv_factor_file = paste('deepSNV_output_summary_', strain, '_factor.txt', sep = '')
#snv_factor_file
snv_summary = read.delim('deepSNV_output_summary_D2_existence_D9_b.txt', row.names=1)
#snv_summary
snv_factor = read.delim('deepSNV_output_summary_D2_factor_D9_b.txt', row.names=1)
#snv_factor
# Import data
#snv_summary = read.delim(snv_summary_file, row.names=1)
#snv_factor = read.delim(snv_factor_file, row.names=1)
# transpose dataframe as vegan wants variables (SNVs) as columns
snv_summary_t = t(snv_summary)
#snv_summary_t
# calculate relative abundance
# snv_summary_t_norm = snv_summary_t / rowSums(snv_summary_t) * 100
## Calculate dissimilarity matrix, type ?vegdist for other options (e.g. Jaccard is useful for presence/absence analyses (composition only)).
snv_summary_t_ja = vegdist(snv_summary_t, method = vegdist_method)
#snv_summary_t_ja
#class(snv_summary_t_ja)
#?vegdist
# PERMANOVA analysis (http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/adonis.html)
snv_summary
snv_summary_t_ja
adonis2(snv_summary_t_ja ~ snv_factor$Species, snv_summary, method = vegdist_method)
?adonis2
#rm(list=ls())
library(vegan) # for metaMDS
library(ecodist) # for nmds,  https://www.rdocumentation.org/packages/ecodist/versions/2.0.1/topics/nmds
# set working directory
setwd("/Users/songweizhi/Dropbox/Research/Flow_cell/deepSNV_output_summary/days")
#strain = 'mono210' # mono210 or monoD2 or Coculture
#data_type = 'existence' # frequency or existence
#data_type = 'frequency' # frequency or existence
vegdist_method = 'euclidean' # euclidean or jaccard or bray
# get input file name
#snv_summary_file = paste('deepSNV_output_summary_', strain, '_', data_type, '.txt', sep = '')
#snv_summary_file
#snv_factor_file = paste('deepSNV_output_summary_', strain, '_factor.txt', sep = '')
#snv_factor_file
snv_summary = read.delim('deepSNV_output_summary_D2_existence_D9_b.txt', row.names=1)
#snv_summary
snv_factor = read.delim('deepSNV_output_summary_D2_factor_D9_b.txt', row.names=1)
#snv_factor
# Import data
#snv_summary = read.delim(snv_summary_file, row.names=1)
#snv_factor = read.delim(snv_factor_file, row.names=1)
# transpose dataframe as vegan wants variables (SNVs) as columns
snv_summary_t = t(snv_summary)
#snv_summary_t
# calculate relative abundance
# snv_summary_t_norm = snv_summary_t / rowSums(snv_summary_t) * 100
## Calculate dissimilarity matrix, type ?vegdist for other options (e.g. Jaccard is useful for presence/absence analyses (composition only)).
snv_summary_t_ja = vegdist(snv_summary_t, method = vegdist_method)
#snv_summary_t_ja
#class(snv_summary_t_ja)
#?vegdist
# PERMANOVA analysis (http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/adonis.html)
snv_summary
snv_summary_t_ja
adonis2(snv_summary_t_ja ~ snv_factor$Species, snv_summary, method = vegdist_method)
?adonis2
rm(list=ls())
library(vegan) # for metaMDS
library(ecodist) # for nmds,  https://www.rdocumentation.org/packages/ecodist/versions/2.0.1/topics/nmds
# set working directory
setwd("/Users/songweizhi/Dropbox/Research/Flow_cell/deepSNV_output_summary/days")
#strain = 'mono210' # mono210 or monoD2 or Coculture
#data_type = 'existence' # frequency or existence
#data_type = 'frequency' # frequency or existence
vegdist_method = 'euclidean' # euclidean or jaccard or bray
# get input file name
#snv_summary_file = paste('deepSNV_output_summary_', strain, '_', data_type, '.txt', sep = '')
#snv_summary_file
#snv_factor_file = paste('deepSNV_output_summary_', strain, '_factor.txt', sep = '')
#snv_factor_file
snv_summary = read.delim('deepSNV_output_summary_D2_existence_D9_b.txt', row.names=1)
#snv_summary
snv_factor = read.delim('deepSNV_output_summary_D2_factor_D9_b.txt', row.names=1)
#snv_factor
# Import data
#snv_summary = read.delim(snv_summary_file, row.names=1)
#snv_factor = read.delim(snv_factor_file, row.names=1)
# transpose dataframe as vegan wants variables (SNVs) as columns
snv_summary_t = t(snv_summary)
#snv_summary_t
# calculate relative abundance
# snv_summary_t_norm = snv_summary_t / rowSums(snv_summary_t) * 100
## Calculate dissimilarity matrix, type ?vegdist for other options (e.g. Jaccard is useful for presence/absence analyses (composition only)).
snv_summary_t_ja = vegdist(snv_summary_t, method = vegdist_method)
#snv_summary_t_ja
#class(snv_summary_t_ja)
#?vegdist
# PERMANOVA analysis (http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/adonis.html)
snv_summary
snv_summary_t_ja
adonis2(snv_summary_t_ja ~ snv_factor$Species, snv_summary, method = vegdist_method)
?adonis2
rm(list=ls())
library(vegan) # for metaMDS
library(ecodist) # for nmds,  https://www.rdocumentation.org/packages/ecodist/versions/2.0.1/topics/nmds
# set working directory
setwd("/Users/songweizhi/Dropbox/Research/Flow_cell/deepSNV_output_summary/days")
#strain = 'mono210' # mono210 or monoD2 or Coculture
#data_type = 'existence' # frequency or existence
#data_type = 'frequency' # frequency or existence
vegdist_method = 'euclidean' # euclidean or jaccard or bray
# get input file name
#snv_summary_file = paste('deepSNV_output_summary_', strain, '_', data_type, '.txt', sep = '')
#snv_summary_file
#snv_factor_file = paste('deepSNV_output_summary_', strain, '_factor.txt', sep = '')
#snv_factor_file
snv_summary = read.delim('deepSNV_output_summary_D2_existence_D9_b.txt', row.names=1)
#snv_summary
snv_factor = read.delim('deepSNV_output_summary_D2_factor_D9_b.txt', row.names=1)
#snv_factor
# Import data
#snv_summary = read.delim(snv_summary_file, row.names=1)
#snv_factor = read.delim(snv_factor_file, row.names=1)
# transpose dataframe as vegan wants variables (SNVs) as columns
snv_summary_t = t(snv_summary)
#snv_summary_t
# calculate relative abundance
# snv_summary_t_norm = snv_summary_t / rowSums(snv_summary_t) * 100
## Calculate dissimilarity matrix, type ?vegdist for other options (e.g. Jaccard is useful for presence/absence analyses (composition only)).
snv_summary_t_ja = vegdist(snv_summary_t, method = vegdist_method)
#snv_summary_t_ja
#class(snv_summary_t_ja)
#?vegdist
# PERMANOVA analysis (http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/adonis.html)
snv_summary
snv_summary_t_ja
adonis2(snv_summary_t_ja ~ snv_factor$Species, snv_summary, method = vegdist_method)
?adonis2
rm(list=ls())
?adonis2
library(vegan) # for metaMDS
library(ecodist) # for nmds,  https://www.rdocumentation.org/packages/ecodist/versions/2.0.1/topics/nmds
# set working directory
setwd("/Users/songweizhi/Dropbox/Research/Flow_cell/deepSNV_output_summary/days")
#strain = 'mono210' # mono210 or monoD2 or Coculture
#data_type = 'existence' # frequency or existence
#data_type = 'frequency' # frequency or existence
vegdist_method = 'euclidean' # euclidean or jaccard or bray
# get input file name
#snv_summary_file = paste('deepSNV_output_summary_', strain, '_', data_type, '.txt', sep = '')
#snv_summary_file
#snv_factor_file = paste('deepSNV_output_summary_', strain, '_factor.txt', sep = '')
#snv_factor_file
snv_summary = read.delim('deepSNV_output_summary_D2_existence_D9_b.txt', row.names=1)
#snv_summary
snv_factor = read.delim('deepSNV_output_summary_D2_factor_D9_b.txt', row.names=1)
#snv_factor
# Import data
#snv_summary = read.delim(snv_summary_file, row.names=1)
#snv_factor = read.delim(snv_factor_file, row.names=1)
# transpose dataframe as vegan wants variables (SNVs) as columns
snv_summary_t = t(snv_summary)
#snv_summary_t
# calculate relative abundance
# snv_summary_t_norm = snv_summary_t / rowSums(snv_summary_t) * 100
## Calculate dissimilarity matrix, type ?vegdist for other options (e.g. Jaccard is useful for presence/absence analyses (composition only)).
snv_summary_t_ja = vegdist(snv_summary_t, method = vegdist_method)
#snv_summary_t_ja
#class(snv_summary_t_ja)
#?vegdist
# PERMANOVA analysis (http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/adonis.html)
snv_summary
snv_summary_t_ja
adonis2(snv_summary_t_ja ~ snv_factor$Species, snv_summary, method = vegdist_method)
?adonis2
rm(list=ls())
snv_summary_t_ja
library(vegan) # for metaMDS
library(ecodist) # for nmds,  https://www.rdocumentation.org/packages/ecodist/versions/2.0.1/topics/nmds
# set working directory
setwd("/Users/songweizhi/Dropbox/Research/Flow_cell/deepSNV_output_summary/days")
#strain = 'mono210' # mono210 or monoD2 or Coculture
#data_type = 'existence' # frequency or existence
#data_type = 'frequency' # frequency or existence
vegdist_method = 'euclidean' # euclidean or jaccard or bray
# get input file name
#snv_summary_file = paste('deepSNV_output_summary_', strain, '_', data_type, '.txt', sep = '')
#snv_summary_file
#snv_factor_file = paste('deepSNV_output_summary_', strain, '_factor.txt', sep = '')
#snv_factor_file
snv_summary = read.delim('deepSNV_output_summary_D2_existence_D9_b.txt', row.names=1)
#snv_summary
snv_factor = read.delim('deepSNV_output_summary_D2_factor_D9_b.txt', row.names=1)
#snv_factor
# Import data
#snv_summary = read.delim(snv_summary_file, row.names=1)
#snv_factor = read.delim(snv_factor_file, row.names=1)
# transpose dataframe as vegan wants variables (SNVs) as columns
snv_summary_t = t(snv_summary)
#snv_summary_t
# calculate relative abundance
# snv_summary_t_norm = snv_summary_t / rowSums(snv_summary_t) * 100
## Calculate dissimilarity matrix, type ?vegdist for other options (e.g. Jaccard is useful for presence/absence analyses (composition only)).
snv_summary_t_ja = vegdist(snv_summary_t, method = vegdist_method)
#snv_summary_t_ja
#class(snv_summary_t_ja)
#?vegdist
# PERMANOVA analysis (http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/adonis.html)
snv_summary
snv_summary_t_ja
adonis2(snv_summary_t_ja ~ snv_factor$Species, snv_summary, method = vegdist_method)
?adonis2
rm(list=ls())
library(vegan) # for metaMDS
library(ecodist) # for nmds,  https://www.rdocumentation.org/packages/ecodist/versions/2.0.1/topics/nmds
# set working directory
setwd("/Users/songweizhi/Dropbox/Research/Flow_cell/deepSNV_output_summary/days")
#strain = 'mono210' # mono210 or monoD2 or Coculture
#data_type = 'existence' # frequency or existence
#data_type = 'frequency' # frequency or existence
vegdist_method = 'euclidean' # euclidean or jaccard or bray
# get input file name
#snv_summary_file = paste('deepSNV_output_summary_', strain, '_', data_type, '.txt', sep = '')
#snv_summary_file
#snv_factor_file = paste('deepSNV_output_summary_', strain, '_factor.txt', sep = '')
#snv_factor_file
snv_summary = read.delim('deepSNV_output_summary_D2_existence_D9_b.txt', row.names=1)
#snv_summary
snv_factor = read.delim('deepSNV_output_summary_D2_factor_D9_b.txt', row.names=1)
#snv_factor
# Import data
#snv_summary = read.delim(snv_summary_file, row.names=1)
#snv_factor = read.delim(snv_factor_file, row.names=1)
# transpose dataframe as vegan wants variables (SNVs) as columns
snv_summary_t = t(snv_summary)
#snv_summary_t
# calculate relative abundance
# snv_summary_t_norm = snv_summary_t / rowSums(snv_summary_t) * 100
## Calculate dissimilarity matrix, type ?vegdist for other options (e.g. Jaccard is useful for presence/absence analyses (composition only)).
snv_summary_t_ja = vegdist(snv_summary_t, method = vegdist_method)
#snv_summary_t_ja
#class(snv_summary_t_ja)
#?vegdist
# PERMANOVA analysis (http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/adonis.html)
snv_summary
snv_summary_t_ja
adonis2(snv_summary_t_ja ~ snv_factor$Species, snv_summary, method = vegdist_method)
?adonis2
rm(list=ls())
library(vegan) # for metaMDS
library(ecodist) # for nmds,  https://www.rdocumentation.org/packages/ecodist/versions/2.0.1/topics/nmds
# set working directory
setwd("/Users/songweizhi/Dropbox/Research/Flow_cell/deepSNV_output_summary/days")
#strain = 'mono210' # mono210 or monoD2 or Coculture
#data_type = 'existence' # frequency or existence
#data_type = 'frequency' # frequency or existence
vegdist_method = 'euclidean' # euclidean or jaccard or bray
# get input file name
#snv_summary_file = paste('deepSNV_output_summary_', strain, '_', data_type, '.txt', sep = '')
#snv_summary_file
#snv_factor_file = paste('deepSNV_output_summary_', strain, '_factor.txt', sep = '')
#snv_factor_file
snv_summary = read.delim('deepSNV_output_summary_D2_existence_D9_b.txt', row.names=1)
#snv_summary
snv_factor = read.delim('deepSNV_output_summary_D2_factor_D9_b.txt', row.names=1)
#snv_factor
# Import data
#snv_summary = read.delim(snv_summary_file, row.names=1)
#snv_factor = read.delim(snv_factor_file, row.names=1)
# transpose dataframe as vegan wants variables (SNVs) as columns
snv_summary_t = t(snv_summary)
#snv_summary_t
# calculate relative abundance
# snv_summary_t_norm = snv_summary_t / rowSums(snv_summary_t) * 100
## Calculate dissimilarity matrix, type ?vegdist for other options (e.g. Jaccard is useful for presence/absence analyses (composition only)).
snv_summary_t_ja = vegdist(snv_summary_t, method = vegdist_method)
#snv_summary_t_ja
#class(snv_summary_t_ja)
#?vegdist
# PERMANOVA analysis (http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/adonis.html)
snv_summary
snv_summary_t_ja
adonis2(snv_summary_t_ja ~ snv_factor$Species, snv_summary, method = vegdist_method)
?adonis2
rm(list=ls())
?adonis2
library(vegan) # for metaMDS
library(ecodist) # for nmds,  https://www.rdocumentation.org/packages/ecodist/versions/2.0.1/topics/nmds
# set working directory
setwd("/Users/songweizhi/Dropbox/Research/Flow_cell/deepSNV_output_summary/days")
#strain = 'mono210' # mono210 or monoD2 or Coculture
#data_type = 'existence' # frequency or existence
#data_type = 'frequency' # frequency or existence
vegdist_method = 'euclidean' # euclidean or jaccard or bray
# get input file name
#snv_summary_file = paste('deepSNV_output_summary_', strain, '_', data_type, '.txt', sep = '')
#snv_summary_file
#snv_factor_file = paste('deepSNV_output_summary_', strain, '_factor.txt', sep = '')
#snv_factor_file
snv_summary = read.delim('deepSNV_output_summary_D2_existence_D9_b.txt', row.names=1)
#snv_summary
snv_factor = read.delim('deepSNV_output_summary_D2_factor_D9_b.txt', row.names=1)
#snv_factor
# Import data
#snv_summary = read.delim(snv_summary_file, row.names=1)
#snv_factor = read.delim(snv_factor_file, row.names=1)
# transpose dataframe as vegan wants variables (SNVs) as columns
snv_summary_t = t(snv_summary)
#snv_summary_t
# calculate relative abundance
# snv_summary_t_norm = snv_summary_t / rowSums(snv_summary_t) * 100
## Calculate dissimilarity matrix, type ?vegdist for other options (e.g. Jaccard is useful for presence/absence analyses (composition only)).
snv_summary_t_ja = vegdist(snv_summary_t, method = vegdist_method)
#snv_summary_t_ja
#class(snv_summary_t_ja)
#?vegdist
# PERMANOVA analysis (http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/adonis.html)
snv_summary
snv_summary_t_ja
adonis2(snv_summary_t_ja ~ snv_factor$Species, snv_summary, method = vegdist_method)
?adonis2
rm(list=ls())
library(vegan) # for metaMDS
library(ecodist) # for nmds,  https://www.rdocumentation.org/packages/ecodist/versions/2.0.1/topics/nmds
# set working directory
setwd("/Users/songweizhi/Dropbox/Research/Flow_cell/deepSNV_output_summary/days")
#strain = 'mono210' # mono210 or monoD2 or Coculture
#data_type = 'existence' # frequency or existence
#data_type = 'frequency' # frequency or existence
vegdist_method = 'euclidean' # euclidean or jaccard or bray
# get input file name
#snv_summary_file = paste('deepSNV_output_summary_', strain, '_', data_type, '.txt', sep = '')
#snv_summary_file
#snv_factor_file = paste('deepSNV_output_summary_', strain, '_factor.txt', sep = '')
#snv_factor_file
snv_summary = read.delim('deepSNV_output_summary_D2_existence_D9_b.txt', row.names=1)
#snv_summary
snv_factor = read.delim('deepSNV_output_summary_D2_factor_D9_b.txt', row.names=1)
#snv_factor
# Import data
#snv_summary = read.delim(snv_summary_file, row.names=1)
#snv_factor = read.delim(snv_factor_file, row.names=1)
# transpose dataframe as vegan wants variables (SNVs) as columns
snv_summary_t = t(snv_summary)
#snv_summary_t
# calculate relative abundance
# snv_summary_t_norm = snv_summary_t / rowSums(snv_summary_t) * 100
## Calculate dissimilarity matrix, type ?vegdist for other options (e.g. Jaccard is useful for presence/absence analyses (composition only)).
snv_summary_t_ja = vegdist(snv_summary_t, method = vegdist_method)
#snv_summary_t_ja
#class(snv_summary_t_ja)
#?vegdist
# PERMANOVA analysis (http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/adonis.html)
snv_summary
snv_summary_t_ja
adonis2(snv_summary_t_ja ~ snv_factor$Species, snv_summary, method = vegdist_method)
?adonis2
rm(list=ls())
