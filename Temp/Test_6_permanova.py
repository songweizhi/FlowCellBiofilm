import os
import numpy as np
import pandas as pd
from skbio import DistanceMatrix
from skbio.stats.distance import anosim
from skbio.stats.distance import permanova
from sklearn.metrics.pairwise import pairwise_distances


# input
wd = "/Users/songweizhi/Dropbox/Research/Flow_cell/deepSNV_output_summary/days"
infile_data = 'deepSNV_output_summary_D2_existence_D9_b.txt'
infile_factor = 'deepSNV_output_summary_D2_factor_D9_b.txt'
distance_metric = 'euclidean'


# set wd
os.chdir(wd)


# get sample_id_list and sample_group_list from input factor file
sample_id_list = []
sample_group_list = []
for each_sample in open(infile_factor):
    if not each_sample.startswith('Sample'):
        each_sample_split = each_sample.strip().split('\t')
        sample_id = each_sample_split[0]
        sample_group = each_sample_split[1]
        sample_id_list.append(sample_id)
        sample_group_list.append(sample_group)


# read in data as dataframe
df = pd.read_csv(infile_data, sep='\t')


# get list of list from dataframe
lol_data_in = []
for col_id in sample_id_list:
    column_num_list = (df[col_id].values).tolist()
    lol_data_in.append(column_num_list)


# calculate distance matrix
dist_arrary = pairwise_distances(lol_data_in, lol_data_in, metric=distance_metric)


# add sample id to distance matrix
dist_matrix = DistanceMatrix(dist_arrary, sample_id_list)


# perform anosim test
anosim_test = anosim(dist_matrix, sample_group_list, permutations=999)
print(anosim_test)
print()

# perform permanova test
permanova_test = permanova(dist_matrix, sample_group_list, permutations=999)
print(permanova_test)

