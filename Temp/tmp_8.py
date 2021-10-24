import os
import argparse
from Bio import SeqIO
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def cluster(data, maxgap, min_size):

    data.sort()
    groups = [[data[0]]]
    for x in data[1:]:
        if abs(x - groups[-1][-1]) <= maxgap:
            groups[-1].append(x)
        else:
            groups.append([x])

    groups_qualified = []
    for each_cluster in groups:
        if len(each_cluster) >= min_size:
            groups_qualified.append(each_cluster)

    return groups_qualified


paired_x_list = [18, 20, 21, 21, 18, 20, 21, 21, 39, 40, 41, 42, 39, 40, 41, 42, 59, 60, 61, 59, 60, 61, 39, 40, 41, 42, 39, 40, 41, 42]
paired_y_list = [29, 30, 31, 32, 32, 31, 29, 30, 49, 50, 51, 52, 50, 52, 49, 51, 69, 70, 71, 70, 69, 71, 69, 70, 71, 72, 70, 72, 69, 71]

clipping_x_list = [23, 25, 26, 26, 23, 25, 26, 26, 39, 40, 41, 40, 39, 40, 41, 42, 54, 55, 56, 57, 54, 55, 56, 57, 74, 75, 76, 74, 75, 76, 39, 40, 41, 40, 39, 40, 41, 42]
clipping_y_list = [14, 15, 16, 17, 17, 16, 14, 15, 51, 52, 51, 49, 50, 52, 49, 51, 34, 35, 36, 37, 35, 37, 34, 36, 54, 55, 56, 55, 54, 56, 81, 82, 81, 79, 80, 82, 79, 81]

paired_cluster_max_gap         = 3
paired_min_cluster_reads_num   = 5
clipping_cluster_max_gap       = 3
clipping_min_cluster_reads_num = 5

pwd_figure_paired   = '/Users/songweizhi/Desktop/plot_paired.pdf'
pwd_figure_clipping = '/Users/songweizhi/Desktop/plot_clipping.pdf'


########################################################################################################################


paired_x_list_clustered   = cluster(sorted(paired_x_list), paired_cluster_max_gap, paired_min_cluster_reads_num)
paired_y_list_clustered   = cluster(sorted(paired_y_list), paired_cluster_max_gap, paired_min_cluster_reads_num)
clipping_x_list_clustered = cluster(sorted(clipping_x_list), clipping_cluster_max_gap, clipping_min_cluster_reads_num)
clipping_y_list_clustered = cluster(sorted(clipping_y_list), clipping_cluster_max_gap, clipping_min_cluster_reads_num)

print('paired_x_list')
print(paired_x_list)
print(paired_x_list_clustered)
print()

print('paired_y_list')
print(paired_y_list)
print(paired_y_list_clustered)
print()

print('clipping_x_list')
print(clipping_x_list)
print(clipping_x_list_clustered)
print()

print('clipping_y_list')
print(clipping_y_list)
print(clipping_y_list_clustered)
print()


########################################################################################################################

print('Plotting')

fig = plt.figure(figsize=(10, 10))
plt.scatter(paired_x_list, paired_y_list, s=32, linewidths=0, marker='.', c='blue', alpha=1)

# set axis range
plt.xlim(0, 100)
plt.ylim(0, 100)

# add axis label
plt.xlabel('2.10 chromosome (bp)')
plt.ylabel('D2 chromosome (bp)')

plt.tight_layout()
plt.savefig(pwd_figure_paired)
plt.close()


fig2 = plt.figure(figsize=(10, 10))
plt.scatter(clipping_x_list, clipping_y_list, s=32, linewidths=0, marker='.', c='red', alpha=1)

# set axis range
plt.xlim(0, 100)
plt.ylim(0, 100)

# add axis label
plt.xlabel('2.10 chromosome (bp)')
plt.ylabel('D2 chromosome (bp)')

plt.tight_layout()
plt.savefig(pwd_figure_clipping)
plt.close()

print('Done!')