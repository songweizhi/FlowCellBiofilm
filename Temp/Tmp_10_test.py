import argparse
import numpy as np
from Bio import SeqIO
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec


pwd_figure = '/Users/songweizhi/Desktop/test.pdf'

ref_length_dict = {'2.10_chromosome': 10, 'D2_c': 10}

clipping_pos_list_210_uniq = [1,2,3,4,5,6,7,8,9,0]
clipping_pos_list_210_uniq_count = [1,2,3,4,5,6,7,8,9,0]
max_count_210 = 10

clipping_pos_list_D2_uniq = [1,2,3,4,5,6,7,8,9,0]
clipping_pos_list_D2_uniq_count = [1,2,3,4,5,6,7,8,9,0]
max_count_D2 = 10

plot_num_list_dict_no_16s = {'2.10_chromosome': [1,2,3,4,5,6,7,8,9], 'D2_c': [1,2,3,4,5,6,7,8,9]}





print('Plotting')

fig = plt.figure(figsize=(10, 10))
gs = gridspec.GridSpec(2, 2, height_ratios=[1, 9], width_ratios=[9, 1])

ax0 = plt.subplot(gs[0])

ax0.bar(clipping_pos_list_210_uniq, clipping_pos_list_210_uniq_count, color='skyblue', width=1, linewidth=0)
plt.xlim(0, ref_length_dict['2.10_chromosome'])
plt.ylim((0, max_count_210 + 1))
# hide x-axis
x_axis = ax0.axes.get_xaxis()
x_axis.set_visible(False)

ax2 = plt.subplot(gs[2])
ax2.scatter(plot_num_list_dict_no_16s['2.10_chromosome'], plot_num_list_dict_no_16s['D2_c'], s=5, linewidths=0, marker='o')
plt.xlim(0, ref_length_dict['2.10_chromosome'])
plt.ylim(0, ref_length_dict['D2_c'])


ax3 = plt.subplot(gs[3])
ax3.barh(clipping_pos_list_D2_uniq, clipping_pos_list_D2_uniq_count, color='skyblue', linewidth=0)
plt.ylim(0, ref_length_dict['D2_c'])
plt.xlim((0, max_count_D2 + 1))
# hide x-axis
y_axis = ax3.axes.get_yaxis()
y_axis.set_visible(False)

plt.tight_layout()
plt.savefig(pwd_figure)
plt.close()

print('Done!')
