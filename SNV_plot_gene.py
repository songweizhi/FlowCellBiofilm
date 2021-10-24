import math
import random
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns


def get_color_list(color_num):

    if color_num <= 8:
        color_list_combined = ['#3787c0', '#39399f', '#ffb939', '#399f39', '#9f399f', '#fb694a', '#9f9f39', '#959595']

    elif 8 < color_num <= 16:
        color_list_combined = ['#2b7bba', '#89bedc', '#2e2e99', '#8a8acc', '#ffa500', '#ffc55c', '#2e992e', '#8acc8a', '#992e99', '#cc8acc', '#d52221', '#fc8161', '#99992e', '#cccc8a', '#5c5c5c', '#adadad']

    else:
        color_num_each = math.ceil(color_num/8) + 2

        color_list_1 = sns.color_palette('Blues',  n_colors=color_num_each).as_hex()
        color_list_2 = sns.light_palette('navy',   n_colors=color_num_each).as_hex()
        color_list_3 = sns.light_palette('orange', n_colors=color_num_each).as_hex()
        color_list_4 = sns.light_palette('green',  n_colors=color_num_each).as_hex()
        color_list_5 = sns.light_palette('purple', n_colors=color_num_each).as_hex()
        color_list_6 = sns.color_palette('Reds',   n_colors=color_num_each).as_hex()
        color_list_7 = sns.light_palette('olive',  n_colors=color_num_each).as_hex()
        color_list_8 = sns.color_palette('Greys', n_colors=color_num_each).as_hex()

        color_list_combined = []
        for color_list in [color_list_1, color_list_2, color_list_3, color_list_4, color_list_5, color_list_6, color_list_7, color_list_8]:
            for color in color_list[2:][::-1]:
                color_list_combined.append(color)

    color_list_to_return = random.sample(color_list_combined, color_num)

    color_list_to_return_sorted = []
    for color_to_return in color_list_combined:
        if color_to_return in color_list_to_return:
            color_list_to_return_sorted.append(color_to_return)

    return color_list_to_return_sorted


# file in
deepSNV_output_summary_210_frequency_file   = '/Users/songweizhi/FC_Biofilm/backup/Flow_cell_datasets/deepSNV_summary_subsampled/deepSNV_output_summary_210_frequency.txt'
deepSNV_output_summary_D2_frequency_file    = '/Users/songweizhi/FC_Biofilm/backup/Flow_cell_datasets/deepSNV_summary_subsampled/deepSNV_output_summary_D2_frequency.txt'
ref                                         = '210'

# file out
output_plot_210                             = '/Users/songweizhi/Desktop/000/snp_freq_plot_210.pdf'
output_plot_D2                              = '/Users/songweizhi/Desktop/000/snp_freq_plot_D2.pdf'


# plot parameters
default_color       = 'lightgrey'
sample_color_A      = 'deepskyblue'
sample_color_B      = 'coral'
sample_color_C      = 'limegreen'
label_list          = ['D0', 'D9', 'D18', 'D27', 'D42']
highlight_color     = 'skyblue'
dafult_linewidth    = 0.5


deepSNV_output_summary_file = ''
output_plot = ''
if ref == '210':
    deepSNV_output_summary_file = deepSNV_output_summary_210_frequency_file
    output_plot = output_plot_210
if ref == 'D2':
    deepSNV_output_summary_file = deepSNV_output_summary_D2_frequency_file
    output_plot = output_plot_D2


# create a blank graph
fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(nrows=3, ncols=2, figsize=(15, 10))  # figsize=(width, height)
color_list = get_color_list(8)

num_index = list(range(1, len(label_list) + 1))
color_index = 0
for snp in open(deepSNV_output_summary_file):
    if not snp.startswith('\t'):
        snp_210_split = snp.strip().split('\t')
        snp_id = snp_210_split[0]

        print(snp_id)

        # monoculture
        mono_snp_freq_rep_1 = [0.0] + [float(i) for i in snp_210_split[1:5]]
        mono_snp_freq_rep_2 = [0.0] + [float(i) for i in snp_210_split[5:9]]
        mono_snp_freq_rep_3 = [0.0] + [float(i) for i in snp_210_split[9:13]]

        # coculture
        co_snp_freq_rep_1 = [0.0] + [float(i) for i in snp_210_split[13:17]]
        co_snp_freq_rep_2 = [0.0] + [float(i) for i in snp_210_split[17:21]]
        co_snp_freq_rep_3 = [0.0] + [float(i) for i in snp_210_split[21:25]]

        freq_lol        = [mono_snp_freq_rep_1, mono_snp_freq_rep_2, mono_snp_freq_rep_3, co_snp_freq_rep_1, co_snp_freq_rep_2, co_snp_freq_rep_3]
        freq_mono_lol   = [mono_snp_freq_rep_1, mono_snp_freq_rep_2, mono_snp_freq_rep_3]
        freq_co_lol     = [co_snp_freq_rep_1, co_snp_freq_rep_2, co_snp_freq_rep_3]

        if snp_id == '2.10_chromosome|3347922|C|T':

            # monoculture
            ax1.plot(num_index, mono_snp_freq_rep_1, color=default_color, linewidth=dafult_linewidth)
            ax3.plot(num_index, mono_snp_freq_rep_2, color=default_color, linewidth=dafult_linewidth)
            ax5.plot(num_index, mono_snp_freq_rep_3, color=default_color, linewidth=dafult_linewidth)

            # coculture
            ax2.plot(num_index, co_snp_freq_rep_1, color=default_color, linewidth=dafult_linewidth)
            ax4.plot(num_index, co_snp_freq_rep_2, color=default_color, linewidth=dafult_linewidth)
            ax6.plot(num_index, co_snp_freq_rep_3, color=default_color, linewidth=dafult_linewidth)


# add title for subplot
ax1.title.set_text('%s Monoculture A' % ref)
ax2.title.set_text('%s Coculture A' % ref)
ax3.title.set_text('%s Monoculture B' % ref)
ax4.title.set_text('%s Coculture B' % ref)
ax5.title.set_text('%s Monoculture C' % ref)
ax6.title.set_text('%s Coculture C' % ref)

plt.setp(((ax1, ax2), (ax3, ax4), (ax5, ax6)), xticks=num_index, xticklabels=label_list)
plt.tight_layout()
fig.savefig(output_plot, dpi=300)
plt.close()

