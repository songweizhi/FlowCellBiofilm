import argparse
import numpy as np
from Bio import SeqIO
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec


def get_16S_pos_list(gbk_file, flk_length):

    ribosomal_RNA_16S_regions = set()
    for seq_record in SeqIO.parse(gbk_file, 'genbank'):

        for gene_record in seq_record.features:
            if 'product' in gene_record.qualifiers:
                gene_product = gene_record.qualifiers['product'][0]
                if gene_product == '16S ribosomal RNA':
                    start_pos = gene_record.location.start - flk_length
                    end_pos = gene_record.location.end + flk_length

                    # add 16S pos to list
                    for i in range(start_pos, end_pos):
                        ribosomal_RNA_16S_regions.add(i)

    return ribosomal_RNA_16S_regions


def cluster(data, maxgap):

    data.sort()
    groups = [[data[0]]]
    for x in data[1:]:
        if abs(x - groups[-1][-1]) <= maxgap:
            groups[-1].append(x)
        else:
            groups.append([x])

    return groups


parser = argparse.ArgumentParser()

parser.add_argument('-r',        required=True,                          help='reference')
parser.add_argument('-s',        required=True,                          help='sam file')
parser.add_argument('-g',        required=True, type=int,                help='maximum gap')
parser.add_argument('-n',        required=True, type=int,                help='minimum number of reads in a cluster')
parser.add_argument('-o',        required=True,                          help='output')
parser.add_argument('-flk16s',   required=False, type=int, default=1000, help='ignore paired reads in 16s regions')
parser.add_argument('-no16s',    required=False, action='store_true',    help='ignore paired reads in 16s regions')
parser.add_argument('-hotspots', required=True,                          help='ignore paired reads in 16s regions')


args = vars(parser.parse_args())
ref_file                = args['r']
sam_file                = args['s']
cluster_max_gap         = args['g']
min_cluster_reads_num   = args['n']
pwd_figure              = args['o']
no16s                   = args['no16s']
flk_length_16S          = args['flk16s']
hotspots_location_file  = args['hotspots']
plot_hotspot = False
skip_parse_sam = True
sam_file_basename = '.'.join(sam_file.split('/')[-1].split('.')[:-1])

gbk_file_210            = '/Users/songweizhi/Desktop/mapHGT_biofilm/file_in/210.gbk'
gbk_file_D2             = '/Users/songweizhi/Desktop/mapHGT_biofilm/file_in/D2.gbk'
paired_reads_ref_file   = '/Users/songweizhi/Desktop/mapHGT_biofilm/%s_paired_reads_ref.txt' % sam_file_basename


# get 16S regions
print('get 16S regions')
ribosomal_RNA_16S_pos_list_210 = get_16S_pos_list(gbk_file_210, flk_length_16S)
ribosomal_RNA_16S_pos_list_D2  = get_16S_pos_list(gbk_file_D2, flk_length_16S)
ribosomal_RNA_16S_pos_list_dict = {'2.10_chromosome': ribosomal_RNA_16S_pos_list_210, 'D2_c': ribosomal_RNA_16S_pos_list_D2}


# get length of reference genomes
ref_length_dict = {}
ref_id_list = []
for seq_record in SeqIO.parse(ref_file, 'fasta'):
    ref_length_dict[seq_record.id] = len(seq_record.seq)
    ref_id_list.append(seq_record.id)


# read in sam file
if skip_parse_sam is False:
    print('read in sam file')
    paired_reads_ref_dict = {}
    for each_read in open(sam_file):

        if not each_read.startswith('@'):
            each_read_split = each_read.strip().split('\t')
            read_id         = each_read_split[0]
            read_id_base    = '_'.join(read_id.split('_')[:-1])
            read_strand     = read_id.split('_')[-1]
            ref_id          = each_read_split[2]
            ref_pos         = int(each_read_split[3])
            cigar           = each_read_split[5]
            read_seq        = each_read_split[9]
            read_len        = len(read_seq)

            if read_id_base not in paired_reads_ref_dict:
                if read_strand == '1':
                    paired_reads_ref_dict[read_id_base] = [ref_id, str(ref_pos), '', '']
                if read_strand == '2':
                    paired_reads_ref_dict[read_id_base] = ['', '', ref_id, str(ref_pos)]
            else:
                if read_strand == '1':
                    paired_reads_ref_dict[read_id_base][0] = ref_id
                    paired_reads_ref_dict[read_id_base][1] = str(ref_pos)
                if read_strand == '2':
                    paired_reads_ref_dict[read_id_base][2] = ref_id
                    paired_reads_ref_dict[read_id_base][3] = str(ref_pos)

    # write out to file
    paired_reads_ref_file_handle = open(paired_reads_ref_file, 'w')
    for paired_reads_ref in paired_reads_ref_dict:
        paired_reads_ref_value = paired_reads_ref_dict[paired_reads_ref]
        if paired_reads_ref_value[0] != paired_reads_ref_value[2]:
            paired_reads_ref_file_handle.write('%s\t%s\n' % (paired_reads_ref, '\t'.join(paired_reads_ref_dict[paired_reads_ref])))
    paired_reads_ref_file_handle.close()


print('store paired reads in dict')
ref_pos_dict = {'2.10_chromosome': [], 'D2_c': []}
for paired_reads in open(paired_reads_ref_file):
    paired_reads_split    = paired_reads.strip().split('\t')
    paired_read_1_ref     = paired_reads_split[1]
    paired_read_1_ref_pos = int(paired_reads_split[2])
    paired_read_2_ref     = paired_reads_split[3]
    paired_read_2_ref_pos = int(paired_reads_split[4])
    if paired_read_1_ref != paired_read_2_ref:
        ref_pos_dict[paired_read_1_ref].append(paired_read_1_ref_pos)
        ref_pos_dict[paired_read_2_ref].append(paired_read_2_ref_pos)


# do clustering and write out clustered locations
print('do clustering')
ref_pos_dict_clustered = {'2.10_chromosome': [], 'D2_c': []}
hotspots_location_handle = open(hotspots_location_file, 'w')
for ref in ref_pos_dict:

    # cluster pos list
    ref_pos_list = sorted(ref_pos_dict[ref])
    ref_pos_list_clustered = cluster(ref_pos_list, cluster_max_gap)

    # only keep the clustered pos
    for each_cluster in ref_pos_list_clustered:
        if len(each_cluster) >= min_cluster_reads_num:
            hotspots_location_handle.write('%s\t%s\n' % (ref, ','.join([str(i) for i in each_cluster])))
            for each_pos in each_cluster:
                ref_pos_dict_clustered[ref].append(each_pos)
hotspots_location_handle.close()


print('store paired reads in dict')
plot_num_list_dict = {'2.10_chromosome': [], 'D2_c': []}
for paired_reads in open(paired_reads_ref_file):
    paired_reads_split    = paired_reads.strip().split('\t')
    paired_read_1_ref     = paired_reads_split[1]
    paired_read_1_ref_pos = int(paired_reads_split[2])
    paired_read_2_ref     = paired_reads_split[3]
    paired_read_2_ref_pos = int(paired_reads_split[4])
    if (paired_read_1_ref != paired_read_2_ref) and (paired_read_1_ref_pos in ref_pos_dict_clustered[paired_read_1_ref]) and (paired_read_2_ref_pos in ref_pos_dict_clustered[paired_read_2_ref]):
        plot_num_list_dict[paired_read_1_ref].append(paired_read_1_ref_pos)
        plot_num_list_dict[paired_read_2_ref].append(paired_read_2_ref_pos)


# remove 16s regions
print('remove 16s regions')
plot_num_list_dict_no_16s = {'2.10_chromosome': [], 'D2_c': []}
for (x, y) in zip(plot_num_list_dict['2.10_chromosome'], plot_num_list_dict['D2_c']):
    if (x not in ribosomal_RNA_16S_pos_list_dict['2.10_chromosome']) and (y not in ribosomal_RNA_16S_pos_list_dict['D2_c']):
        plot_num_list_dict_no_16s['2.10_chromosome'].append(x)
        plot_num_list_dict_no_16s['D2_c'].append(y)


########################################################################################################################

clipping_pos_file_210   = '/Users/songweizhi/Desktop/mapHGT_biofilm/clipping_pos_210_VIP.txt'
clipping_pos_file_D2   = '/Users/songweizhi/Desktop/mapHGT_biofilm/clipping_pos_D2_VIP.txt'

clipping_pos_list_210 = []
for each_pos in open(clipping_pos_file_210):
    clipping_pos_list_210.append(int(each_pos.strip()))

clipping_pos_list_D2 = []
for each_pos in open(clipping_pos_file_D2):
    clipping_pos_list_D2.append(int(each_pos.strip()))


clipping_pos_list_210_uniq = []
clipping_pos_list_210_uniq_count = []
max_count_210 = 0
for i in sorted(clipping_pos_list_210):
    if i not in clipping_pos_list_210_uniq:
        clipping_pos_list_210_uniq.append(i)
        clipping_pos_list_210_uniq_count.append(clipping_pos_list_210.count(i))
        if clipping_pos_list_210.count(i) > max_count_210:
            max_count_210 = clipping_pos_list_210.count(i)


clipping_pos_list_D2_uniq = []
clipping_pos_list_D2_uniq_count = []
max_count_D2 = 0
for i in sorted(clipping_pos_list_D2):
    if i not in clipping_pos_list_D2_uniq:
        clipping_pos_list_D2_uniq.append(i)
        clipping_pos_list_D2_uniq_count.append(clipping_pos_list_D2.count(i))
        if clipping_pos_list_D2.count(i) > max_count_D2:
            max_count_D2 = clipping_pos_list_D2.count(i)


########################################################################################################################

print('Plotting')

fig = plt.figure(figsize=(10, 10))
gs = gridspec.GridSpec(2, 2, height_ratios=[1, 9], width_ratios=[9, 1])
ax0 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[2])
ax3 = plt.subplot(gs[3])


# ax0
ax0.bar(clipping_pos_list_210_uniq, clipping_pos_list_210_uniq_count, color='skyblue', width=1, linewidth=0)
plt.xlim(0, ref_length_dict['2.10_chromosome'])
plt.ylim((0, max_count_210 + 1))
# hide x-axis
x_axis = ax0.axes.get_xaxis()
x_axis.set_visible(False)

# ax2
ax2.scatter(plot_num_list_dict_no_16s['2.10_chromosome'], plot_num_list_dict_no_16s['D2_c'], s=5, linewidths=0, marker='o')
plt.xlim(0, ref_length_dict['2.10_chromosome'])
plt.ylim(0, ref_length_dict['D2_c'])

# ax3
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


########################################################################################################################

# get the following two files
# clipping_pos_file_210   = '/Users/songweizhi/Desktop/mapHGT_biofilm/clipping_pos_210_VIP.txt'
# clipping_pos_file_D2   = '/Users/songweizhi/Desktop/mapHGT_biofilm/clipping_pos_D2_VIP.txt'


# def cigar_splitter(cigar):
#
#     # get the position of letters
#     letter_pos_list = []
#     n = 0
#     for each_element in cigar:
#         if each_element.isalpha() is True:
#             letter_pos_list.append(n)
#         n += 1
#
#     # split cigar
#     index = 0
#     cigar_splitted = []
#     while index <= len(letter_pos_list) - 1:
#         if index == 0:
#             cigar_splitted.append(cigar[:(letter_pos_list[index] + 1)])
#         else:
#             cigar_splitted.append(cigar[(letter_pos_list[index - 1] + 1):(letter_pos_list[index] + 1)])
#         index += 1
#
#     return cigar_splitted
#
#
# def get_clipping_pos_list(sam_in):
#
#     clipping_pos_list_dict = {}
#
#     for each in open(sam_in):
#         if not each.startswith('@'):
#             each_split = each.strip().split('\t')
#             ref_id = each_split[2]
#             reads_seq = each_split[9]
#             reads_len = len(reads_seq)
#             full_length_M = '%sM' % str(reads_len)
#             cigar = each_split[5]
#
#             if cigar != full_length_M:
#                 if cigar == '*':
#                     pass
#                 else:
#                     cigar_splitted = cigar_splitter(cigar)
#
#                     # get the length of S region
#                     if (cigar_splitted[0][-1] == 'M') and (cigar_splitted[-1][-1] == 'S') and (len(cigar_splitted) == 2):
#                         len_S = int(cigar_splitted[-1][:-1])
#                         if len_S >= 10:
#                             clipping_pos = int(each_split[3]) + int(cigar_splitted[0][:-1])
#                             if ref_id not in clipping_pos_list_dict:
#                                 clipping_pos_list_dict[ref_id] = [clipping_pos]
#                             else:
#                                 clipping_pos_list_dict[ref_id].append(clipping_pos)
#
#                     if (cigar_splitted[0][-1] == 'S') and (cigar_splitted[-1][-1] == 'M') and (len(cigar_splitted) == 2):
#                         len_S = int(cigar_splitted[0][:-1])
#                         if len_S >= 10:
#                             clipping_pos = int(each_split[3])
#                             if ref_id not in clipping_pos_list_dict:
#                                 clipping_pos_list_dict[ref_id] = [clipping_pos]
#                             else:
#                                 clipping_pos_list_dict[ref_id].append(clipping_pos)
#
#     return clipping_pos_list_dict
#
#
# sam_in                  = '/Users/songweizhi/Desktop/mapHGT_biofilm/file_in/4D42.sam'
# clipping_pos_file_210   = '/Users/songweizhi/Desktop/mapHGT_biofilm/clipping_pos_210.txt'
# clipping_pos_file_D2    = '/Users/songweizhi/Desktop/mapHGT_biofilm/clipping_pos_D2.txt'
#
#
# clipping_pos_dict = get_clipping_pos_list(sam_in)
#
# clipping_pos_file_210_handle = open(clipping_pos_file_210, 'w')
# clipping_pos_file_D2_handle = open(clipping_pos_file_D2, 'w')
# for ref in clipping_pos_dict:
#     if ref == '2.10_chromosome':
#         clipping_pos_file_210_handle.write('\n'.join([str(i) for i in clipping_pos_dict[ref]]))
#     if ref == 'D2_c':
#         clipping_pos_file_D2_handle.write('\n'.join([str(i) for i in clipping_pos_dict[ref]]))
# clipping_pos_file_210_handle.close()
# clipping_pos_file_D2_handle.close()

