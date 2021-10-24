import argparse
import numpy as np
from Bio import SeqIO
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def get_16S_pos_list(gbk_file, flk_length):

    ribosomal_RNA_16S_regions = []
    for seq_record in SeqIO.parse(gbk_file, 'genbank'):

        for gene_record in seq_record.features:
            if 'product' in gene_record.qualifiers:
                gene_product = gene_record.qualifiers['product'][0]
                if gene_product == '16S ribosomal RNA':
                    start_pos = gene_record.location.start - flk_length
                    end_pos = gene_record.location.end + flk_length

                    # add 16S pos to list
                    for i in range(start_pos, end_pos):
                        if i not in ribosomal_RNA_16S_regions:
                            ribosomal_RNA_16S_regions.append(i)

    return ribosomal_RNA_16S_regions


def cluster(data, maxgap):

    '''
    Arrange data into groups where successive elements differ by no more than *maxgap*

        >>> cluster([1, 6, 9, 100, 102, 105, 109, 134, 139], 5)
        [[1, 6, 9], [100, 102, 105, 109], [134, 139]]

        >>> cluster([1, 6, 9, 100, 102, 105, 109, 134, 139], 4)
        [[1], [6, 9], [100, 102, 105, 109], [134], [139]]
    '''

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
parser.add_argument('-flk16s',   required=False, type=int, default=2000, help='ignore paired reads in 16s regions')
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


gbk_file_210     = '/Users/songweizhi/Desktop/mapHGT_biofilm/file_in/210.gbk'
gbk_file_D2      = '/Users/songweizhi/Desktop/mapHGT_biofilm/file_in/D2.gbk'

plot_region_dict = {'2.10_chromosome': [3340000, 3350000], 'D2_c': [3930000, 3940000]}
if plot_hotspot is False:
    plot_region_dict = {'2.10_chromosome': [1, 4160809], 'D2_c': [1, 5020543]}


print('get 16S regions')
ribosomal_RNA_16S_pos_list_210 = get_16S_pos_list(gbk_file_210, flk_length_16S)
ribosomal_RNA_16S_pos_list_D2  = get_16S_pos_list(gbk_file_D2, flk_length_16S)
ribosomal_RNA_16S_pos_list_dict = {'2.10_chromosome': ribosomal_RNA_16S_pos_list_210, 'D2_c': ribosomal_RNA_16S_pos_list_D2}


ref_length_dict = {}
ref_id_list = []
for seq_record in SeqIO.parse(ref_file, 'fasta'):
    ref_length_dict[seq_record.id] = len(seq_record.seq)
    ref_id_list.append(seq_record.id)

print(ref_length_dict)

ref_id_list_sorted = sorted(ref_id_list)
ref_index_dict = {}
ref_index_dict_reverse = {}
ref_index = 1
for ref in ref_id_list_sorted:
    ref_index_dict[ref] = ref_index
    ref_index_dict_reverse[ref_index] = ref
    ref_index += 1


print('read in sam file')
paired_reads_ref_dict = {}
read_ref_pos_dict = {}
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

        read_ref_pos_dict[read_id] = ref_pos

        if read_id_base not in paired_reads_ref_dict:
            if read_strand == '1':
                paired_reads_ref_dict[read_id_base] = [ref_id, '']
            if read_strand == '2':
                paired_reads_ref_dict[read_id_base] = ['', ref_id]
        else:
            if read_strand == '1':
                paired_reads_ref_dict[read_id_base][0] = ref_id
            if read_strand == '2':
                paired_reads_ref_dict[read_id_base][1] = ref_id


print('store paired reads in dict')
ref_pos_dict = {}
paired_reads_link_loc_dict = {}
test_dict = {}
for paired_reads in paired_reads_ref_dict:
    paired_reads_ref_list = paired_reads_ref_dict[paired_reads]
    paired_read_1_ref     = paired_reads_ref_list[0]
    paired_read_2_ref     = paired_reads_ref_list[1]
    paired_read_1_ref_pos = read_ref_pos_dict['%s_1' % paired_reads]
    paired_read_2_ref_pos = read_ref_pos_dict['%s_2' % paired_reads]

    if paired_read_1_ref != paired_read_2_ref:

        test_dict[paired_reads] = {paired_read_1_ref:paired_read_1_ref_pos, paired_read_2_ref:paired_read_2_ref_pos}
        paired_reads_link_loc_dict[paired_reads] = [[paired_read_1_ref_pos, paired_read_2_ref_pos], [ref_index_dict[paired_read_1_ref], ref_index_dict[paired_read_2_ref]]]

        if paired_read_1_ref not in ref_pos_dict:
            ref_pos_dict[paired_read_1_ref] = [paired_read_1_ref_pos]
        else:
            ref_pos_dict[paired_read_1_ref].append(paired_read_1_ref_pos)

        if paired_read_2_ref not in ref_pos_dict:
            ref_pos_dict[paired_read_2_ref] = [paired_read_2_ref_pos]
        else:
            ref_pos_dict[paired_read_2_ref].append(paired_read_2_ref_pos)


# remove 16S ribosomal RNA regions from ref_pos_dict
print('remove 16S ribosomal RNA regions from ref_pos_dict')
ref_pos_dict_no_16S = {}
for ref in ref_pos_dict:
    ref_pos_list = ref_pos_dict[ref]
    ref_pos_list_no_16S = []
    for i in ref_pos_list:
        if i not in ribosomal_RNA_16S_pos_list_dict[ref]:
            ref_pos_list_no_16S.append(i)

    ref_pos_dict_no_16S[ref] = ref_pos_list_no_16S


if no16s is True:
    ref_pos_dict = ref_pos_dict_no_16S


# do clustering and write out clustered locations
print('do clustering')
ref_clustered_pos_dict = {}
hotspots_location_handle = open(hotspots_location_file, 'w')
for ref in ref_pos_dict:
    ref_pos_list = sorted(ref_pos_dict[ref])
    ref_pos_list_clustered = cluster(ref_pos_list, cluster_max_gap)
    for each_cluster in ref_pos_list_clustered:
        if len(each_cluster) >= min_cluster_reads_num:
            hotspots_location_handle.write('%s\t%s\n' % (ref, ','.join([str(i) for i in each_cluster])))
            for each_pos in each_cluster:
                if ref not in ref_clustered_pos_dict:
                    ref_clustered_pos_dict[ref] = [each_pos]
                else:
                    ref_clustered_pos_dict[ref].append(each_pos)
hotspots_location_handle.close()


paired_reads_link_loc_dict_clustered = {}
for each_test in test_dict:

    been_clustered_list = []
    for ref in ref_clustered_pos_dict:
        if test_dict[each_test][ref] in ref_clustered_pos_dict[ref]:
            been_clustered_list.append(True)
        else:
            been_clustered_list.append(False)
    if been_clustered_list == [True, True]:
        value_lol = [[], []]
        for i in test_dict[each_test]:
            value_lol[0].append(test_dict[each_test][i])
            value_lol[1].append(ref_index_dict[i])
        paired_reads_link_loc_dict_clustered[each_test] = value_lol

paired_reads_link_loc_dict = paired_reads_link_loc_dict_clustered


####################################### adjust pos to only plot specified regions ######################################

plot_region_len_dict = {}
for i in plot_region_dict:
    plot_region_len_dict[i] = plot_region_dict[i][1] - plot_region_dict[i][0]

ribosomal_RNA_16S_pos_list_dict_adjusted = {}
for genome in ribosomal_RNA_16S_pos_list_dict:
    genome_16s_pos_list = ribosomal_RNA_16S_pos_list_dict[genome]
    genome_16s_pos_list_updated = []
    for pos in genome_16s_pos_list:
        if plot_region_dict[genome][0] < pos < plot_region_dict[genome][1]:
            pos_updated = pos - plot_region_dict[genome][0]
            genome_16s_pos_list_updated.append(pos_updated)
    ribosomal_RNA_16S_pos_list_dict_adjusted[genome] = genome_16s_pos_list_updated

paired_reads_link_loc_dict_updated = {}
for each in paired_reads_link_loc_dict:
    if (plot_region_dict[ref_index_dict_reverse[paired_reads_link_loc_dict[each][1][0]]][0] < paired_reads_link_loc_dict[each][0][0] < plot_region_dict[ref_index_dict_reverse[paired_reads_link_loc_dict[each][1][0]]][1]) and (plot_region_dict[ref_index_dict_reverse[paired_reads_link_loc_dict[each][1][1]]][0] < paired_reads_link_loc_dict[each][0][1] < plot_region_dict[ref_index_dict_reverse[paired_reads_link_loc_dict[each][1][1]]][1]):
        paired_reads_link_loc_dict_updated[each] = paired_reads_link_loc_dict[each]
        paired_reads_link_loc_dict_updated[each] = [[paired_reads_link_loc_dict[each][0][0] - plot_region_dict[ref_index_dict_reverse[paired_reads_link_loc_dict[each][1][0]]][0], paired_reads_link_loc_dict[each][0][1] - plot_region_dict[ref_index_dict_reverse[paired_reads_link_loc_dict[each][1][1]]][0]], paired_reads_link_loc_dict[each][1]]


########################################################################################################################


print('Plotting')
# get plot
fig = plt.figure(1, figsize=(100, 5))
yticks_list = ['']
y_axis_pos = 1
for i in ref_id_list_sorted:

    # add point
    plt.plot(ribosomal_RNA_16S_pos_list_dict_adjusted[i], [y_axis_pos for num in ribosomal_RNA_16S_pos_list_dict_adjusted[i]], marker='|', markersize=2, color='orange', markeredgewidth=1,linewidth=0)
    plt.axis([0, plot_region_len_dict[i], 0, len(ribosomal_RNA_16S_pos_list_dict_adjusted) + 1])

    # add line
    plt.plot([0, plot_region_len_dict[i]], [y_axis_pos, y_axis_pos], marker='|', markersize=5, color='black', linewidth=0.2, linestyle='dashed')

    # get yticks label
    yticks_list.append(i)

    y_axis_pos += 1


# add links between paired reads
for link in paired_reads_link_loc_dict_updated:
    plt.plot(paired_reads_link_loc_dict_updated[link][0], paired_reads_link_loc_dict_updated[link][1], marker='o', markersize=0, color='skyblue', linewidth=0.1)


# add an extra line to avoid overlap between the top genome track and the plot outline
yticks_list.append('')

plt.yticks(np.arange(len(yticks_list)), yticks_list)
plt.tight_layout()
#fig.savefig(pwd_figure, bbox_inches='tight', dpi=300)
fig.savefig(pwd_figure)
plt.close()

print('Done!')
