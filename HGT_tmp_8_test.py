import numpy as np
from Bio import SeqIO
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


# file in
ref_file = '/Users/songweizhi/Desktop/mapHGT_test/combined_refs.fna'
sam_file = '/Users/songweizhi/Desktop/mapHGT_test/TwoHGTs.sam'
#sam_file = '/Users/songweizhi/Desktop/000/TwoHGTs.sam'
sam_file = '/Users/songweizhi/Desktop/mapHGT_test/TwoHGTs.sam'


# mind line 45,


# file out
pwd_figure = '/Users/songweizhi/Desktop/mapHGT_test/TwoHGTs.png'


ref_length_dict = {}
ref_id_list = []
for seq_record in SeqIO.parse(ref_file, 'fasta'):
    ref_length_dict[seq_record.id] = len(seq_record.seq)
    ref_id_list.append(seq_record.id)


ref_id_list_sorted = sorted(ref_id_list)
ref_index_dict = {}
ref_index = 1
for ref in ref_id_list_sorted:
    ref_index_dict[ref] = ref_index
    ref_index += 1


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

        #print(read_id_base)
        read_ref_pos_dict[read_id] = ref_pos

        if read_id_base not in paired_reads_ref_dict:
            if read_strand == 'r1':
                paired_reads_ref_dict[read_id_base] = [ref_id, '']
            if read_strand == 'r2':
                paired_reads_ref_dict[read_id_base] = ['', ref_id]
        else:
            if read_strand == 'r1':
                paired_reads_ref_dict[read_id_base][0] = ref_id
            if read_strand == 'r2':
                paired_reads_ref_dict[read_id_base][1] = ref_id


ref_pos_dict = {}
paired_reads_link_loc_dict = {}
for paired_reads in paired_reads_ref_dict:
    paired_reads_ref_list = paired_reads_ref_dict[paired_reads]
    paired_read_1_ref     = paired_reads_ref_list[0]
    paired_read_1_ref_pos = read_ref_pos_dict['%s_r1' % paired_reads]
    paired_read_2_ref     = paired_reads_ref_list[1]
    paired_read_2_ref_pos = read_ref_pos_dict['%s_r2' % paired_reads]

    if paired_read_1_ref != paired_read_2_ref:
        print('%s\t%s\t%s\t%s\t%s' % (paired_reads, paired_read_1_ref, paired_read_1_ref_pos, paired_read_2_ref, paired_read_2_ref_pos))

        paired_reads_link_loc_dict[paired_reads] = [[paired_read_1_ref_pos, paired_read_2_ref_pos], [ref_index_dict[paired_read_1_ref], ref_index_dict[paired_read_2_ref]]]

        if paired_read_1_ref not in ref_pos_dict:
            ref_pos_dict[paired_read_1_ref] = [paired_read_1_ref_pos]
        else:
            ref_pos_dict[paired_read_1_ref].append(paired_read_1_ref_pos)

        if paired_read_2_ref not in ref_pos_dict:
            ref_pos_dict[paired_read_2_ref] = [paired_read_2_ref_pos]
        else:
            ref_pos_dict[paired_read_2_ref].append(paired_read_2_ref_pos)


# get plot
fig = plt.figure(1, figsize=(10, 5))
yticks_list = ['']
y_axis_pos = 1
for i in ref_id_list_sorted:

    # add point
    #plt.plot(ref_pos_dict[i], [y_axis_pos for num in ref_pos_dict[i]], marker='o', markersize=3, color='orange', markeredgewidth=0,linewidth=0)
    #plt.axis([0, ref_length_dict[i], 0, len(ref_pos_dict) + 1])

    # add line
    plt.plot([0, ref_length_dict[i]], [y_axis_pos, y_axis_pos], marker='|', markersize=5, color='black', linewidth=0.2, linestyle='dashed')

    # add links between paired reads
    for link in paired_reads_link_loc_dict:
        plt.plot(paired_reads_link_loc_dict[link][0], paired_reads_link_loc_dict[link][1], marker='o', markersize=0, color='skyblue', linewidth=0.1)

    # get yticks label
    yticks_list.append(i)

    y_axis_pos += 1

# add an extra line to avoid overlap between the top genome track and the plot outline
yticks_list.append('')

plt.yticks(np.arange(len(yticks_list)), yticks_list)
plt.tight_layout()
fig.savefig(pwd_figure, bbox_inches='tight', dpi=300)
plt.close()


print(paired_reads_link_loc_dict)




