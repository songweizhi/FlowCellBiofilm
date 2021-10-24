import os
import argparse
from Bio import SeqIO
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def cigar_splitter(cigar):

    # get the position of letters
    letter_pos_list = []
    n = 0
    for each_element in cigar:
        if each_element.isalpha() is True:
            letter_pos_list.append(n)
        n += 1

    # split cigar
    index = 0
    cigar_splitted = []
    while index <= len(letter_pos_list) - 1:
        if index == 0:
            cigar_splitted.append(cigar[:(letter_pos_list[index] + 1)])
        else:
            cigar_splitted.append(cigar[(letter_pos_list[index - 1] + 1):(letter_pos_list[index] + 1)])
        index += 1

    return cigar_splitted


def get_ref_mean_depth(sam_depth_file):
    ref_id_list = []
    ref_len_dict = {}
    ref_depth_dict = {}
    for each_line in open(sam_depth_file):

        each_line_split = each_line.strip().split('\t')
        ref_id = each_line_split[0]
        pos = int(each_line_split[1])
        pos_depth = int(each_line_split[2])

        # get ref_id_list
        if ref_id not in ref_id_list:
            ref_id_list.append(ref_id)

        # get ref_len_dict
        if ref_id not in ref_len_dict:
            ref_len_dict[ref_id] = 1
        else:
            ref_len_dict[ref_id] += 1

        # get ref_depth_dict
        if ref_id not in ref_depth_dict:
            ref_depth_dict[ref_id] = pos_depth
        else:
            ref_depth_dict[ref_id] += pos_depth

    mean_depth_dict = {}
    for ref in sorted(ref_id_list):
        ref_len = ref_len_dict[ref]
        ref_depth = ref_depth_dict[ref]
        ref_depth_mean = float("{0:.2f}".format(ref_depth / ref_len))
        mean_depth_dict[ref] = ref_depth_mean

    return mean_depth_dict


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


parser = argparse.ArgumentParser()

parser.add_argument('-r',        required=True,                          help='reference')
parser.add_argument('-s',        required=True,                          help='sam file')
parser.add_argument('-pg',       required=True, type=int, default=5000,  help='paired reads maximum gap')
parser.add_argument('-pn',       required=True, type=int, default=5,     help='minimum number of paired reads in a cluster')
parser.add_argument('-cg',       required=True, type=int, default=3,     help='clipping reads maximum gap')
parser.add_argument('-cn',       required=True, type=int, default=3,     help='minimum number of clipping reads in a cluster')
parser.add_argument('-p',        required=True,                          help='output prefix')
parser.add_argument('-flk16s',   required=False, type=int, default=5000, help='ignore paired reads in 16s regions')
parser.add_argument('-no16s',    required=False, action='store_true',    help='ignore 16s regions')
parser.add_argument('-skip_sam', required=False, action='store_true',    help='skip parsing sam file')
parser.add_argument('-gbk',      required=True,                          help='gbk folder, 210.gbk and D2.gbk')


args = vars(parser.parse_args())
ref_file                       = args['r']
sam_file                       = args['s']
output_prefix                  = args['p']
paired_cluster_max_gap         = args['pg']
paired_min_cluster_reads_num   = args['pn']
clipping_cluster_max_gap       = args['cg']
clipping_min_cluster_reads_num = args['cn']
no16s                          = args['no16s']
flk_length_16S                 = args['flk16s']
skip_parse_sam                 = args['skip_sam']
gbk_folder                     = args['gbk']

pwd_figure_paired                   = '%s_plot_paired_pg%s_pn%s_cg%s_cn%s.pdf'              % (output_prefix, paired_cluster_max_gap, paired_min_cluster_reads_num, clipping_cluster_max_gap, clipping_min_cluster_reads_num)
pwd_figure_clipping                 = '%s_plot_clipping_pg%s_pn%s_cg%s_cn%s.pdf'            % (output_prefix, paired_cluster_max_gap, paired_min_cluster_reads_num, clipping_cluster_max_gap, clipping_min_cluster_reads_num)
hotspots_location_file              = '%s_hotspots_paired_pg%s_pn%s_cg%s_cn%s.txt'          % (output_prefix, paired_cluster_max_gap, paired_min_cluster_reads_num, clipping_cluster_max_gap, clipping_min_cluster_reads_num)
if no16s is True:
    pwd_figure_paired               = '%s_plot_paired_pg%s_pn%s_cg%s_cn%s_no16s.pdf'        % (output_prefix, paired_cluster_max_gap, paired_min_cluster_reads_num, clipping_cluster_max_gap, clipping_min_cluster_reads_num)
    pwd_figure_clipping             = '%s_plot_clipping_pg%s_pn%s_cg%s_cn%s_no16s.pdf'      % (output_prefix, paired_cluster_max_gap, paired_min_cluster_reads_num, clipping_cluster_max_gap, clipping_min_cluster_reads_num)
    hotspots_location_file          = '%s_hotspots_paired_pg%s_pn%s_cg%s_cn%s_no16s.txt'    % (output_prefix, paired_cluster_max_gap, paired_min_cluster_reads_num, clipping_cluster_max_gap, clipping_min_cluster_reads_num)

gbk_file_210                        = '%s/210.gbk' % gbk_folder
gbk_file_D2                         = '%s/D2.gbk'  % gbk_folder
min_cigar_M                         = 30
min_cigar_S                         = 30
paired_reads_ref_file               = '%s_paired_reads_ref.txt' % output_prefix


# get 16S regions
print('get 16S regions')
ribosomal_RNA_16S_pos_list_210  = get_16S_pos_list(gbk_file_210, flk_length_16S)
ribosomal_RNA_16S_pos_list_D2   = get_16S_pos_list(gbk_file_D2, flk_length_16S)
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
    ref_pos_list_clustered = cluster(sorted(ref_pos_dict[ref]), paired_cluster_max_gap, paired_min_cluster_reads_num)

    # only keep the clustered pos
    for each_cluster in ref_pos_list_clustered:
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
########################################################################################################################
########################################################################################################################

# tmp files
clipping_reads_matched_part                         = '%s_hotspots_clipping_reads_matched_part.txt'                      % (output_prefix)
clipping_reads_not_matched_part_seq                 = '%s_hotspots_clipping_reads_not_matched_part.fasta'                % (output_prefix)
clipping_reads_not_matched_part_seq_blast           = '%s_hotspots_clipping_reads_not_matched_part_blast.tab'            % (output_prefix)
clipping_reads_not_matched_part_seq_blast_BestHit   = '%s_hotspots_clipping_reads_not_matched_part_blast_BestHit.tab'    % (output_prefix)
clipping_reads_not_matched_part                     = '%s_hotspots_clipping_reads_not_matched_part.txt'                  % (output_prefix)
clipping_reads_matched_info_combined                = '%s_hotspots_clipping_matched_info_combined.txt'                   % (output_prefix)

clipping_reads_matched_to_diff_refs                 = '%s_hotspots_clipping_pg%s_pn%s_cg%s_cn%s.txt'                     % (output_prefix, paired_cluster_max_gap, paired_min_cluster_reads_num, clipping_cluster_max_gap, clipping_min_cluster_reads_num)
if no16s is True:
    clipping_reads_matched_to_diff_refs             = '%s_hotspots_clipping_pg%s_pn%s_cg%s_cn%s_no16s.txt'               % (output_prefix, paired_cluster_max_gap, paired_min_cluster_reads_num, clipping_cluster_max_gap, clipping_min_cluster_reads_num)


clipping_reads_matched_part_handle = open(clipping_reads_matched_part, 'w')
clipping_reads_not_matched_part_seq_handle = open(clipping_reads_not_matched_part_seq, 'w')
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
        cigar_splitted  = cigar_splitter(cigar)

        if ('S' in cigar) and (len(cigar_splitted) == 2):

            cigar_M_len = 0
            cigar_S_len = 0
            split_pos = 0
            if cigar_splitted[0][-1] == 'M':
                cigar_M_len = int(cigar_splitted[0][:-1])
                cigar_S_len = int(cigar_splitted[1][:-1])
                split_pos = ref_pos + cigar_M_len
            if cigar_splitted[1][-1] == 'M':
                cigar_M_len = int(cigar_splitted[1][:-1])
                cigar_S_len = int(cigar_splitted[0][:-1])
                split_pos = ref_pos

            if (cigar_M_len >= min_cigar_M) and (cigar_S_len >= min_cigar_S):

                read_seq_left = read_seq[: int(cigar_splitted[0][:-1])]
                read_seq_right = read_seq[-int(cigar_splitted[1][:-1]):]

                if cigar_splitted[0][-1] == 'M':
                    clipping_reads_matched_part_handle.write('%s_l\t%s\t%s\n' % (read_id, ref_id, split_pos))
                    clipping_reads_not_matched_part_seq_handle.write('>%s_r\n' % read_id)
                    clipping_reads_not_matched_part_seq_handle.write(read_seq_right + '\n')
                if cigar_splitted[1][-1] == 'M':
                    clipping_reads_matched_part_handle.write('%s_r\t%s\t%s\n' % (read_id, ref_id, split_pos))
                    clipping_reads_not_matched_part_seq_handle.write('>%s_l\n' % read_id)
                    clipping_reads_not_matched_part_seq_handle.write(read_seq_left + '\n')

clipping_reads_matched_part_handle.close()
clipping_reads_not_matched_part_seq_handle.close()


#makeblastdb_cmd = 'makeblastdb -in %s -dbtype nucl -parse_seqids -logfile /dev/null' % combined_refs
#os.system(makeblastdb_cmd)

blast_parameters = '-evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn -num_threads %s' % 1
blastn_cmd = 'blastn -query %s -db %s -out %s %s' % (clipping_reads_not_matched_part_seq, ref_file, clipping_reads_not_matched_part_seq_blast, blast_parameters)
os.system(blastn_cmd)

keep_BestHit_cmd = 'BioSAK BestHit -i %s -o %s' % (clipping_reads_not_matched_part_seq_blast, clipping_reads_not_matched_part_seq_blast_BestHit)
os.system(keep_BestHit_cmd)


clipping_reads_not_matched_part_handle = open(clipping_reads_not_matched_part, 'w')
for match in open(clipping_reads_not_matched_part_seq_blast_BestHit):

    match_split = match.strip().split('\t')
    query       = match_split[0]
    subject     = match_split[1]
    iden        = float(match_split[2])
    align_len   = int(match_split[3])
    sstart      = int(match_split[8])
    send        = int(match_split[9])
    query_len   = int(match_split[12])

    if (iden == 100) and (align_len == query_len):
        split_pos = 0
        if query[-1] == 'l':
            split_pos = send
        if query[-1] == 'r':
            split_pos = sstart

        clipping_reads_not_matched_part_handle.write('%s\t%s\t%s\n' % (query, subject, split_pos))

clipping_reads_not_matched_part_handle.close()


os.system('cat %s %s > %s' % (clipping_reads_matched_part, clipping_reads_not_matched_part, clipping_reads_matched_info_combined))


clipping_reads_match_info_dict = {}
for each_read_segment in open(clipping_reads_matched_info_combined):
    each_read_segment_split = each_read_segment.strip().split('\t')

    read_id = each_read_segment_split[0][:-2]
    matched_ref = each_read_segment_split[1]
    matched_ref_pos = int(each_read_segment_split[2])

    if each_read_segment_split[0][-1] == 'l':
        if read_id not in clipping_reads_match_info_dict:
            clipping_reads_match_info_dict[read_id] = [matched_ref, matched_ref_pos, '', '']
        else:
            clipping_reads_match_info_dict[read_id][0] = matched_ref
            clipping_reads_match_info_dict[read_id][1] = matched_ref_pos

    if each_read_segment_split[0][-1] == 'r':
        if read_id not in clipping_reads_match_info_dict:
            clipping_reads_match_info_dict[read_id] = ['', '', matched_ref, matched_ref_pos]
        else:
            clipping_reads_match_info_dict[read_id][2] = matched_ref
            clipping_reads_match_info_dict[read_id][3] = matched_ref_pos


clipping_reads_matched_to_diff_refs_handle = open(clipping_reads_matched_to_diff_refs, 'w')
for i in clipping_reads_match_info_dict:
    i_value = clipping_reads_match_info_dict[i]
    if (i_value[0] != '') and (i_value[2] != '') and (i_value[0] != i_value[2]):
        clipping_reads_matched_to_diff_refs_handle.write('%s\t%s\n' % (i, '\t'.join([str(j) for j in i_value])))
clipping_reads_matched_to_diff_refs_handle.close()


clipping_reads_pos_210 = []
clipping_reads_pos_D2 = []
for each_clipping_read in open(clipping_reads_matched_to_diff_refs):
    each_clipping_read_split = each_clipping_read.strip().split('\t')

    if (each_clipping_read_split[1] == '2.10_chromosome') and (each_clipping_read_split[3] == 'D2_c'):
        clipping_reads_pos_210.append(int(each_clipping_read_split[2]))
        clipping_reads_pos_D2.append(int(each_clipping_read_split[4]))

    if (each_clipping_read_split[1] == 'D2_c') and (each_clipping_read_split[3] == '2.10_chromosome'):
        clipping_reads_pos_210.append(int(each_clipping_read_split[4]))
        clipping_reads_pos_D2.append(int(each_clipping_read_split[2]))


clipping_reads_pos_210_sorted_clustered = cluster(sorted(clipping_reads_pos_210), clipping_cluster_max_gap, clipping_min_cluster_reads_num)
clipping_reads_pos_D2_sorted_clustered  = cluster(sorted(clipping_reads_pos_D2), clipping_cluster_max_gap, clipping_min_cluster_reads_num)


clipping_reads_pos_210_kept_pos_set = set()
for each_cluster_210 in clipping_reads_pos_210_sorted_clustered:
    for each_cluster_210_element in each_cluster_210:
        clipping_reads_pos_210_kept_pos_set.add(each_cluster_210_element)

clipping_reads_pos_D2_kept_pos_set = set()
for each_cluster_D2 in clipping_reads_pos_D2_sorted_clustered:
    for each_cluster_D2_element in each_cluster_D2:
        clipping_reads_pos_D2_kept_pos_set.add(each_cluster_D2_element)


clipping_reads_pos_210_filtered = []
clipping_reads_pos_D2_filtered = []
for (x, y) in zip(clipping_reads_pos_210, clipping_reads_pos_D2):
    if (x in clipping_reads_pos_210_kept_pos_set) and (y in clipping_reads_pos_D2_kept_pos_set):
        clipping_reads_pos_210_filtered.append(x)
        clipping_reads_pos_D2_filtered.append(y)


# remove 16s regions
print('remove 16s regions')
clipping_reads_pos_210_filtered_no_16s = []
clipping_reads_pos_D2_filtered_no_16s = []
for (m, n) in zip(clipping_reads_pos_210_filtered, clipping_reads_pos_D2_filtered):
    if (m not in ribosomal_RNA_16S_pos_list_dict['2.10_chromosome']) and (n not in ribosomal_RNA_16S_pos_list_dict['D2_c']):
        clipping_reads_pos_210_filtered_no_16s.append(m)
        clipping_reads_pos_D2_filtered_no_16s.append(n)


########################################################################################################################
########################################################################################################################
########################################################################################################################

print('Plotting')

fig = plt.figure(figsize=(10, 10))
if no16s is True:
    plt.scatter(plot_num_list_dict_no_16s['2.10_chromosome'], plot_num_list_dict_no_16s['D2_c'], s=12, linewidths=0, marker='.', c='blue', alpha=0.3)
else:
    plt.scatter(plot_num_list_dict['2.10_chromosome'], plot_num_list_dict['D2_c'], s=18, linewidths=0, marker='.', c='blue', alpha=0.3)

# set axis range
plt.xlim(0, ref_length_dict['2.10_chromosome'])
plt.ylim(0, ref_length_dict['D2_c'])

# add axis label
plt.xlabel('2.10 chromosome (bp)')
plt.ylabel('D2 chromosome (bp)')

plt.tight_layout()
plt.savefig(pwd_figure_paired)
plt.close()


fig2 = plt.figure(figsize=(10, 10))
if no16s is True:
    plt.scatter(clipping_reads_pos_210_filtered_no_16s, clipping_reads_pos_D2_filtered_no_16s, s=12, linewidths=0, marker='.', c='red', alpha=0.3)
else:
    plt.scatter(clipping_reads_pos_210_filtered, clipping_reads_pos_D2_filtered, s=18, linewidths=0, marker='.', c='red', alpha=0.3)

# set axis range
plt.xlim(0, ref_length_dict['2.10_chromosome'])
plt.ylim(0, ref_length_dict['D2_c'])

# add axis label
plt.xlabel('2.10 chromosome (bp)')
plt.ylabel('D2 chromosome (bp)')

plt.tight_layout()
plt.savefig(pwd_figure_clipping)
plt.close()

print('Done!')


# remove tmp files
os.remove(clipping_reads_matched_part)
os.remove(clipping_reads_not_matched_part_seq)
os.remove(clipping_reads_not_matched_part_seq_blast)
os.remove(clipping_reads_not_matched_part_seq_blast_BestHit)
os.remove(clipping_reads_not_matched_part)
os.remove(clipping_reads_matched_info_combined)
# os.remove(clipping_reads_matched_to_diff_refs)
# os.remove(hotspots_location_file)


''''

cd /Users/songweizhi/Desktop/111
python3 ~/PycharmProjects/FlowCellBiofilm/MapHGT_2020_08_03.py -r 0_ref_wt/combined_refs.fasta -s OneHGT_50_to_50.sam -p OneHGT_50_to_50 -pg 5000 -pn 5 -cg 3 -cn 3 -no16s -skip_sam -gbk /Users/songweizhi/Desktop/mapHGT_biofilm/file_in
python3 ~/PycharmProjects/FlowCellBiofilm/MapHGT_2020_08_03.py -r 0_ref_wt/combined_refs.fasta -s OneHGT_75_to_25.sam -p OneHGT_75_to_25 -pg 5000 -pn 5 -cg 3 -cn 3 -no16s -skip_sam -gbk /Users/songweizhi/Desktop/mapHGT_biofilm/file_in
python3 ~/PycharmProjects/FlowCellBiofilm/MapHGT_2020_08_03.py -r 0_ref_wt/combined_refs.fasta -s OneHGT_90_to_10.sam -p OneHGT_90_to_10 -pg 5000 -pn 5 -cg 3 -cn 3 -no16s -skip_sam -gbk /Users/songweizhi/Desktop/mapHGT_biofilm/file_in
python3 ~/PycharmProjects/FlowCellBiofilm/MapHGT_2020_08_03.py -r 0_ref_wt/combined_refs.fasta -s OneHGT_95_to_05.sam -p OneHGT_95_to_05 -pg 5000 -pn 5 -cg 3 -cn 3 -no16s -skip_sam -gbk /Users/songweizhi/Desktop/mapHGT_biofilm/file_in
python3 ~/PycharmProjects/FlowCellBiofilm/MapHGT_2020_08_03.py -r 0_ref_wt/combined_refs.fasta -s OneHGT_99_to_01.sam -p OneHGT_99_to_01 -pg 5000 -pn 5 -cg 3 -cn 3 -no16s -skip_sam -gbk /Users/songweizhi/Desktop/mapHGT_biofilm/file_in

python3 ~/PycharmProjects/FlowCellBiofilm/MapHGT_2020_08_03.py -r 0_ref_wt/combined_refs.fasta -s 4D42.sam -p 4D42 -pg 5000 -pn 5 -cg 3 -cn 1 -no16s -skip_sam -gbk /Users/songweizhi/Desktop/mapHGT_biofilm/file_in
python3 ~/PycharmProjects/FlowCellBiofilm/MapHGT_2020_08_03.py -r 0_ref_wt/combined_refs.fasta -s 4D42.sam -p 4D42 -pg 5000 -pn 5 -cg 3 -cn 1 -no16s -skip_sam -gbk /Users/songweizhi/Desktop/mapHGT_biofilm/file_in
python3 ~/PycharmProjects/FlowCellBiofilm/MapHGT_2020_08_03.py -r 0_ref_wt/combined_refs.fasta -s 8D42.sam -p 8D42 -pg 5000 -pn 5 -cg 3 -cn 1 -no16s -skip_sam -gbk /Users/songweizhi/Desktop/mapHGT_biofilm/file_in
python3 ~/PycharmProjects/FlowCellBiofilm/MapHGT_2020_08_03.py -r 0_ref_wt/combined_refs.fasta -s 8D42.sam -p 8D42 -pg 5000 -pn 5 -cg 3 -cn 1 -no16s -skip_sam -gbk /Users/songweizhi/Desktop/mapHGT_biofilm/file_in
python3 ~/PycharmProjects/FlowCellBiofilm/MapHGT_2020_08_03.py -r 0_ref_wt/combined_refs.fasta -s 12D42.sam -p 12D42 -pg 5000 -pn 5 -cg 3 -cn 1 -no16s -skip_sam -gbk /Users/songweizhi/Desktop/mapHGT_biofilm/file_in
python3 ~/PycharmProjects/FlowCellBiofilm/MapHGT_2020_08_03.py -r 0_ref_wt/combined_refs.fasta -s 12D42.sam -p 12D42 -pg 5000 -pn 5 -cg 3 -cn 1 -no16s -skip_sam -gbk /Users/songweizhi/Desktop/mapHGT_biofilm/file_in

50_to_50：{'2.10_chromosome': 55.29, 'D2_c': 45.67}
75_to_25：{'2.10_chromosome': 82.91, 'D2_c': 22.84}
90_to_10：{'2.10_chromosome': 99.47, 'D2_c': 9.14}
95_to_05: {'2.10_chromosome': 105.0, 'D2_c': 4.62}
99_to_01: {'2.10_chromosome': 109.41, 'D2_c': 1.53}


module load python/3.7.3
source ~/mypython3env/bin/activate
cd /srv/scratch/z5039045/Flow_cell_biofilm/3_mapping_mapHGT

python3 MapHGT_2020_08_03.py -r ../0_References_for_mapHGT/combined_refs.fasta -s 4D9.sam -p 4D9 -pg 5000 -pn 5 -cg 3 -cn 3 -no16s -gbk gbk_files
python3 MapHGT_2020_08_03.py -r ../0_References_for_mapHGT/combined_refs.fasta -s 4D9.sam -p 4D9 -pg 5000 -pn 5 -cg 3 -cn 1 -no16s -gbk gbk_files -skip_sam
python3 MapHGT_2020_08_03.py -r ../0_References_for_mapHGT/combined_refs.fasta -s 4D18.sam -p 4D18 -pg 5000 -pn 5 -cg 3 -cn 3 -no16s -gbk gbk_files
python3 MapHGT_2020_08_03.py -r ../0_References_for_mapHGT/combined_refs.fasta -s 4D18.sam -p 4D18 -pg 5000 -pn 5 -cg 3 -cn 1 -no16s -gbk gbk_files -skip_sam
python3 MapHGT_2020_08_03.py -r ../0_References_for_mapHGT/combined_refs.fasta -s 4D27.sam -p 4D27 -pg 5000 -pn 5 -cg 3 -cn 3 -no16s -gbk gbk_files
python3 MapHGT_2020_08_03.py -r ../0_References_for_mapHGT/combined_refs.fasta -s 4D27.sam -p 4D27 -pg 5000 -pn 5 -cg 3 -cn 1 -no16s -gbk gbk_files -skip_sam
python3 MapHGT_2020_08_03.py -r ../0_References_for_mapHGT/combined_refs.fasta -s 4D42.sam -p 4D42 -pg 5000 -pn 5 -cg 3 -cn 3 -no16s -gbk gbk_files
python3 MapHGT_2020_08_03.py -r ../0_References_for_mapHGT/combined_refs.fasta -s 4D42.sam -p 4D42 -pg 5000 -pn 5 -cg 3 -cn 1 -no16s -gbk gbk_files -skip_sam

python3 MapHGT_2020_08_03.py -r ../0_References_for_mapHGT/combined_refs.fasta -s 8D9.sam -p 8D9 -pg 5000 -pn 5 -cg 3 -cn 3 -no16s -gbk gbk_files
python3 MapHGT_2020_08_03.py -r ../0_References_for_mapHGT/combined_refs.fasta -s 8D9.sam -p 8D9 -pg 5000 -pn 5 -cg 3 -cn 1 -no16s -gbk gbk_files -skip_sam
python3 MapHGT_2020_08_03.py -r ../0_References_for_mapHGT/combined_refs.fasta -s 8D18.sam -p 8D18 -pg 5000 -pn 5 -cg 3 -cn 3 -no16s -gbk gbk_files
python3 MapHGT_2020_08_03.py -r ../0_References_for_mapHGT/combined_refs.fasta -s 8D18.sam -p 8D18 -pg 5000 -pn 5 -cg 3 -cn 1 -no16s -gbk gbk_files -skip_sam
python3 MapHGT_2020_08_03.py -r ../0_References_for_mapHGT/combined_refs.fasta -s 8D27.sam -p 8D27 -pg 5000 -pn 5 -cg 3 -cn 3 -no16s -gbk gbk_files
python3 MapHGT_2020_08_03.py -r ../0_References_for_mapHGT/combined_refs.fasta -s 8D27.sam -p 8D27 -pg 5000 -pn 5 -cg 3 -cn 1 -no16s -gbk gbk_files -skip_sam
python3 MapHGT_2020_08_03.py -r ../0_References_for_mapHGT/combined_refs.fasta -s 8D42.sam -p 8D42 -pg 5000 -pn 5 -cg 3 -cn 3 -no16s -gbk gbk_files
python3 MapHGT_2020_08_03.py -r ../0_References_for_mapHGT/combined_refs.fasta -s 8D42.sam -p 8D42 -pg 5000 -pn 5 -cg 3 -cn 1 -no16s -gbk gbk_files -skip_sam

python3 MapHGT_2020_08_03.py -r ../0_References_for_mapHGT/combined_refs.fasta -s 12D9.sam -p 12D9 -pg 5000 -pn 5 -cg 3 -cn 3 -no16s -gbk gbk_files
python3 MapHGT_2020_08_03.py -r ../0_References_for_mapHGT/combined_refs.fasta -s 12D9.sam -p 12D9 -pg 5000 -pn 5 -cg 3 -cn 1 -no16s -gbk gbk_files -skip_sam
python3 MapHGT_2020_08_03.py -r ../0_References_for_mapHGT/combined_refs.fasta -s 12D18.sam -p 12D18 -pg 5000 -pn 5 -cg 3 -cn 3 -no16s -gbk gbk_files
python3 MapHGT_2020_08_03.py -r ../0_References_for_mapHGT/combined_refs.fasta -s 12D18.sam -p 12D18 -pg 5000 -pn 5 -cg 3 -cn 1 -no16s -gbk gbk_files -skip_sam
python3 MapHGT_2020_08_03.py -r ../0_References_for_mapHGT/combined_refs.fasta -s 12D27.sam -p 12D27 -pg 5000 -pn 5 -cg 3 -cn 3 -no16s -gbk gbk_files
python3 MapHGT_2020_08_03.py -r ../0_References_for_mapHGT/combined_refs.fasta -s 12D27.sam -p 12D27 -pg 5000 -pn 5 -cg 3 -cn 1 -no16s -gbk gbk_files -skip_sam
python3 MapHGT_2020_08_03.py -r ../0_References_for_mapHGT/combined_refs.fasta -s 12D42.sam -p 12D42 -pg 5000 -pn 5 -cg 3 -cn 3 -no16s -gbk gbk_files
python3 MapHGT_2020_08_03.py -r ../0_References_for_mapHGT/combined_refs.fasta -s 12D42.sam -p 12D42 -pg 5000 -pn 5 -cg 3 -cn 1 -no16s -gbk gbk_files -skip_sam

'''
