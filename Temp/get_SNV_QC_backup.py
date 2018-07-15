import os
import glob
import shutil
import argparse
from Bio import SeqIO
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def take_kmer_mean(num_list, k_mer):

    k_mer_average_list = []
    n = 0

    while n < len(num_list):

        if n <= len(num_list) - k_mer:

            to_add_value = list(range(1, k_mer + 1))

            pos_list = []
            for each in to_add_value:
                pos_list.append(num_list[n + each - 1])

            # get average
            current_average = sum(pos_list) / k_mer
            k_mer_average_list.append(float("{0:.2f}".format(current_average)))

            n += 1

        elif (len(num_list) - k_mer) < n < len(num_list):
            k_mer_average_list.append(num_list[n])

            n += 1

    return k_mer_average_list


def plot_sam_depth(depth_file, seq_to_plot, start_pos, end_pos, each_snv_pos, mean_depth_start, mean_depth_end, k_mer, bps_to_marker, plot_filename):
    #print('Extracting absolute depth from input file')
    x = []
    y = []
    bp_num = 1
    current_pos = 0
    total_depth_left_side = 0
    total_depth_right_side = 0

    for each_base in open(depth_file):
        each_base_split = each_base.strip().split('\t')
        seq_id = each_base_split[0]
        pos = int(each_base_split[1])
        depth = int(each_base_split[2])

        if seq_id == seq_to_plot:

            # get total_depth_left_side
            if mean_depth_start <= pos < each_snv_pos:
                total_depth_left_side += depth

            # get total_depth_right_side
            if each_snv_pos < pos <= mean_depth_end:
                total_depth_right_side += depth

            if pos < start_pos:
                pass

            elif pos == start_pos:
                x.append(pos)
                y.append(depth)
                current_pos = pos

            else:
                start_0 = None
                end_0 = None
                to_add = []

                if (pos == current_pos + 1) and (pos <= end_pos):
                    x.append(pos)
                    y.append(depth)
                    current_pos = pos

                elif (pos > current_pos + 1) and (pos <= end_pos):

                    # add zero
                    start_0 = current_pos + 1
                    end_0 = pos - 1
                    to_add = list(range(start_0, end_0 + 1))
                    for each_0 in to_add:
                        x.append(each_0)
                        y.append(0)

                    x.append(pos)
                    y.append(depth)
                    current_pos = pos

                elif (pos > current_pos + 1) and (pos > end_pos):

                    # add zero
                    start_0 = current_pos + 1
                    end_0 = end_pos
                    to_add = list(range(start_0, end_0 + 1))
                    for each_0 in to_add:
                        x.append(each_0)
                        y.append(0)

                    current_pos = pos

    #print('Calculating k-mer means')
    y = take_kmer_mean(y, k_mer)

    # Change the color and its transparency
    plt.plot(x, y, color="skyblue", alpha=0.7, linewidth=0.7)

    # titles
    plt.title(plot_filename,fontsize=7)
    plt.xlabel('Position (bp)', fontsize=10)
    plt.ylabel('Depth (X)', fontsize=10)

    plt.xticks(fontsize=7)
    plt.yticks(fontsize=7)

    ymax = 0
    if max(y) <= 5000:
        ymax = 5000
    elif (5000 < max(y) <= 10000):
        ymax = 10000
    else:
        ymax = max(y)

    plt.ylim(ymin=0, ymax=ymax)

    # add lines to specified positions
    if bps_to_marker != None:
        bps_to_marker_list = bps_to_marker.strip().split(',')
        for each_line in bps_to_marker_list:
            plt.axvline(x=int(each_line), c='red', linewidth=0.3)

    # Get plot
    plt.savefig('%s.png' % plot_filename, dpi=300)
    plt.close()

    return total_depth_left_side, total_depth_right_side


def combine_continuous_deletion(file_in, file_out):

    file_out_handle = open(file_out, 'w')
    current_seq = ''
    current_pos_start = 0
    current_profile_start = ''
    current_pos = 0
    deleted_ncs = ''

    for each_snv in open(file_in):

        if each_snv.startswith('\t'):
            file_out_handle.write(each_snv)

        else:
            each_snv_seq = each_snv.strip().split('\t')[0].split('|')[0]
            each_snv_pos = int(each_snv.strip().split('\t')[0].split('|')[1])
            each_snv_pos_wt = each_snv.strip().split('\t')[0].split('|')[2]
            each_snv_pos_v = each_snv.strip().split('\t')[0].split('|')[3]

            if (each_snv_pos_v == '-') and (current_seq == '') and (current_pos == 0) and (deleted_ncs == ''):
                current_seq = each_snv_seq
                current_pos_start = each_snv_pos
                current_profile_start = '\t'.join(each_snv.strip().split('\t')[1:])
                current_pos = each_snv_pos
                deleted_ncs = each_snv_pos_wt
            elif (each_snv_seq == current_seq) and (each_snv_pos == (current_pos + 1)) and (each_snv_pos_v != '-'):
                for_out = '%s|%s-%s|%s|-\t%s\n' % (
                current_seq, current_pos_start, current_pos, deleted_ncs, current_profile_start)
                file_out_handle.write(for_out)
                file_out_handle.write(each_snv)
                current_seq = ''
                current_pos_start = 0
                current_profile_start = ''
                current_pos = 0
                deleted_ncs = ''

            elif (each_snv_seq == current_seq) and (each_snv_pos == (current_pos + 1)) and (each_snv_pos_v == '-'):
                current_pos = each_snv_pos
                deleted_ncs += each_snv_pos_wt
            elif each_snv_pos != (current_pos + 1):
                if (current_seq != '') and (current_pos != 0) and (deleted_ncs != ''):
                    for_out = ''
                    if len(deleted_ncs) == 1:
                        for_out = '%s|%s|%s|-\t%s\n' % (current_seq, current_pos, deleted_ncs, current_profile_start)
                    elif len(deleted_ncs) > 1:
                        for_out = '%s|%s-%s|%s|-\t%s\n' % (
                        current_seq, current_pos_start, current_pos, deleted_ncs, current_profile_start)
                    file_out_handle.write(for_out)

                    if each_snv_pos_v == '-':
                        current_seq = each_snv_seq
                        current_pos_start = each_snv_pos
                        current_profile_start = '\t'.join(each_snv.strip().split('\t')[1:])
                        current_pos = each_snv_pos
                        deleted_ncs = each_snv_pos_wt
                    else:
                        file_out_handle.write(each_snv)
                        current_seq = ''
                        current_pos_start = 0
                        current_profile_start = ''
                        current_pos = 0
                        deleted_ncs = ''

                elif (current_seq == '') and (current_pos == 0) and (deleted_ncs == '') and (each_snv_pos_v != '-'):
                    file_out_handle.write(each_snv)

    file_out_handle.close()


########################################## specify input files and parameters ##########################################

parser = argparse.ArgumentParser(description='', add_help=False)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')

optional.add_argument('-h', action='help', help='Show this help message and exit')
required.add_argument('-snv', dest='SNV', nargs='?', required=True,  type=str, help='deepSNV output folder')
required.add_argument('-depth', dest='DEPTH', nargs='?', required=True,  type=str, help='depth file folder')
required.add_argument('-ref', dest='REF', nargs='?', required=True, type=str, help='reference genome folder')
required.add_argument('-deplen', dest='DEPLEN', nargs='?', required=True, type=int, help='flanking length for mean depth calculation')
required.add_argument('-depdiff', dest='DEPDIFF', nargs='?', required=True, type=int, help='depth difference cutoff')

optional.add_argument('-k', dest='Kmer', nargs='?', required=False, type=int, default=100, help='k-mer mean depth')
optional.add_argument('-o', dest='Out', nargs='?', required=False, type=str, default=None, help='output plot name')
optional.add_argument('-l', dest='Lines', nargs='?', required=False, type=str, default=None, help='output plot name')


args = vars(parser.parse_args())
deepSNV_output_folder = args['SNV']
depth_file_folder = args['DEPTH']
reference_genome_folder = args['REF']
flanking_length_for_mean_depth = args['DEPLEN']
depth_difference_cutoff = args['DEPDIFF']


wd = '/Users/songweizhi/Desktop/666666'
deepSNV_output_folder = '/Users/songweizhi/Desktop/666666/4_1_deepSNV'
reference_genome_folder = '/Users/songweizhi/Desktop/666666/0_References'
depth_file_folder = '/Users/songweizhi/Desktop/666666/depth_files'

min_var_reads_num = 10
min_at_each_direction = 1
strand_bias_cutoff = 20
#freq_difference_cutoff = 10
depth_flanking_length = 50000
depth_kmer = 1000
flanking_length_for_mean_depth = 5000
depth_difference_cutoff = 20

os.chdir(wd)

# example cmd
# cd /Users/songweizhi/Desktop/666666
# python3 ~/PycharmProjects/FlowCellBiofilm/get_SNV_QC.py -snv 4_1_deepSNV -depth depth_files -ref 0_References -deplen 5000 -depdiff 20
# python3 ~/PycharmProjects/FlowCellBiofilm/get_SNV_QC.py -snv -depth -ref -deplen -depdiff


################################################### experiment design ##################################################

matrix_header_210 = ['1D9', '1D18', '1D27', '1D42', '5D9', '5D18', '5D27', '5D42', '9D9', '9D18', '9D27', '9D42', '4D9', '4D18', '4D27', '4D42', '8D9', '8D18', '8D27', '8D42', '12D9', '12D18', '12D27', '12D42']
matrix_header_D2 = ['2D9', '2D18', '2D27', '2D42', '6D9', '6D18', '6D27', '6D42', '10D9', '10D18', '10D27', '10D42', '4D9', '4D18', '4D27', '4D42', '8D9', '8D18', '8D27', '8D42', '12D9', '12D18', '12D27', '12D42']

ref_length_dict = {'2.10_chromosome': 3758219,
                   '2.10_plasmid1': 237747,
                   '2.10_plasmid2': 94490,
                   '2.10_plasmid3': 70353,
                   'D2_c': 4010148,
                   'D2_p': 1010395}


############################################ define the name of output files ###########################################


output_folder = 'output_f%sbp_%smer_md%sbp_c%s' % (depth_flanking_length, depth_kmer, flanking_length_for_mean_depth, depth_difference_cutoff)

deepSNV_output_combined = 'deepSNV_output_combined.txt'
deepSNV_output_combined_QC = 'deepSNV_output_combined_QC.txt'
combined_existence_210 = 'deepSNV_output_summary_210_existence.txt'
combined_existence_D2 =  'deepSNV_output_summary_D2_existence.txt'
combined_frequency_210 = 'deepSNV_output_summary_210_frequency.txt'
combined_frequency_D2 =  'deepSNV_output_summary_D2_frequency.txt'
combined_existence_210_cdc = 'deepSNV_output_summary_210_existence_cdc.txt'
combined_existence_D2_cdc =  'deepSNV_output_summary_D2_existence_cdc.txt'
combined_frequency_210_cdc = 'deepSNV_output_summary_210_frequency_cdc.txt'
combined_frequency_D2_cdc =  'deepSNV_output_summary_D2_frequency_cdc.txt'
SNV_with_even_depth = 'SNV_with_even_depth_f%sbp_%smer_md%sbp' % (depth_flanking_length, depth_kmer, flanking_length_for_mean_depth)
SNV_with_diff_depth = 'SNV_with_diff_depth_f%sbp_%smer_md%sbp' % (depth_flanking_length, depth_kmer, flanking_length_for_mean_depth)

pwd_deepSNV_output_combined = '%s/%s' % (output_folder, deepSNV_output_combined)
pwd_deepSNV_output_combined_QC = '%s/%s' % (output_folder, deepSNV_output_combined_QC)
pwd_combined_existence_210 = '%s/%s' % (output_folder, combined_existence_210)
pwd_combined_existence_D2 =  '%s/%s' % (output_folder, combined_existence_D2)
pwd_combined_frequency_210 = '%s/%s' % (output_folder, combined_frequency_210)
pwd_combined_frequency_D2 =  '%s/%s' % (output_folder, combined_frequency_D2)
pwd_combined_existence_210_cdc = '%s/%s' % (output_folder, combined_existence_210_cdc)
pwd_combined_existence_D2_cdc =  '%s/%s' % (output_folder, combined_existence_D2_cdc)
pwd_combined_frequency_210_cdc = '%s/%s' % (output_folder, combined_frequency_210_cdc)
pwd_combined_frequency_D2_cdc =  '%s/%s' % (output_folder, combined_frequency_D2_cdc)
pwd_SNV_with_even_depth = '%s/%s' % (output_folder, SNV_with_even_depth)
pwd_SNV_with_diff_depth = '%s/%s' % (output_folder, SNV_with_diff_depth)


# create depth_plot_with_coverage_change_folder
if os.path.isdir(output_folder):
    shutil.rmtree(output_folder, ignore_errors=True)
    if os.path.isdir(output_folder):
        shutil.rmtree(output_folder, ignore_errors=True)
    os.makedirs(output_folder)
    os.makedirs('%s/%s' % (output_folder, SNV_with_even_depth))
    os.makedirs('%s/%s' % (output_folder, SNV_with_diff_depth))

else:
    os.makedirs(output_folder)
    os.makedirs('%s/%s' % (output_folder, SNV_with_even_depth))
    os.makedirs('%s/%s' % (output_folder, SNV_with_diff_depth))


################################################ combine deepSNV output ################################################

# get deepSNV output file list
deepSNV_output_file_re = '%s/*.txt' % deepSNV_output_folder
deepSNV_output_file_list = [os.path.basename(file_name) for file_name in glob.glob(deepSNV_output_file_re)]

# combine all deepSNV output in one file
deepSNV_output_combined_handle = open(pwd_deepSNV_output_combined, 'w')
for each_deepSNV_output in deepSNV_output_file_list:
    pwd_each_deepSNV_output = '%s/%s' % (deepSNV_output_folder, each_deepSNV_output)
    treatment_id = each_deepSNV_output.split('_vs_')[0]
    for each_snv in open(pwd_each_deepSNV_output):
        # ignore the first line
        if not each_snv.startswith('chr,pos,ref,var'):
            deepSNV_output_combined_handle.write('%s,%s' % (treatment_id, each_snv))
deepSNV_output_combined_handle.close()


########################################## quality filtering of detected SNVs ##########################################

snv_num_overall = 0
snv_freq_dict = {}
snv_cov_dict = {}
qualified_snv_uniq = []
snv_occurrence_dict = {}
snv_edge_value_dict = {}
snv_edge_category_dict = {}
deepSNV_output_combined_QC_handle = open(pwd_deepSNV_output_combined_QC, 'w')

deepSNV_output_combined_QC_handle.write('sample,chr,pos,ref,var,p_val,freq_var,n_tst_b,n_tst_fw,n_tst_bw,strand_bias,mean_depth_len,mean_depth_left,mean_depth_right,mean_depth_diff,mean_depth_cate\n')
for each_snv in open(pwd_deepSNV_output_combined):
    snv_num_overall += 1
    each_snv_split = each_snv.strip().split(',')
    each_snv_sample = each_snv_split[0]
    each_snv_chr = each_snv_split[1]
    each_snv_pos = int(each_snv_split[2])
    each_snv_ref = each_snv_split[3]
    each_snv_var = each_snv_split[4]
    each_snv_p_val = float("{0:.3f}".format(float(each_snv_split[5])))
    each_snv_freq = float("{0:.3f}".format(float(each_snv_split[6])))
    n_tst_fw = int(each_snv_split[8])
    cov_tst_fw = int(each_snv_split[9])
    n_tst_bw = int(each_snv_split[10])
    cov_tst_bw = int(each_snv_split[11])

    each_snv_key_no_sample_id = '%s|%s|%s|%s' % (each_snv_chr, each_snv_pos, each_snv_ref, each_snv_var)
    each_snv_key_with_sample_id = '%s|%s|%s|%s|%s' % (each_snv_sample, each_snv_chr, each_snv_pos, each_snv_ref, each_snv_var)
    each_snv_cov_both_direction = cov_tst_fw + cov_tst_bw
    n_test_both_direction = n_tst_fw + n_tst_bw
    strand_bias_value = abs((n_tst_fw / n_test_both_direction) - (cov_tst_fw / each_snv_cov_both_direction)) * 100
    strand_bias_value = float("{0:.2f}".format(strand_bias_value))

    if (n_test_both_direction >= min_var_reads_num) and (n_tst_fw >= 1) and (n_tst_bw >= 1) and (strand_bias_value <= strand_bias_cutoff):

        # get snv_freq_dict
        snv_freq_dict[each_snv_key_with_sample_id] = each_snv_freq

        # get snv_cov_dict
        snv_cov_dict[each_snv_key_with_sample_id] = each_snv_cov_both_direction

        # get snv_occurrence_dict
        if each_snv_key_no_sample_id not in snv_occurrence_dict:
            snv_occurrence_dict[each_snv_key_no_sample_id] = [each_snv_sample]
        else:
            snv_occurrence_dict[each_snv_key_no_sample_id].append(each_snv_sample)

        # get qualified_snv_uniq list
        if each_snv_key_no_sample_id not in qualified_snv_uniq:
            qualified_snv_uniq.append(each_snv_key_no_sample_id)


############################################## plot SNV depth and position #############################################

    # plot depth
    pwd_depth_file = '%s/%s.depth' % (depth_file_folder, each_snv_sample)

    # get start position to plot
    plot_start = each_snv_pos - depth_flanking_length
    if plot_start <= 0:
        plot_start = 1

    # get end position to plot
    plot_end = each_snv_pos + depth_flanking_length
    if plot_end >= ref_length_dict[each_snv_chr]:
        plot_end = ref_length_dict[each_snv_chr]

    # get start position to calculate mean depth
    mean_depth_start = each_snv_pos - flanking_length_for_mean_depth
    if mean_depth_start <= 0:
        mean_depth_start = 1

    # get end position to calculate mean depth
    mean_depth_end = each_snv_pos + flanking_length_for_mean_depth
    if mean_depth_end >= ref_length_dict[each_snv_chr]:
        mean_depth_end = ref_length_dict[each_snv_chr]

    # png file name
    each_snv_id_modified = '_'.join(each_snv_key_no_sample_id.split('|'))
    plot_file_name = '%s_%s_f%sbp_%smer' % (each_snv_id_modified, each_snv_sample, depth_flanking_length, depth_kmer)
    pwd_plot_file_name = '%s/%s' % (output_folder, plot_file_name)

    # plot depth
    total_depth_left_side, total_depth_right_side = plot_sam_depth(pwd_depth_file, each_snv_chr, plot_start, plot_end, each_snv_pos, mean_depth_start, mean_depth_end, depth_kmer, str(each_snv_pos), pwd_plot_file_name)

    # mean depth of the two sides
    mean_depth_left = total_depth_left_side/flanking_length_for_mean_depth
    mean_depth_right = total_depth_right_side/flanking_length_for_mean_depth

    # get mean depth difference
    mean_depth_small = sorted([mean_depth_left, mean_depth_right])[0]
    mean_depth_big = sorted([mean_depth_left, mean_depth_right])[1]
    mean_depth_difference = (1 - (mean_depth_small/mean_depth_big))*100
    mean_depth_difference = float("{0:.2f}".format(mean_depth_difference))

    plot_file_name_new = '%s_%s_f%sbp_%smer_md%sbp_%s' % (each_snv_id_modified, each_snv_sample, depth_flanking_length, depth_kmer, flanking_length_for_mean_depth, mean_depth_difference)
    at_edge = 0
    if mean_depth_difference > depth_difference_cutoff:
        at_edge = 1
        os.system('mv %s.png %s/%s.png' % (pwd_plot_file_name, pwd_SNV_with_diff_depth, plot_file_name_new))
    else:
        os.system('mv %s.png %s/%s.png' % (pwd_plot_file_name, pwd_SNV_with_even_depth, plot_file_name_new))

    snv_edge_value_dict[each_snv_key_with_sample_id] = mean_depth_difference
    snv_edge_category_dict[each_snv_key_with_sample_id] = at_edge

    deepSNV_output_combined_QC_handle.write('%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' % (each_snv_sample,
                                                                                                   each_snv_chr,
                                                                                                   each_snv_pos,
                                                                                                   each_snv_ref,
                                                                                                   each_snv_var,
                                                                                                   each_snv_p_val,
                                                                                                   each_snv_freq,
                                                                                                   n_test_both_direction,
                                                                                                   n_tst_fw,
                                                                                                   n_tst_bw,
                                                                                                   strand_bias_value,
                                                                                                   flanking_length_for_mean_depth,
                                                                                                   mean_depth_left,
                                                                                                   mean_depth_right,
                                                                                                   mean_depth_difference,
                                                                                                   at_edge))

deepSNV_output_combined_QC_handle.close()


############################################### get SNV occurrence matrix ##############################################

# sort qualified_snv_uniq list
qualified_snv_uniq_sorted = sorted(qualified_snv_uniq)

# get SNV occurrence matrix for 2.10 and D2
combined_existence_210_handle = open(pwd_combined_existence_210, 'w')
combined_frequency_210_handle = open(pwd_combined_frequency_210, 'w')
combined_existence_D2_handle = open(pwd_combined_existence_D2, 'w')
combined_frequency_D2_handle = open(pwd_combined_frequency_D2, 'w')
combined_existence_210_handle.write('\t%s\n' % '\t'.join(matrix_header_210))
combined_frequency_210_handle.write('\t%s\n' % '\t'.join(matrix_header_210))
combined_existence_D2_handle.write('\t%s\n' % '\t'.join(matrix_header_D2))
combined_frequency_D2_handle.write('\t%s\n' % '\t'.join(matrix_header_D2))

for each_uniq_snv in qualified_snv_uniq_sorted:
    # get matrix for 2.10
    if each_uniq_snv.startswith('2.10'):
        occurrence_list_existence_210 = []
        occurrence_list_frequency_210 = []
        for each_210_sample in matrix_header_210:
            if each_210_sample in snv_occurrence_dict[each_uniq_snv]:
                occurrence_list_existence_210.append('1')
                occurrence_list_frequency_210.append(str(snv_freq_dict['%s|%s' % (each_210_sample, each_uniq_snv)]))
            else:
                occurrence_list_existence_210.append('0')
                occurrence_list_frequency_210.append('0')

        for_write_existence_210 = '%s\t%s\n' % (each_uniq_snv, '\t'.join(occurrence_list_existence_210))
        for_write_frequency_210 = '%s\t%s\n' % (each_uniq_snv, '\t'.join(occurrence_list_frequency_210))
        combined_existence_210_handle.write(for_write_existence_210)
        combined_frequency_210_handle.write(for_write_frequency_210)

    # get matrix for D2
    if each_uniq_snv.startswith('D2'):
        occurrence_list_existence_D2 = []
        occurrence_list_frequency_D2 = []
        for each_D2_sample in matrix_header_D2:
            if each_D2_sample in snv_occurrence_dict[each_uniq_snv]:
                occurrence_list_existence_D2.append('1')
                occurrence_list_frequency_D2.append(str(snv_freq_dict['%s|%s' % (each_D2_sample, each_uniq_snv)]))
            else:
                occurrence_list_existence_D2.append('0')
                occurrence_list_frequency_D2.append('0')

        for_write_existence_D2 = '%s\t%s\n' % (each_uniq_snv, '\t'.join(occurrence_list_existence_D2))
        for_write_frequency_D2 = '%s\t%s\n' % (each_uniq_snv, '\t'.join(occurrence_list_frequency_D2))
        combined_existence_D2_handle.write(for_write_existence_D2)
        combined_frequency_D2_handle.write(for_write_frequency_D2)

combined_existence_210_handle.close()
combined_frequency_210_handle.close()
combined_existence_D2_handle.close()
combined_frequency_D2_handle.close()


############################################### combine continuous deletion ############################################

combine_continuous_deletion(pwd_combined_existence_210, pwd_combined_existence_210_cdc)
combine_continuous_deletion(pwd_combined_frequency_210, pwd_combined_frequency_210_cdc)
combine_continuous_deletion(pwd_combined_existence_D2, pwd_combined_existence_D2_cdc)
combine_continuous_deletion(pwd_combined_frequency_D2, pwd_combined_frequency_D2_cdc)


######################################################## report ########################################################

print('The total number deepSNV detected SNVs: %s' % snv_num_overall)
print('Further validation with the following criteria:')
print('\t1. Significant level set to 0.05')
print('\t2. At least %s reads harbour each SNV' % min_var_reads_num)
print('\t3. At least %s read in each direction' % min_at_each_direction)
print('\t4. Strand bias cut-off set to %s' % strand_bias_cutoff)
print('The number of qualified SNVs after quality filtering: %s' % len(snv_freq_dict))
print('The number of uniq SNVs among all samples: %s' % len(qualified_snv_uniq))
print('The number of uniq SNVs among all samples after combination of continuous deletion: %s' % 'XXX')
