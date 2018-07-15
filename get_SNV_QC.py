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


########################################## specify input files and parameters ##########################################

parser = argparse.ArgumentParser(description='', add_help=False)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')

optional.add_argument('-h', action='help', help='Show this help message and exit')
required.add_argument('-snv', dest='SNV', nargs='?', required=True,  type=str, help='deepSNV output folder')
required.add_argument('-depth', dest='DEPTH', nargs='?', required=True,  type=str, help='depth file folder')
required.add_argument('-ref', dest='REF', nargs='?', required=True, type=str, help='reference genome folder')
required.add_argument('-deplen', dest='DEPLEN', nargs='?', required=True, type=int, help='flanking length for mean depth calculation')

args = vars(parser.parse_args())
deepSNV_output_folder = args['SNV']
depth_file_folder = args['DEPTH']
reference_genome_folder = args['REF']
flanking_length_for_mean_depth = args['DEPLEN']

depth_flanking_length = 50000
depth_kmer = 1000


# example cmd
# cd /Users/songweizhi/Desktop/666666
# python3 ~/PycharmProjects/FlowCellBiofilm/get_SNV_QC.py -snv 4_1_deepSNV -depth depth_files -ref 0_References -deplen 5000
# python3 ~/PycharmProjects/FlowCellBiofilm/get_SNV_QC.py -snv -depth -ref -deplen -depdiff


################################################### experiment design ##################################################

ref_length_dict = {'2.10_chromosome': 3758219,
                   '2.10_plasmid1': 237747,
                   '2.10_plasmid2': 94490,
                   '2.10_plasmid3': 70353,
                   'D2_c': 4010148,
                   'D2_p': 1010395}


############################################ define the name of output files ###########################################

output_folder = 'output_f%sbp_%smer_dl%sbp' % (depth_flanking_length, depth_kmer, flanking_length_for_mean_depth)

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
SNV_depth_plot = 'SNV_depth_plot'

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
pwd_SNV_depth_plot = '%s/%s' % (output_folder, SNV_depth_plot)


# create depth_plot_with_coverage_change_folder
if os.path.isdir(output_folder):
    shutil.rmtree(output_folder, ignore_errors=True)
    if os.path.isdir(output_folder):
        shutil.rmtree(output_folder, ignore_errors=True)
    os.makedirs(output_folder)
    os.makedirs(pwd_SNV_depth_plot)

else:
    os.makedirs(output_folder)
    os.makedirs(pwd_SNV_depth_plot)


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
    pwd_plot_file_name = '%s/%s' % (pwd_SNV_depth_plot, plot_file_name)

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
    pwd_plot_file_name_new = '%s/%s' % (pwd_SNV_depth_plot, plot_file_name_new)

    os.system('mv %s.png %s.png' % (pwd_plot_file_name, pwd_plot_file_name_new))

    deepSNV_output_combined_QC_handle.write('%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' % (each_snv_sample,
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
                                                                                                   mean_depth_difference))
deepSNV_output_combined_QC_handle.close()
