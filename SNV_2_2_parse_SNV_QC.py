import os
import shutil
import argparse
import numpy as np
from scipy import stats
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_ext = os.path.splitext(file_name)

    return file_path, file_basename, file_ext


def get_affect(frequency_list, plot_effect, plot_folder, plot_filename):
    time_point = [9, 18, 27, 42]
    time_point_rescaled = []
    for each_tp in time_point:
        each_tp_rescaled = each_tp / 42
        time_point_rescaled.append(each_tp_rescaled)
    time_point_rescaled_arrary = np.array(time_point_rescaled)

    frequency_list_percentage = []
    for each_p in frequency_list:
        frequency_list_percentage.append(each_p*100)

    frequency_arrary = np.array(frequency_list_percentage)
    min_frequency = frequency_arrary.min()
    max_frequency = frequency_arrary.max()
    slope, intercept, r_value, p_value, std_err = stats.linregress(time_point_rescaled_arrary, frequency_arrary)

    if frequency_arrary[-1] >= 25:
        affect = 'Positive'
    elif (slope > 0) and (max_frequency >= 10):
        affect = 'Positive'
    elif (slope < 0) and (max_frequency <= 10):
        affect = 'Negative'
    else:
        affect = 'Neutral'

    # plot
    if plot_effect == 1:
        plt.plot(time_point_rescaled_arrary, frequency_arrary, 'o')
        plt.plot(time_point_rescaled_arrary, intercept + slope*time_point_rescaled_arrary, 'r')

        # add text
        y_min = plt.ylim()[0]  # get the y-axes minimum value
        y_max = plt.ylim()[1]  # get the y-axes maximum value

        # set text position
        text_x = 0.2
        text_y_slope = y_min + (y_max - y_min) / 5 * 4.4
        text_y_p_value = y_min + (y_max - y_min) / 5 * 4.1
        text_y_affect = y_min + (y_max - y_min) / 5 * 3.8
        plt.text(text_x, text_y_slope, 'Slope: %s' % float("{0:.2f}".format(slope)))
        plt.text(text_x, text_y_p_value, 'P_value: %s' % float("{0:.2f}".format(p_value)))
        plt.text(text_x, text_y_affect, 'Affect: %s' % affect)

        plt.xlabel('Time point')
        plt.ylabel('Frequency (%)')
        plt.savefig('%s/%s.png' % (plot_folder, plot_filename), dpi=300)
        plt.close()

    return affect


def get_mutation_effect(SNV_matrix_cdc, plot_effect):

    # get input file basename
    SNV_matrix_cdc_path, SNV_matrix_cdc_basename, SNV_matrix_cdc_ext = sep_path_basename_ext(SNV_matrix_cdc)
    output_mutation_effect = '%s_mutation_effect.txt' % SNV_matrix_cdc_basename
    plot_folder = '%s_mutation_effect_plot' % SNV_matrix_cdc_basename

    sample_dict = {'1': 'Mono210_A', '5': 'Mono210_B', '9': 'Mono210_C',
                   '2': 'MonoD2_A', '6': 'MonoD2_B', '10': 'MonoD2_C',
                   '4': 'Coculture_A', '8': 'Coculture_B', '12': 'Coculture_C'}

    # create outputs folder
    if plot_effect == 1:
        if os.path.isdir(plot_folder):
            shutil.rmtree(plot_folder, ignore_errors=True)
            if os.path.isdir(plot_folder):
                shutil.rmtree(plot_folder, ignore_errors=True)
            os.makedirs(plot_folder)
        else:
            os.makedirs(plot_folder)

    # get mutation affect
    output_mutation_effect_handle = open(output_mutation_effect, 'w')
    header_no_tp_list = []
    for each_snv in open(SNV_matrix_cdc):
        each_snv_split = each_snv.strip().split('\t')

        # get header list
        if each_snv.startswith('\t'):
            for each in each_snv_split:
                header_no_tp = each.split('D')[0]
                header_no_tp_list.append(header_no_tp)

        if not each_snv.startswith('\t'):

            snv_id = each_snv_split[0]
            occur_profile = []

            if (each_snv_split[1] != '0') or (each_snv_split[2] != '0') or (each_snv_split[3] != '0') or (
                    each_snv_split[4] != '0'):
                frequency_list = [float(each_snv_split[1]), float(each_snv_split[2]), float(each_snv_split[3]),
                                  float(each_snv_split[4])]
                png_filename = '%s_%s' % (snv_id, sample_dict[header_no_tp_list[1 - 1]])
                affect = get_affect(frequency_list, plot_effect, plot_folder, png_filename)
                occur_profile.append('%s_%s' % (header_no_tp_list[1 - 1], affect))

            if (each_snv_split[5] != '0') or (each_snv_split[6] != '0') or (each_snv_split[7] != '0') or (
                    each_snv_split[8] != '0'):
                frequency_list = [float(each_snv_split[5]), float(each_snv_split[6]), float(each_snv_split[7]),
                                  float(each_snv_split[8])]
                png_filename = '%s_%s' % (snv_id, sample_dict[header_no_tp_list[5 - 1]])
                affect = get_affect(frequency_list, plot_effect, plot_folder, png_filename)
                occur_profile.append('%s_%s' % (header_no_tp_list[5 - 1], affect))

            if (each_snv_split[9] != '0') or (each_snv_split[10] != '0') or (each_snv_split[11] != '0') or (
                    each_snv_split[12] != '0'):
                frequency_list = [float(each_snv_split[9]), float(each_snv_split[10]), float(each_snv_split[11]),
                                  float(each_snv_split[12])]
                png_filename = '%s_%s' % (snv_id, sample_dict[header_no_tp_list[9 - 1]])
                affect = get_affect(frequency_list, plot_effect, plot_folder, png_filename)
                occur_profile.append('%s_%s' % (header_no_tp_list[9 - 1], affect))

            if (each_snv_split[13] != '0') or (each_snv_split[14] != '0') or (each_snv_split[15] != '0') or (
                    each_snv_split[16] != '0'):
                frequency_list = [float(each_snv_split[13]), float(each_snv_split[14]), float(each_snv_split[15]),
                                  float(each_snv_split[16])]
                png_filename = '%s_%s' % (snv_id, sample_dict[header_no_tp_list[13 - 1]])
                affect = get_affect(frequency_list, plot_effect, plot_folder, png_filename)
                occur_profile.append('%s_%s' % (header_no_tp_list[13 - 1], affect))

            if (each_snv_split[17] != '0') or (each_snv_split[18] != '0') or (each_snv_split[19] != '0') or (
                    each_snv_split[20] != '0'):
                frequency_list = [float(each_snv_split[17]), float(each_snv_split[18]), float(each_snv_split[19]),
                                  float(each_snv_split[20])]
                png_filename = '%s_%s' % (snv_id, sample_dict[header_no_tp_list[17 - 1]])
                affect = get_affect(frequency_list, plot_effect, plot_folder, png_filename)
                occur_profile.append('%s_%s' % (header_no_tp_list[17 - 1], affect))

            if (each_snv_split[21] != '0') or (each_snv_split[22] != '0') or (each_snv_split[23] != '0') or (
                    each_snv_split[24] != '0'):
                frequency_list = [float(each_snv_split[21]), float(each_snv_split[22]), float(each_snv_split[23]),
                                  float(each_snv_split[24])]
                png_filename = '%s_%s' % (snv_id, sample_dict[header_no_tp_list[21 - 1]])
                affect = get_affect(frequency_list, plot_effect, plot_folder, png_filename)
                occur_profile.append('%s_%s' % (header_no_tp_list[21 - 1], affect))

            output_mutation_effect_handle.write('%s\t%s\n' % (snv_id, '|'.join(occur_profile)))


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


def get_summary_matrix(SNV_file_in, plot_mutation_effect):
    uniq_snv_list = []
    snv_freq_dict = {}
    snv_occurrence_dict = {}
    for snv in open(SNV_file_in):
        if not snv.startswith('sample'):
            snv_split = snv.strip().split(',')
            snv_sample = snv_split[0]
            snv_chr = snv_split[1]
            snv_pos = int(snv_split[2])
            snv_ref = snv_split[3]
            snv_var = snv_split[4]
            snv_freq = float(snv_split[6])

            snv_key_no_sample_id = '%s|%s|%s|%s' % (snv_chr, snv_pos, snv_ref, snv_var)
            snv_key_with_sample_id = '%s|%s|%s|%s|%s' % (snv_sample, snv_chr, snv_pos, snv_ref, snv_var)

            # get snv_freq_dict
            snv_freq_dict[snv_key_with_sample_id] = snv_freq

            # get uniq_snv_list
            if snv_key_no_sample_id not in uniq_snv_list:
                uniq_snv_list.append(snv_key_no_sample_id)

            # get snv_occurrence_dict
            if snv_key_no_sample_id not in snv_occurrence_dict:
                snv_occurrence_dict[snv_key_no_sample_id] = [snv_sample]
            else:
                snv_occurrence_dict[snv_key_no_sample_id].append(snv_sample)

    # sort qualified_snv_uniq list
    uniq_snv_list_sorted = sorted(uniq_snv_list)

    SNV_file_in_basename, SNV_file_in_ext = os.path.splitext(SNV_file_in)

    SNV_file_out_existence_210 = '%s_matrix_210_existence.txt' % SNV_file_in_basename
    SNV_file_out_frequency_210 = '%s_matrix_210_frequency.txt' % SNV_file_in_basename
    SNV_file_out_existence_D2 = '%s_matrix_D2_existence.txt' % SNV_file_in_basename
    SNV_file_out_frequency_D2 = '%s_matrix_D2_frequency.txt' % SNV_file_in_basename

    SNV_file_out_existence_210_cdc = '%s_matrix_210_existence_cdc.txt' % SNV_file_in_basename
    SNV_file_out_frequency_210_cdc = '%s_matrix_210_frequency_cdc.txt' % SNV_file_in_basename
    SNV_file_out_existence_D2_cdc = '%s_matrix_D2_existence_cdc.txt' % SNV_file_in_basename
    SNV_file_out_frequency_D2_cdc = '%s_matrix_D2_frequency_cdc.txt' % SNV_file_in_basename

    # get SNV occurrence matrix for 2.10 and D2
    SNV_file_out_existence_210_handle = open(SNV_file_out_existence_210, 'w')
    SNV_file_out_frequency_210_handle = open(SNV_file_out_frequency_210, 'w')
    SNV_file_out_existence_D2_handle = open(SNV_file_out_existence_D2, 'w')
    SNV_file_out_frequency_D2_handle = open(SNV_file_out_frequency_D2, 'w')

    SNV_file_out_existence_210_handle.write('\t%s\n' % '\t'.join(matrix_header_210))
    SNV_file_out_frequency_210_handle.write('\t%s\n' % '\t'.join(matrix_header_210))
    SNV_file_out_existence_D2_handle.write('\t%s\n' % '\t'.join(matrix_header_D2))
    SNV_file_out_frequency_D2_handle.write('\t%s\n' % '\t'.join(matrix_header_D2))

    # get matrix
    for each_uniq_snv in uniq_snv_list_sorted:
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
            SNV_file_out_existence_210_handle.write(for_write_existence_210)
            SNV_file_out_frequency_210_handle.write(for_write_frequency_210)

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
            SNV_file_out_existence_D2_handle.write(for_write_existence_D2)
            SNV_file_out_frequency_D2_handle.write(for_write_frequency_D2)

    SNV_file_out_existence_210_handle.close()
    SNV_file_out_frequency_210_handle.close()
    SNV_file_out_existence_D2_handle.close()
    SNV_file_out_frequency_D2_handle.close()

    # combine continuous deletion
    combine_continuous_deletion(SNV_file_out_existence_210, SNV_file_out_existence_210_cdc)
    combine_continuous_deletion(SNV_file_out_frequency_210, SNV_file_out_frequency_210_cdc)
    combine_continuous_deletion(SNV_file_out_existence_D2, SNV_file_out_existence_D2_cdc)
    combine_continuous_deletion(SNV_file_out_frequency_D2, SNV_file_out_frequency_D2_cdc)

    # delete tmp files
    os.system('rm %s' % SNV_file_out_existence_210)
    os.system('rm %s' % SNV_file_out_frequency_210)
    os.system('rm %s' % SNV_file_out_existence_D2)
    os.system('rm %s' % SNV_file_out_frequency_D2)

    # get_mutation_effect
    print('Get mutation effect for %s' % SNV_file_out_frequency_210_cdc)
    get_mutation_effect(SNV_file_out_frequency_210_cdc, plot_mutation_effect)

    print('Get mutation effect for %s' % SNV_file_out_frequency_D2_cdc)
    get_mutation_effect(SNV_file_out_frequency_D2_cdc, plot_mutation_effect)


########################################## specify input files and parameters ##########################################

parser = argparse.ArgumentParser(description='', add_help=False)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')

optional.add_argument('-h', action='help', help='Show this help message and exit')
required.add_argument('-snv_qc', dest='SNV_QC', nargs='?', required=True,  type=str, help='deepSNV QC file')
required.add_argument('-min_both', dest='MIN_BOTH', nargs='?', required=True, type=int, help='The minimum number of reads harboring SNV')
required.add_argument('-min_each', dest='MIN_EACH', nargs='?', required=True, type=int, help='The minimum number of reads harboring SNV at each direction')
required.add_argument('-strand_bias', dest='STRAND_BIAS', nargs='?', required=True, type=int, help='strand_bias cutoff')
required.add_argument('-depth_diff', dest='DEPTH_DIFF', nargs='?', required=True,  type=int, help='depth difference cutoff')
required.add_argument('-deplen', dest='DEPLEN', nargs='?', required=True, type=int, help='flanking length for mean depth calculation')
required.add_argument('-plot', dest='PLOT', nargs='?', required=True, type=int, help='plot mutation effect')

args = vars(parser.parse_args())
SNV_quality_file = args['SNV_QC']
depth_difference_cutoff = args['DEPTH_DIFF']
strand_bias_cutoff = args['STRAND_BIAS']
mean_depth_len = args['DEPLEN']
min_var_reads_num = args['MIN_BOTH']
min_at_each_direction = args['MIN_EACH']
plot_mutation_effect = args['PLOT']


############################################ define the name of output files ###########################################

matrix_header_210 = ['1D9', '1D18', '1D27', '1D42', '5D9', '5D18', '5D27', '5D42', '9D9', '9D18', '9D27', '9D42', '4D9', '4D18', '4D27', '4D42', '8D9', '8D18', '8D27', '8D42', '12D9', '12D18', '12D27', '12D42']
matrix_header_D2 = ['2D9', '2D18', '2D27', '2D42', '6D9', '6D18', '6D27', '6D42', '10D9', '10D18', '10D27', '10D42', '4D9', '4D18', '4D27', '4D42', '8D9', '8D18', '8D27', '8D42', '12D9', '12D18', '12D27', '12D42']

SNV_quality_file_path, SNV_quality_file_name = os.path.split(SNV_quality_file)
SNV_quality_file_basename, SNV_quality_file_ext = os.path.splitext(SNV_quality_file_name)

qualified_SNVs_even_flanking_depth_file = '%s_even_depth.txt' % SNV_quality_file_basename
qualified_SNVs_diff_flanking_depth_file = '%s_diff_depth.txt' % SNV_quality_file_basename


############################################ parse quality of detected SNVs ############################################

n_tst_total_unqualified = []
n_tst_each_unqualified = []
strand_bias_unqualified = []
depth_diff_unqualified = []
qualified_SNVs_even_flanking_depth = []
qualified_SNVs_diff_flanking_depth = []
total_num = 0

qualified_SNVs_even_flanking_depth_file_handle = open(qualified_SNVs_even_flanking_depth_file, 'w')
qualified_SNVs_diff_flanking_depth_file_handle = open(qualified_SNVs_diff_flanking_depth_file, 'w')

for each_snv in open(SNV_quality_file):

    if each_snv.startswith('sample'):
        qualified_SNVs_even_flanking_depth_file_handle.write(each_snv)
        qualified_SNVs_diff_flanking_depth_file_handle.write(each_snv)

    if not each_snv.startswith('sample'):
        each_snv_split = each_snv.strip().split(',')
        each_snv_sample = each_snv_split[0]
        each_snv_chr = each_snv_split[1]
        each_snv_pos = int(each_snv_split[2])
        each_snv_ref = each_snv_split[3]
        each_snv_var = each_snv_split[4]
        each_snv_p_val = float(each_snv_split[5])
        each_snv_freq = float(each_snv_split[6])
        each_snv_n_tst_b = int(each_snv_split[7])
        each_snv_n_tst_fw = int(each_snv_split[8])
        each_snv_n_tst_bw = int(each_snv_split[9])
        each_snv_strand_bias = float(each_snv_split[10])
        each_snv_mean_depth_diff = float(each_snv_split[14])

        each_snv_key_no_sample_id = '%s|%s|%s|%s' % (each_snv_chr, each_snv_pos, each_snv_ref, each_snv_var)
        each_snv_key_with_sample_id = '%s|%s|%s|%s|%s' % (each_snv_sample, each_snv_chr, each_snv_pos, each_snv_ref, each_snv_var)

        # get SNVs with unqualified number of reads harbouring SNV
        if each_snv_n_tst_b < min_var_reads_num:
            n_tst_total_unqualified.append(each_snv_key_with_sample_id)

        # get SNVs with unqualified number of reads in each direction
        if (each_snv_n_tst_fw < min_at_each_direction) or (each_snv_n_tst_bw < min_at_each_direction):
            n_tst_each_unqualified.append(each_snv_key_with_sample_id)

        # get SNVs with unqualified strand bias
        if each_snv_strand_bias > strand_bias_cutoff:
            strand_bias_unqualified.append(each_snv_key_with_sample_id)

        # get SNVs with unqualified flanking depth difference
        if each_snv_mean_depth_diff > depth_difference_cutoff:
            depth_diff_unqualified.append(each_snv_key_with_sample_id)

        # get qualified SNV with similar flanking depth
        if (each_snv_n_tst_b >= min_var_reads_num) and (each_snv_n_tst_fw >= min_at_each_direction) and (each_snv_n_tst_bw >= min_at_each_direction) and (each_snv_strand_bias <= strand_bias_cutoff) and (each_snv_mean_depth_diff <= depth_difference_cutoff):
            qualified_SNVs_even_flanking_depth.append(each_snv_key_with_sample_id)
            qualified_SNVs_even_flanking_depth_file_handle.write(each_snv)

        # get qualified SNV with different flanking depth
        if (each_snv_n_tst_b >= min_var_reads_num) and (each_snv_n_tst_fw >= min_at_each_direction) and (each_snv_n_tst_bw >= min_at_each_direction) and (each_snv_strand_bias <= strand_bias_cutoff) and (each_snv_mean_depth_diff > depth_difference_cutoff):
            qualified_SNVs_diff_flanking_depth.append(each_snv_key_with_sample_id)
            qualified_SNVs_diff_flanking_depth_file_handle.write(each_snv)

        total_num += 1

qualified_SNVs_even_flanking_depth_file_handle.close()
qualified_SNVs_diff_flanking_depth_file_handle.close()


# report
print('The total number of detected SNVs: %s' % total_num)
print('The number of SNVs with less than %s reads harboring it: %s' % (min_var_reads_num, len(n_tst_total_unqualified)))
print('The number of SNVs with reads only from one direction: %s' % len(n_tst_each_unqualified))
print('The number of SNVs with strand bias higher than %s%s: %s' % (strand_bias_cutoff, '%',len(strand_bias_unqualified)))
print('The number of SNVs with flanking depth (%sbp) difference higher than %s%s: %s' % (mean_depth_len, depth_difference_cutoff, '%', len(depth_diff_unqualified)))
print('The number of qualified SNVs with similar flanking depth: %s' % len(qualified_SNVs_even_flanking_depth))
print('The number of qualified SNVs with different flanking depth: %s' % len(qualified_SNVs_diff_flanking_depth))
print()

############################## get matrix for SNVs with similar/different flanking depth ###############################

get_summary_matrix(qualified_SNVs_even_flanking_depth_file, plot_mutation_effect)
get_summary_matrix(qualified_SNVs_diff_flanking_depth_file, plot_mutation_effect)

# os.system('rm *_existence_cdc.txt')
# os.system('rm %s' % qualified_SNVs_even_flanking_depth_file)
# os.system('rm %s' % qualified_SNVs_diff_flanking_depth_file)
