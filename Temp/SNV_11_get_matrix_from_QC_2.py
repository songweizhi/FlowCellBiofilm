import os
import shutil
import argparse
import numpy as np
from scipy import stats
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna


def merge_two_dict(dict_a, dict_b):
    dict_c = dict_a.copy()
    dict_c.update(dict_b)

    return dict_c


def unique_list_elements(list_input):

    list_output = []
    for each_element in list_input:
        if each_element not in list_output:
            list_output.append(each_element)

    return list_output


def remove_l2_elements_from_l1(l1, l2):

    l1_new = []
    for each in l1:
        if each not in l2:
            l1_new.append(each)
    return l1_new


def write_out(op_file, list):

    op_file_handle = open(op_file, 'w')

    for each in list:
        op_file_handle.write('%s\n' % each)

    op_file_handle.close()


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_ext = os.path.splitext(file_name)

    return file_path, file_basename, file_ext


def get_rf_pos_list(gene_len):

    rf_pos_list = []
    n = 0
    while n < (gene_len / 3):
        start_pos = 1 + 3 * n
        rf_pos_list.append([start_pos, start_pos + 1, start_pos + 2])
        n += 1

    return rf_pos_list


def get_mutation_cate_summary(file_in, file_out):

    mutation_cate_dict = {}
    for each_snv in open(file_in):
        each_snv_split = each_snv.strip().split('\t')
        mutated_gene = each_snv_split[3]
        mutation_cate = each_snv_split[5]

        if mutated_gene != 'NA':
            if mutated_gene not in mutation_cate_dict:
                mutation_cate_dict[mutated_gene] = [mutation_cate]
            else:
                mutation_cate_dict[mutated_gene].append(mutation_cate)

    mutation_cate_list = ['Missense', 'Nonsense', 'Silent', 'Fragment_deletion', 'Frameshift']

    file_out_handle = open(file_out, 'w')
    file_out_handle.write('Gene\tMis\tNon\tSilen\tFD\tFS\n')
    for each_gene in mutation_cate_dict:

        mutation_cate_occurence = []
        for each_cate in mutation_cate_list:
            mutation_cate_occurence.append(str(mutation_cate_dict[each_gene].count(each_cate)))
        file_out_handle.write('%s\t%s\n' % (each_gene, '\t'.join(mutation_cate_occurence)))

    file_out_handle.close()

    return mutation_cate_dict


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


def get_matrix(SNV_file_in, CCD):
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

    # SNV_file_out_existence_210 = '%s_matrix_210_existence.txt' % SNV_file_in_basename
    SNV_file_out_frequency_210 = '%s_matrix_210_frequency.txt' % SNV_file_in_basename
    # SNV_file_out_existence_D2 = '%s_matrix_D2_existence.txt' % SNV_file_in_basename
    SNV_file_out_frequency_D2 = '%s_matrix_D2_frequency.txt' % SNV_file_in_basename

    #SNV_file_out_existence_210_cdc = '%s_matrix_210_existence_cdc.txt' % SNV_file_in_basename
    SNV_file_out_frequency_210_cdc = '%s_matrix_210_frequency_cdc.txt' % SNV_file_in_basename
    #SNV_file_out_existence_D2_cdc = '%s_matrix_D2_existence_cdc.txt' % SNV_file_in_basename
    SNV_file_out_frequency_D2_cdc = '%s_matrix_D2_frequency_cdc.txt' % SNV_file_in_basename

    # get SNV occurrence matrix for 2.10 and D2
    # SNV_file_out_existence_210_handle = open(SNV_file_out_existence_210, 'w')
    SNV_file_out_frequency_210_handle = open(SNV_file_out_frequency_210, 'w')
    # SNV_file_out_existence_D2_handle = open(SNV_file_out_existence_D2, 'w')
    SNV_file_out_frequency_D2_handle = open(SNV_file_out_frequency_D2, 'w')

    # SNV_file_out_existence_210_handle.write('\t%s\n' % '\t'.join(matrix_header_210))
    SNV_file_out_frequency_210_handle.write('\t%s\n' % '\t'.join(matrix_header_210))
    # SNV_file_out_existence_D2_handle.write('\t%s\n' % '\t'.join(matrix_header_D2))
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
            # SNV_file_out_existence_210_handle.write(for_write_existence_210)
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
            # SNV_file_out_existence_D2_handle.write(for_write_existence_D2)
            SNV_file_out_frequency_D2_handle.write(for_write_frequency_D2)

    # SNV_file_out_existence_210_handle.close()
    SNV_file_out_frequency_210_handle.close()
    # SNV_file_out_existence_D2_handle.close()
    SNV_file_out_frequency_D2_handle.close()

    # combine continuous deletion
    if CCD == 1:
        # combine_continuous_deletion(SNV_file_out_existence_210, SNV_file_out_existence_210_cdc)
        combine_continuous_deletion(SNV_file_out_frequency_210, SNV_file_out_frequency_210_cdc)
        # combine_continuous_deletion(SNV_file_out_existence_D2, SNV_file_out_existence_D2_cdc)
        combine_continuous_deletion(SNV_file_out_frequency_D2, SNV_file_out_frequency_D2_cdc)

        # delete tmp files
        # os.system('rm %s' % SNV_file_out_existence_210)
        os.system('rm %s' % SNV_file_out_frequency_210)
        # os.system('rm %s' % SNV_file_out_existence_D2)
        os.system('rm %s' % SNV_file_out_frequency_D2)


def get_summary_matrix(SNV_file_in, CCD, plot_mutation_effect):

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

    #SNV_file_out_existence_210 = '%s_matrix_210_existence.txt' % SNV_file_in_basename
    SNV_file_out_frequency_210 = '%s_matrix_210_frequency.txt' % SNV_file_in_basename
    #SNV_file_out_existence_D2 = '%s_matrix_D2_existence.txt' % SNV_file_in_basename
    SNV_file_out_frequency_D2 = '%s_matrix_D2_frequency.txt' % SNV_file_in_basename

    SNV_file_out_existence_210_cdc = '%s_matrix_210_existence_cdc.txt' % SNV_file_in_basename
    SNV_file_out_frequency_210_cdc = '%s_matrix_210_frequency_cdc.txt' % SNV_file_in_basename
    SNV_file_out_existence_D2_cdc = '%s_matrix_D2_existence_cdc.txt' % SNV_file_in_basename
    SNV_file_out_frequency_D2_cdc = '%s_matrix_D2_frequency_cdc.txt' % SNV_file_in_basename

    # get SNV occurrence matrix for 2.10 and D2
    #SNV_file_out_existence_210_handle = open(SNV_file_out_existence_210, 'w')
    SNV_file_out_frequency_210_handle = open(SNV_file_out_frequency_210, 'w')
    #SNV_file_out_existence_D2_handle = open(SNV_file_out_existence_D2, 'w')
    SNV_file_out_frequency_D2_handle = open(SNV_file_out_frequency_D2, 'w')

    #SNV_file_out_existence_210_handle.write('\t%s\n' % '\t'.join(matrix_header_210))
    SNV_file_out_frequency_210_handle.write('\t%s\n' % '\t'.join(matrix_header_210))
    #SNV_file_out_existence_D2_handle.write('\t%s\n' % '\t'.join(matrix_header_D2))
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
            #SNV_file_out_existence_210_handle.write(for_write_existence_210)
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
            #SNV_file_out_existence_D2_handle.write(for_write_existence_D2)
            SNV_file_out_frequency_D2_handle.write(for_write_frequency_D2)

    #SNV_file_out_existence_210_handle.close()
    SNV_file_out_frequency_210_handle.close()
    #SNV_file_out_existence_D2_handle.close()
    SNV_file_out_frequency_D2_handle.close()

    # combine continuous deletion
    if CCD == 1:
        #combine_continuous_deletion(SNV_file_out_existence_210, SNV_file_out_existence_210_cdc)
        combine_continuous_deletion(SNV_file_out_frequency_210, SNV_file_out_frequency_210_cdc)
        #combine_continuous_deletion(SNV_file_out_existence_D2, SNV_file_out_existence_D2_cdc)
        combine_continuous_deletion(SNV_file_out_frequency_D2, SNV_file_out_frequency_D2_cdc)

        # delete tmp files
        #os.system('rm %s' % SNV_file_out_existence_210)
        os.system('rm %s' % SNV_file_out_frequency_210)
        #os.system('rm %s' % SNV_file_out_existence_D2)
        os.system('rm %s' % SNV_file_out_frequency_D2)

    # get_mutation_effect
    if CCD == 0:
        # print('Get mutation effect for %s' % SNV_file_out_frequency_210_cdc)
        get_mutation_effect(SNV_file_out_frequency_210, plot_mutation_effect)

        # print('Get mutation effect for %s' % SNV_file_out_frequency_D2_cdc)
        get_mutation_effect(SNV_file_out_frequency_D2, plot_mutation_effect)

    if CCD == 1:
        # print('Get mutation effect for %s' % SNV_file_out_frequency_210_cdc)
        get_mutation_effect(SNV_file_out_frequency_210_cdc, plot_mutation_effect)

        # print('Get mutation effect for %s' % SNV_file_out_frequency_D2_cdc)
        get_mutation_effect(SNV_file_out_frequency_D2_cdc, plot_mutation_effect)


def check_existence(num_list_1, num_list_2):

    existence_list = []
    for each_num in num_list_1:
        if each_num in num_list_2:
            existence_list.append(1)
        else:
            existence_list.append(0)

    if 1 in existence_list:
        return True
    else:
        return False


def check_occurrence (list_in):

    current_occurrence = 0
    for each in list_in:
        if each != '0':
            current_occurrence = 1

    return current_occurrence


def remove_0_from_list(list_in):

    no_0_out_list = []
    for each in list_in:
        if each != 0:
            no_0_out_list.append(each)

    return no_0_out_list


def check_parallel(matrix_file):
    parallel_dict = {}
    for snv in open(matrix_file):
        if not snv.startswith('\t'):
            snv_split = snv.strip().split('\t')
            snv_id = snv_split[0]
            mono_A_occurrence = check_occurrence([snv_split[1], snv_split[2], snv_split[3], snv_split[4]])
            mono_B_occurrence = check_occurrence([snv_split[5], snv_split[6], snv_split[7], snv_split[8]])
            mono_C_occurrence = check_occurrence([snv_split[9], snv_split[10], snv_split[11], snv_split[12]])
            co_A_occurrence = check_occurrence([snv_split[13], snv_split[14], snv_split[15], snv_split[16]])
            co_B_occurrence = check_occurrence([snv_split[17], snv_split[18], snv_split[19], snv_split[20]])
            co_C_occurrence = check_occurrence([snv_split[21], snv_split[22], snv_split[23], snv_split[24]])

            mono_occurrence = [mono_A_occurrence, mono_B_occurrence, mono_C_occurrence]
            co_occurrence = [co_A_occurrence, co_B_occurrence, co_C_occurrence]

            # get parallel profile
            parallel_profile = 'NA'

            if (mono_occurrence.count(1) >= 1) and (co_occurrence.count(1) >= 1):
                parallel_profile = 'PMC'

            elif (mono_occurrence.count(1) > 1) and (co_occurrence.count(1) == 0):
                parallel_profile = 'PM'

            elif (mono_occurrence.count(1) == 0) and (co_occurrence.count(1) > 1):
                parallel_profile = 'PC'

            elif (mono_occurrence.count(1) == 1) and (co_occurrence.count(1) == 0):
                parallel_profile = 'NPM'

            elif (mono_occurrence.count(1) == 0) and (co_occurrence.count(1) == 1):
                parallel_profile = 'NPC'

            parallel_dict[snv_id] = parallel_profile


            # print('%s\t%s\t\t%s|%s|%s|%s\t\t%s|%s|%s|%s\t\t%s|%s|%s|%s\t\t%s|%s|%s|%s\t\t%s|%s|%s|%s\t\t%s|%s|%s|%s' % (snv_id, parallel_profile,
            #                                                                                                 snv_split[1], snv_split[2], snv_split[3], snv_split[4],
            #                                                                                                 snv_split[5], snv_split[6], snv_split[7], snv_split[8],
            #                                                                                                 snv_split[9], snv_split[10], snv_split[11], snv_split[12],
            #                                                                                                 snv_split[13], snv_split[14], snv_split[15], snv_split[16],
            #                                                                                                 snv_split[17], snv_split[18], snv_split[19], snv_split[20],
            #                                                                                                 snv_split[21], snv_split[22], snv_split[23], snv_split[24]))

    return parallel_dict


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
required.add_argument('-ccd', dest='CCD', nargs='?', required=True, type=int, help='combine continuous deletion')


args = vars(parser.parse_args())
SNV_quality_file = args['SNV_QC']
depth_difference_cutoff = args['DEPTH_DIFF']
strand_bias_cutoff = args['STRAND_BIAS']
mean_depth_len = args['DEPLEN']
min_var_reads_num = args['MIN_BOTH']
min_at_each_direction = args['MIN_EACH']
plot_mutation_effect = args['PLOT']
CCD = args['CCD']


############################################ define the name of output files ###########################################

matrix_header_210 = ['1D9', '1D18', '1D27', '1D42', '5D9', '5D18', '5D27', '5D42', '9D9', '9D18', '9D27', '9D42', '4D9', '4D18', '4D27', '4D42', '8D9', '8D18', '8D27', '8D42', '12D9', '12D18', '12D27', '12D42']
matrix_header_D2 = ['2D9', '2D18', '2D27', '2D42', '6D9', '6D18', '6D27', '6D42', '10D9', '10D18', '10D27', '10D42', '4D9', '4D18', '4D27', '4D42', '8D9', '8D18', '8D27', '8D42', '12D9', '12D18', '12D27', '12D42']

combined_ref_fasta = '/Users/songweizhi/Dropbox/Research/Flow_cell_datasets/reference_files/combined_ref.fasta'
combined_ref_gff = '/Users/songweizhi/Dropbox/Research/Flow_cell_datasets/reference_files/combined_ref.gff'
combined_ref_ffn = '/Users/songweizhi/Dropbox/Research/Flow_cell_datasets/reference_files/combined_ref.ffn'
combined_ref_faa = '/Users/songweizhi/Dropbox/Research/Flow_cell_datasets/reference_files/combined_ref.faa'
annotation_file = '/Users/songweizhi/Dropbox/Research/Flow_cell_datasets/reference_files/combined_ref_aa.faa.COG.arCOG.kegg'
transl_table = 11
pwd_QC_txt = 'SNV_QC.txt'

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


###################################################### get matrix ######################################################

print('Run get matrix')
print(qualified_SNVs_even_flanking_depth_file)
print(CCD)
get_matrix(qualified_SNVs_even_flanking_depth_file, CCD)


#################################################### check parallel ####################################################

# input: matrix file

SNV_file_in_basename, SNV_file_in_ext = os.path.splitext(qualified_SNVs_even_flanking_depth_file)

SNV_parallel_dict_210 = {}
SNV_parallel_dict_D2 = {}

if CCD == 0:
    SNV_file_out_frequency_210 = '%s_matrix_210_frequency.txt' % SNV_file_in_basename
    SNV_file_out_frequency_D2 = '%s_matrix_D2_frequency.txt' % SNV_file_in_basename
    SNV_parallel_dict_210 = check_parallel(SNV_file_out_frequency_210)
    SNV_parallel_dict_D2 = check_parallel(SNV_file_out_frequency_D2)

if CCD == 1:
    SNV_file_out_frequency_210_cdc = '%s_matrix_210_frequency_cdc.txt' % SNV_file_in_basename
    SNV_file_out_frequency_D2_cdc = '%s_matrix_D2_frequency_cdc.txt' % SNV_file_in_basename
    SNV_parallel_dict_210 = check_parallel(SNV_file_out_frequency_210_cdc)
    SNV_parallel_dict_D2 = check_parallel(SNV_file_out_frequency_D2_cdc)

# merge two dict
SNV_parallel_dict_combined = merge_two_dict(SNV_parallel_dict_210, SNV_parallel_dict_D2)


################################################# functional annotation ################################################

# store annotation results in dicts
gene_COG_cat_dict = {}
gene_COG_id_dict = {}
gene_COG_function_dict = {}
for each_snv3 in open(annotation_file):
    each_snv3_split = each_snv3.strip().split('\t')
    gene_id = each_snv3_split[0]
    COG_cat = each_snv3_split[10]
    COG_id = each_snv3_split[8]
    COG_function = each_snv3_split[9]
    gene_COG_cat_dict[gene_id] = COG_cat
    gene_COG_id_dict[gene_id] = COG_id
    gene_COG_function_dict[gene_id] = COG_function

SNV_matrix_cdc_list = []

if CCD == 0:
    SNV_file_out_frequency_210 = '%s_matrix_210_frequency.txt' % SNV_file_in_basename
    SNV_file_out_frequency_D2 = '%s_matrix_D2_frequency.txt' % SNV_file_in_basename
    SNV_matrix_cdc_list = [SNV_file_out_frequency_210, SNV_file_out_frequency_D2]

if CCD == 1:
    SNV_file_out_frequency_210_cdc = '%s_matrix_210_frequency_cdc.txt' % SNV_file_in_basename
    SNV_file_out_frequency_D2_cdc = '%s_matrix_D2_frequency_cdc.txt' % SNV_file_in_basename
    SNV_matrix_cdc_list = [SNV_file_out_frequency_210_cdc, SNV_file_out_frequency_D2_cdc]


for SNV_matrix_cdc in SNV_matrix_cdc_list:

    # output files
    SNV_matrix_cdc_path, SNV_matrix_cdc_basename, SNV_matrix_cdc_ext = sep_path_basename_ext(SNV_matrix_cdc)
    output_mutated_genes = '%s_summary.txt' % SNV_matrix_cdc_basename
    #output_mutated_genes_cate = '%s_mutated_genes_category.txt' % SNV_matrix_cdc_basename
    #output_summary = '%s_summary.txt' % SNV_matrix_cdc_basename
    # output_seq_nc = '%s_mutated_genes_nc.fasta' % SNV_matrix_cdc_basename
    # output_seq_aa = '%s_mutated_genes_aa.fasta' % SNV_matrix_cdc_basename

    # store sequences in dict
    ref_seq_dict = {}
    for each_seq in SeqIO.parse(combined_ref_fasta, 'fasta'):
        ref_seq_dict[each_seq.id] = str(each_seq.seq)

    # get seq id list
    seq_id_list = []
    for cds in open(combined_ref_gff):
        if (cds.startswith('2.10')) or (cds.startswith('D2')):
            seq_ID = cds.strip().split('\t')[0]
            if seq_ID not in seq_id_list:
                seq_id_list.append(seq_ID)

    # get ending positions for all genes
    ORF_ending_pos_dict = {}
    ORF_seq_id_dict = {}
    ORF_strand_dict = {}
    for cds2 in open(combined_ref_gff):
        if (cds2.startswith('2.10')) or (cds2.startswith('D2')):
            cds2_split = cds2.strip().split('\t')
            seq_id2 = cds2_split[0]
            gene_id2 = cds2_split[8].split(';')[0].split('=')[1]
            start_pos2 = int(cds2_split[3])
            end_pos2 = int(cds2_split[4])
            strand = cds2_split[6]
            ORF_ending_pos_dict[gene_id2] = [start_pos2, end_pos2]
            ORF_seq_id_dict[gene_id2] = seq_id2
            ORF_strand_dict[gene_id2] = strand

    # initialize coding region dict
    coding_region_dict = {}
    for each_id in seq_id_list:
        coding_region_dict[each_id] = []

    # add coding pos to coding_region_dict
    for each_cds in open(combined_ref_gff):
        if (each_cds.startswith('2.10')) or (each_cds.startswith('D2')):
            each_cds_split = each_cds.strip().split('\t')
            seq_id = each_cds_split[0]
            start_pos = int(each_cds_split[3])
            end_pos = int(each_cds_split[4])
            pos_list = list(range(start_pos, end_pos + 1))
            for each_pos in pos_list:
                coding_region_dict[seq_id].append(each_pos)

    output_handle = open(output_mutated_genes, 'w')
    for each_snv in open(SNV_matrix_cdc):

        if not each_snv.startswith('\t'):
            each_snv_id = each_snv.strip().split('\t')[0]
            each_snv_seq = each_snv.strip().split('\t')[0].split('|')[0]
            each_snv_pos = each_snv.strip().split('\t')[0].split('|')[1]
            each_snv_pos_wt = each_snv.strip().split('\t')[0].split('|')[2]
            each_snv_pos_v = each_snv.strip().split('\t')[0].split('|')[3]

            # check parallel
            each_snv_parallel = SNV_parallel_dict_combined[each_snv_id]


            # get all affected genes for
            each_snv_pos = int(each_snv_pos)

            # get SNV location
            location = ''
            if each_snv_pos in coding_region_dict[each_snv_seq]:
                location = 'Coding'
            else:
                location = 'Intergenic'
                output_handle.write('%s\t%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\n' % (each_snv.strip().split('\t')[0], each_snv_parallel, location))
                #output_handle.write('%s\t%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\n' % (each_snv.strip(), each_snv_parallel, location))

            # get mutation type
            if location == 'Coding':
                for each_gene in ORF_ending_pos_dict:
                    if (each_snv_seq == ORF_seq_id_dict[each_gene]) and (ORF_ending_pos_dict[each_gene][0] <= each_snv_pos <= ORF_ending_pos_dict[each_gene][1]):

                        start_pos_raw = ORF_ending_pos_dict[each_gene][0]
                        end_pos_raw = ORF_ending_pos_dict[each_gene][1]
                        snv_pos_raw = each_snv_pos
                        start_pos_rescaled = ORF_ending_pos_dict[each_gene][0] - (ORF_ending_pos_dict[each_gene][0] - 1)
                        end_pos_rescaled = ORF_ending_pos_dict[each_gene][1] - (ORF_ending_pos_dict[each_gene][0] - 1)
                        snv_pos_rescaled = each_snv_pos - (ORF_ending_pos_dict[each_gene][0] - 1)

                        # functional annotation
                        current_COG_id2 = ''
                        current_COG_function2 = ''
                        if each_gene in gene_COG_id_dict:
                            current_COG_id2 = gene_COG_id_dict[each_gene]
                            current_COG_function2 = gene_COG_function_dict[each_gene]
                        else:
                            current_COG_id2 = 'NA'
                            current_COG_function2 = 'NA'

                        mutation_type = ''

                        if each_snv_pos_v == '-':
                            mutation_type = 'Frameshift'

                            output_handle.write('%s\t%s\t%s\t%s\t%s\tNA\t%s\t%s\t%s\n' % (each_snv.strip().split('\t')[0], each_snv_parallel, location, ORF_strand_dict[each_gene], each_gene, mutation_type, current_COG_id2, current_COG_function2))
                            #output_handle.write('%s\t%s\t%s\t%s\t%s\tNA\t%s\t%s\t%s\n' % (each_snv.strip(), each_snv_parallel, location, ORF_strand_dict[each_gene], each_gene, mutation_type, current_COG_id2, current_COG_function2))

                        elif each_snv_pos_v in ['A', 'T', 'C', 'G']:

                            # get all reading frame
                            rf_pos_list = get_rf_pos_list(snv_pos_rescaled)

                            # get the mutated reading frame with rescaled position
                            mutated_rf_rescaled = []
                            for each_rf in rf_pos_list:
                                if snv_pos_rescaled in each_rf:
                                    mutated_rf_rescaled = each_rf

                            # get the mutated reading frame with raw position
                            mutated_rf_raw = []
                            for each_rescaled_pos in mutated_rf_rescaled:
                                raw_pos = each_rescaled_pos + ORF_ending_pos_dict[each_gene][0] - 1
                                mutated_rf_raw.append(raw_pos)

                            # get sequence of raw and mutated reading frame
                            rf_seq_raw = ref_seq_dict[each_snv_seq][(mutated_rf_raw[0] - 1):(mutated_rf_raw[0] + 2)]
                            rf_seq_mutated = ''
                            for each_bp in mutated_rf_raw:
                                current_bp = ''
                                if each_bp == each_snv_pos:
                                    current_bp = each_snv_pos_v
                                else:
                                    current_bp = ref_seq_dict[each_snv_seq][each_bp - 1]
                                rf_seq_mutated += current_bp

                            if ORF_strand_dict[each_gene] == '-':
                                rf_seq_raw = str(Seq(rf_seq_raw, generic_dna).reverse_complement())
                                rf_seq_mutated = str(Seq(rf_seq_mutated, generic_dna).reverse_complement())

                            rf_seq_raw_aa = str(SeqRecord(Seq(rf_seq_raw)).seq.translate(table=transl_table))
                            rf_seq_mutated_aa = str(SeqRecord(Seq(rf_seq_mutated)).seq.translate(table=transl_table))

                            mutation_type_term = ''
                            if rf_seq_raw_aa == rf_seq_mutated_aa:
                                mutation_type_term = 'Silent'
                            elif rf_seq_mutated_aa == '*':
                                mutation_type_term = 'Nonsense'
                            else:
                                mutation_type_term = 'Missense'
                            aa_mutation = '%s(%s)->%s(%s)' % (rf_seq_raw, rf_seq_raw_aa, rf_seq_mutated, rf_seq_mutated_aa)

                            # print out
                            output_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (each_snv.strip().split('\t')[0], each_snv_parallel, location, ORF_strand_dict[each_gene], each_gene, aa_mutation, mutation_type_term, current_COG_id2, current_COG_function2))
                            #output_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (each_snv.strip(), each_snv_parallel, location, ORF_strand_dict[each_gene], each_gene, aa_mutation, mutation_type_term, current_COG_id2, current_COG_function2))

    output_handle.close()

    # get mutation_cate_summary
    #mutation_cate_dict = get_mutation_cate_summary(output_mutated_genes, output_mutated_genes_cate)

    # report
    #print('The number of affected genes for %s: %s' % (SNV_matrix_cdc, len(mutation_cate_dict)))


    ################################################### export sequences ###################################################

    # # get the sequence of affected genes
    # output_seq_nc_handle = open(output_seq_nc, 'w')
    # for each_gene_nc in SeqIO.parse(combined_ref_ffn, 'fasta'):
    #     if each_gene_nc.id in mutation_cate_dict:
    #         SeqIO.write(each_gene_nc, output_seq_nc_handle, 'fasta')
    # output_seq_nc_handle.close()
    #
    # output_seq_aa_handle = open(output_seq_aa, 'w')
    # for each_gene_aa in SeqIO.parse(combined_ref_faa, 'fasta'):
    #     if each_gene_aa.id in mutation_cate_dict:
    #         SeqIO.write(each_gene_aa, output_seq_aa_handle, 'fasta')
    # output_seq_aa_handle.close()



    ############################################### combine mutation effect ################################################

    # remove tmp file
    # os.system('rm %s' % output_mutated_genes_cate)
    # os.system('rm %s' % output_seq_nc)
    # os.system('rm %s' % output_seq_aa)




