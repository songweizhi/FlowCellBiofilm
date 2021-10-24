import os
import shutil
import argparse
import numpy as np
from datetime import datetime
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from scipy import stats
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import os
import shutil


def force_create_folder(folder_to_create):
    if os.path.isdir(folder_to_create):
        shutil.rmtree(folder_to_create, ignore_errors=True)
        if os.path.isdir(folder_to_create):
            shutil.rmtree(folder_to_create, ignore_errors=True)
            if os.path.isdir(folder_to_create):
                shutil.rmtree(folder_to_create, ignore_errors=True)
                if os.path.isdir(folder_to_create):
                    shutil.rmtree(folder_to_create, ignore_errors=True)
    os.mkdir(folder_to_create)


def turn_element_to_str(list_in):

    list_out = []
    for each_element in list_in:
        each_element_str = str(each_element)
        list_out.append(each_element_str)

    return list_out


def merge_two_dict(dict_a, dict_b):
    dict_c = dict_a.copy()
    dict_c.update(dict_b)
    return dict_c


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '..'

    # separate file basename and extension
    file_basename, file_ext = os.path.splitext(file_name)

    return file_path, file_basename, file_ext


def remove_l2_elements_from_l1(l1, l2):

    l1_new = []
    for each in l1:
        if each not in l2:
            l1_new.append(each)
    return l1_new


def unique_list_elements(list_input):

    list_output = []
    for each_element in list_input:
        if each_element not in list_output:
            list_output.append(each_element)

    return list_output


def write_out(list, op_file):

    op_file_handle = open(op_file, 'w')
    for each in list:
        op_file_handle.write('%s\n' % each)
    op_file_handle.close()


def combined_continuous_deletions(pwd_QC_txt_cd, pwd_QC_txt_cd_combined):

    current_sample = ''
    current_seq = ''
    current_pos_start = 0
    previous_snv_pos = 0
    deleted_ncs = ''
    deepSNV_output_cdc_handle = open(pwd_QC_txt_cd_combined, 'w')
    for each_snv in open(pwd_QC_txt_cd):

        each_snv_split = each_snv.strip().split(',')
        each_snv_sample =   each_snv_split[0]
        each_snv_seq =      each_snv_split[1]
        each_snv_pos =      int(each_snv_split[2])
        each_snv_pos_wt =   each_snv_split[3]

        if (current_sample == '') and (current_seq == '') and (current_pos_start == 0):
            current_sample = each_snv_sample
            current_seq = each_snv_seq
            current_pos_start = each_snv_pos
            previous_snv_pos = each_snv_pos
            deleted_ncs = each_snv_pos_wt
        elif (current_sample == each_snv_sample) and (current_seq == each_snv_seq) and (each_snv_pos == previous_snv_pos + 1):
            previous_snv_pos = each_snv_pos
            deleted_ncs += each_snv_pos_wt
        elif (current_sample != each_snv_sample) or (current_seq != each_snv_seq) or (each_snv_pos != previous_snv_pos + 1):
            for_out = '%s|%s|%s-%s|%s|%s\n' % (current_sample, current_seq, current_pos_start, previous_snv_pos, deleted_ncs, '-'*len(deleted_ncs))
            deepSNV_output_cdc_handle.write(for_out)

            current_sample = each_snv_sample
            current_seq = each_snv_seq
            current_pos_start = each_snv_pos
            previous_snv_pos = each_snv_pos
            deleted_ncs = each_snv_pos_wt

    for_out = '%s|%s|%s-%s|%s|%s\n' % (current_sample, current_seq, current_pos_start, previous_snv_pos, deleted_ncs, '-'*len(deleted_ncs))
    deepSNV_output_cdc_handle.write(for_out)

    deepSNV_output_cdc_handle.close()


def separate_continuous_deletion(file_in, file_out_cd, file_out_ncd):

    # write out separately
    file_out_cd_handle = open(file_out_cd, 'w')
    file_out_ncd_handle = open(file_out_ncd, 'w')

    current_sample_id = ''
    current_seq = ''
    current_pos = 0
    current_line = ''
    current_line_wrote = 0
    for each_snv in open(file_in):
        if not each_snv.startswith('sample'):

            each_snv_sample_id = each_snv.strip().split(',')[0]
            each_snv_seq = each_snv.strip().split(',')[1]
            each_snv_pos = int(each_snv.strip().split(',')[2])
            each_snv_var = each_snv.strip().split(',')[4]

            # if not deletion, write out previous deletion (if have) and current SNV
            if each_snv_var != '-':

                if current_line == '':
                    file_out_ncd_handle.write(each_snv)
                elif current_line != '':
                    if current_line_wrote == 0:
                        file_out_ncd_handle.write(current_line)
                    file_out_ncd_handle.write(each_snv)

                # reset
                current_sample_id = ''
                current_seq = ''
                current_pos = 0
                current_line = ''
                current_line_wrote = 0

            # if it is a deletion
            elif each_snv_var == '-':

                if current_line == '':
                    current_sample_id = each_snv_sample_id
                    current_seq = each_snv_seq
                    current_pos = each_snv_pos
                    current_line = each_snv
                    current_line_wrote = 0

                elif current_line != '':
                    if (each_snv_sample_id == current_sample_id) and (each_snv_seq == current_seq) and (each_snv_pos == current_pos + 1):
                        if current_line_wrote == 0:
                            file_out_cd_handle.write(current_line)
                        file_out_cd_handle.write(each_snv)
                        current_line = each_snv
                        current_line_wrote = 1
                        current_pos = each_snv_pos

                    else:
                        if current_line_wrote == 0:
                            file_out_ncd_handle.write(current_line)

                        current_sample_id = each_snv_sample_id
                        current_seq = each_snv_seq
                        current_pos = each_snv_pos
                        current_line = each_snv
                        current_line_wrote = 0

    file_out_cd_handle.close()
    file_out_ncd_handle.close()


def get_rf_pos_list(gene_len):

    rf_pos_list = []
    n = 0
    while n < (gene_len / 3):
        start_pos = 1 + 3 * n
        rf_pos_list.append([start_pos, start_pos + 1, start_pos + 2])
        n += 1

    return rf_pos_list


def get_categorical_scatter_plot(num_list, label_list, label_rotation, output_plot):

    num_index = list(range(1, len(num_list) + 1))

    plt.scatter(num_index, num_list)
    plt.xticks(num_index, label_list, rotation=label_rotation)
    plt.margins(0.2)
    plt.subplots_adjust(bottom=0.15)
    plt.tight_layout()
    plt.savefig(output_plot, dpi=300)
    plt.close()


def get_matrix_freq(SNV_file_in, SNV_file_out_frequency_210, SNV_file_out_frequency_D2):

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

    # get SNV occurrence matrix for 2.10 and D2
    SNV_file_out_frequency_210_handle = open(SNV_file_out_frequency_210, 'w')
    SNV_file_out_frequency_D2_handle = open(SNV_file_out_frequency_D2, 'w')

    SNV_file_out_frequency_210_handle.write('\t%s\n' % '\t'.join(matrix_header_210))
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

            for_write_frequency_210 = '%s\t%s\n' % (each_uniq_snv, '\t'.join(occurrence_list_frequency_210))
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

            for_write_frequency_D2 = '%s\t%s\n' % (each_uniq_snv, '\t'.join(occurrence_list_frequency_D2))
            SNV_file_out_frequency_D2_handle.write(for_write_frequency_D2)

    SNV_file_out_frequency_210_handle.close()
    SNV_file_out_frequency_D2_handle.close()


def get_matrix(SNV_file_in, SNV_file_out_frequency_210, SNV_file_out_existence_210, SNV_file_out_frequency_D2, SNV_file_out_existence_D2):

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

    # get SNV occurrence matrix for 2.10 and D2
    SNV_file_out_frequency_210_handle = open(SNV_file_out_frequency_210, 'w')
    SNV_file_out_frequency_D2_handle = open(SNV_file_out_frequency_D2, 'w')
    SNV_file_out_existence_210_handle = open(SNV_file_out_existence_210, 'w')
    SNV_file_out_existence_D2_handle = open(SNV_file_out_existence_D2, 'w')

    SNV_file_out_frequency_210_handle.write('\t%s\n' % '\t'.join(matrix_header_210))
    SNV_file_out_frequency_D2_handle.write('\t%s\n' % '\t'.join(matrix_header_D2))
    SNV_file_out_existence_210_handle.write('\t%s\n' % '\t'.join(matrix_header_210))
    SNV_file_out_existence_D2_handle.write('\t%s\n' % '\t'.join(matrix_header_D2))


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

            for_write_frequency_D2 = '%s\t%s\n' % (each_uniq_snv, '\t'.join(occurrence_list_frequency_D2))
            for_write_existence_D2 = '%s\t%s\n' % (each_uniq_snv, '\t'.join(occurrence_list_existence_D2))

            SNV_file_out_frequency_D2_handle.write(for_write_frequency_D2)
            SNV_file_out_existence_D2_handle.write(for_write_existence_D2)

    SNV_file_out_frequency_210_handle.close()
    SNV_file_out_frequency_D2_handle.close()
    SNV_file_out_existence_210_handle.close()
    SNV_file_out_existence_D2_handle.close()


def check_occurrence(list_in):

    current_occurrence = 0
    for each in list_in:
        if each != '0':
            current_occurrence = 1

    return current_occurrence


def check_parallel(matrix_file):

    parallel_dict = {}
    for snv in open(matrix_file):

        if not snv.startswith('\t'):
            snv_split = snv.strip().split('\t')
            snv_id = snv_split[0]
            mono_A_occurrence = check_occurrence([snv_split[1], snv_split[2], snv_split[3], snv_split[4]])
            mono_B_occurrence = check_occurrence([snv_split[5], snv_split[6], snv_split[7], snv_split[8]])
            mono_C_occurrence = check_occurrence([snv_split[9], snv_split[10], snv_split[11], snv_split[12]])
            co_A_occurrence =   check_occurrence([snv_split[13], snv_split[14], snv_split[15], snv_split[16]])
            co_B_occurrence =   check_occurrence([snv_split[17], snv_split[18], snv_split[19], snv_split[20]])
            co_C_occurrence =   check_occurrence([snv_split[21], snv_split[22], snv_split[23], snv_split[24]])

            mono_occurrence = [mono_A_occurrence, mono_B_occurrence, mono_C_occurrence]
            co_occurrence = [co_A_occurrence, co_B_occurrence, co_C_occurrence]

            # get parallel profile
            parallel_profile = 'NA'
            if (mono_occurrence.count(1) >= 1) and (co_occurrence.count(1) >= 1):
                parallel_profile = ('PMC(%s_%s)' % (''.join(turn_element_to_str(mono_occurrence)), ''.join(turn_element_to_str(co_occurrence))))

            elif (mono_occurrence.count(1) > 1) and (co_occurrence.count(1) == 0):
                parallel_profile = ('PM(%s_%s)' % (''.join(turn_element_to_str(mono_occurrence)), ''.join(turn_element_to_str(co_occurrence))))

            elif (mono_occurrence.count(1) == 0) and (co_occurrence.count(1) > 1):
                parallel_profile = ('PC(%s_%s)' % (''.join(turn_element_to_str(mono_occurrence)), ''.join(turn_element_to_str(co_occurrence))))

            elif (mono_occurrence.count(1) == 1) and (co_occurrence.count(1) == 0):
                parallel_profile = ('NPM(%s_%s)' % (''.join(turn_element_to_str(mono_occurrence)), ''.join(turn_element_to_str(co_occurrence))))

            elif (mono_occurrence.count(1) == 0) and (co_occurrence.count(1) == 1):
                parallel_profile = ('NPC(%s_%s)' % (''.join(turn_element_to_str(mono_occurrence)), ''.join(turn_element_to_str(co_occurrence))))

            parallel_dict[snv_id] = parallel_profile

    return parallel_dict


def get_affect_backup(frequency_list, plot_effect, plot_folder, plot_filename):

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
    elif (slope > 0) and (max_frequency >= 5):
        affect = 'Positive'
    elif (slope < 0) and (max_frequency <= 5):
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


def get_mutation_cate_summary(file_in, file_out):

    mutation_cate_dict = {}
    for each_snv in open(file_in):
        each_snv_split = each_snv.strip().split('\t')
        mutated_gene = each_snv_split[4]
        mutation_cate = each_snv_split[6]

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


def scatter_plotter(frequency_list, label_name, plot_title, plot_filename):

    frequency_list_percentage = [i * 100 for i in frequency_list]

    num_index = list(range(1, len(frequency_list_percentage) + 1))
    plt.scatter(num_index, frequency_list_percentage)
    plt.xticks(num_index, label_name, rotation=0)

    for i, txt in enumerate(frequency_list_percentage):
        plt.annotate(txt, (num_index[i], frequency_list_percentage[i]))

    plt.margins(0.2)
    plt.subplots_adjust(bottom=0.15)
    plt.xlabel('Time point')
    plt.ylabel('Frequency (%)')
    plt.title(plot_title, fontsize=5)
    plt.tight_layout()
    plt.savefig(plot_filename, dpi=300)
    plt.close()


def plot_freq(SNV_matrix_cdc, plot_folder, SNV_to_gene_dict, SNV_to_ko_id_dict, SNV_to_ko_desc_dict):

    sample_dict = {'1': 'Mono210_A', '5': 'Mono210_B', '9': 'Mono210_C',
                   '2': 'MonoD2_A', '6': 'MonoD2_B', '10': 'MonoD2_C',
                   '4': 'Coculture_A', '8': 'Coculture_B', '12': 'Coculture_C'}

    # get mutation affect
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
            affected_gene = SNV_to_gene_dict[snv_id]
            affected_gene_ko_id = SNV_to_ko_id_dict[snv_id]
            affected_gene_ko_desc = SNV_to_ko_desc_dict[snv_id]
            overall_title = '%s__%s__%s__%s' % (snv_id, affected_gene, affected_gene_ko_id, affected_gene_ko_desc)
            num_index = [1, 2, 3, 4]
            label_name = ['D9', 'D18', 'D27', 'D42']
            combined_plot = '%s/%s.png' % (plot_folder, snv_id)

            plt.figure()

            plt.subplot(231)
            subplot_231_frequency_list = [float(each_snv_split[1]) * 100, float(each_snv_split[2]) * 100, float(each_snv_split[3]) * 100, float(each_snv_split[4]) * 100]
            if (each_snv_split[1] != '0') or (each_snv_split[2] != '0') or (each_snv_split[3] != '0') or (each_snv_split[4] != '0'):
                plt.scatter(num_index, subplot_231_frequency_list, s=8)
                for i, txt in enumerate(subplot_231_frequency_list):
                    plt.annotate(txt, (num_index[i], subplot_231_frequency_list[i]), fontsize=5)
                plt.xticks(num_index, label_name, rotation=0, fontsize=6)
                plt.yticks(fontsize=6)
            else:
                plt.xticks([])
                plt.yticks([])
            plt.ylabel('Mono-culture', fontsize=8)
            plt.title('A', fontsize=8)

            plt.subplot(232)
            subplot_232_frequency_list = [float(each_snv_split[5]) * 100, float(each_snv_split[6]) * 100, float(each_snv_split[7]) * 100, float(each_snv_split[8]) * 100]
            if (each_snv_split[5] != '0') or (each_snv_split[6] != '0') or (each_snv_split[7] != '0') or (each_snv_split[8] != '0'):
                plt.scatter(num_index, subplot_232_frequency_list, s=8)
                for i, txt in enumerate(subplot_232_frequency_list):
                    plt.annotate(txt, (num_index[i], subplot_232_frequency_list[i]), fontsize=5)
                plt.xticks(num_index, label_name, rotation=0, fontsize=6)
                plt.yticks(fontsize=6)
            else:
                plt.xticks([])
                plt.yticks([])
            plt.title('B', fontsize=8)

            plt.subplot(233)
            subplot_233_frequency_list = [float(each_snv_split[9]) * 100, float(each_snv_split[10]) * 100, float(each_snv_split[11]) * 100, float(each_snv_split[12]) * 100]
            if (each_snv_split[9] != '0') or (each_snv_split[10] != '0') or (each_snv_split[11] != '0') or (each_snv_split[12] != '0'):
                plt.scatter(num_index, subplot_233_frequency_list, s=8)
                for i, txt in enumerate(subplot_233_frequency_list):
                    plt.annotate(txt, (num_index[i], subplot_233_frequency_list[i]), fontsize=5)
                plt.xticks(num_index, label_name, rotation=0, fontsize=6)
                plt.yticks(fontsize=6)
            else:
                plt.xticks([])
                plt.yticks([])
            plt.title('C', fontsize=8)

            plt.subplot(234)
            subplot_234_frequency_list = [float(each_snv_split[13]) * 100, float(each_snv_split[14]) * 100, float(each_snv_split[15]) * 100, float(each_snv_split[16]) * 100]
            if (each_snv_split[13] != '0') or (each_snv_split[14] != '0') or (each_snv_split[15] != '0') or (each_snv_split[16] != '0'):
                plt.scatter(num_index, subplot_234_frequency_list, s=8)
                for i, txt in enumerate(subplot_234_frequency_list):
                    plt.annotate(txt, (num_index[i], subplot_234_frequency_list[i]), fontsize=5)
                plt.xticks(num_index, label_name, rotation=0, fontsize=6)
                plt.yticks(fontsize=6)
            else:
                plt.xticks([])
                plt.yticks([])
            plt.ylabel('Co-culture', fontsize=8)

            plt.subplot(235)
            subplot_235_frequency_list = [float(each_snv_split[17]) * 100, float(each_snv_split[18]) * 100, float(each_snv_split[19]) * 100, float(each_snv_split[20]) * 100]
            if (each_snv_split[17] != '0') or (each_snv_split[18] != '0') or (each_snv_split[19] != '0') or (each_snv_split[20] != '0'):
                plt.scatter(num_index, subplot_235_frequency_list, s=8)
                for i, txt in enumerate(subplot_235_frequency_list):
                    plt.annotate(txt, (num_index[i], subplot_235_frequency_list[i]), fontsize=5)
                plt.xticks(num_index, label_name, rotation=0, fontsize=6)
                plt.yticks(fontsize=6)
            else:
                plt.xticks([])
                plt.yticks([])

            plt.subplot(236)
            subplot_236_frequency_list = [float(each_snv_split[21]) * 100, float(each_snv_split[22]) * 100, float(each_snv_split[23]) * 100, float(each_snv_split[24]) * 100]
            if (each_snv_split[21] != '0') or (each_snv_split[22] != '0') or (each_snv_split[23] != '0') or (each_snv_split[24] != '0'):
                plt.scatter(num_index, subplot_236_frequency_list, s=8)
                for i, txt in enumerate(subplot_236_frequency_list):
                    plt.annotate(txt, (num_index[i], subplot_236_frequency_list[i]), fontsize=5)
                plt.xticks(num_index, label_name, rotation=0, fontsize=6)
                plt.yticks(fontsize=6)
            else:
                plt.xticks([])
                plt.yticks([])

            plt.suptitle(overall_title, fontsize=5)
            plt.savefig(combined_plot, dpi=300)
            plt.close()


############################################### input file and parameters ##############################################

parser = argparse.ArgumentParser(description='', add_help=False)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')

optional.add_argument('-h', action='help', help='Show this help message and exit')
required.add_argument('-snv_qc',        dest='SNV_QC',      nargs='?', required=True,   type=str, help='deepSNV QC file')
required.add_argument('-min_both',      dest='MIN_BOTH',    nargs='?', required=True,   type=int, help='The minimum number of reads harboring SNV')
required.add_argument('-min_each',      dest='MIN_EACH',    nargs='?', required=True,   type=int, help='The minimum number of reads harboring SNV at each direction')
required.add_argument('-strand_bias',   dest='STRAND_BIAS', nargs='?', required=True,   type=int, help='strand_bias cutoff')
required.add_argument('-depth_diff',    dest='DEPTH_DIFF',  nargs='?', required=True,   type=int, help='depth difference cutoff')
required.add_argument('-deplen',        dest='DEPLEN',      nargs='?', required=True,   type=int, help='flanking length for mean depth calculation')
required.add_argument('-sep_plot',      dest='SEP_PLOT',    nargs='?', required=False,  type=int, help='separate depth plot with provide cutoff')
required.add_argument('-out',           dest='OUT',         nargs='?', required=True,             help='output dir')

args = vars(parser.parse_args())
SNV_quality_file = args['SNV_QC']
min_var_reads_num = args['MIN_BOTH']
min_at_each_direction = args['MIN_EACH']
strand_bias_cutoff = args['STRAND_BIAS']
depth_difference_cutoff = args['DEPTH_DIFF']
mean_depth_len = args['DEPLEN']
separate_plot = args['SEP_PLOT']

transl_table = 11
ending_len_cutoff = 50000
exculding_regions = ''
plot_snv_freq = False
permanova = 1
pwd_plot_folder = 'SNV_depth_plot_output_f50000bp_1000mer_dl2000bp'
pwd_plot_folder = 'SNV_depth_plot'

output_dir = 'OneStep_MinBoth_%s_MinEach_%s_StrandBias_%s_DepthDiff_%s' % (min_var_reads_num,
                                                                           min_at_each_direction,
                                                                           strand_bias_cutoff,
                                                                           depth_difference_cutoff)

################################################ prepare files and list ################################################

reference_seq_folder                        = '/Users/songweizhi/FC_Biofilm/reference_files'
combined_ref_fasta                          = '%s/combined_ref.fasta'                       % reference_seq_folder
combined_ref_gff                            = '%s/combined_ref.gff'                         % reference_seq_folder
combined_ref_ffn                            = '%s/combined_ref.ffn'                         % reference_seq_folder
combined_ref_faa                            = '%s/combined_ref.faa'                         % reference_seq_folder
annotation_file_KEGG                        = '%s/combined_ref_ko_assignment_ABCD.txt'      % reference_seq_folder
annotation_file_COG_for_ko_unknown_genes    = '%s/function_unknown_genes_query_to_cog.txt'  % reference_seq_folder

matrix_header_210 = ['1D9', '1D18', '1D27', '1D42', '5D9', '5D18', '5D27', '5D42', '9D9', '9D18', '9D27', '9D42', '4D9', '4D18', '4D27', '4D42', '8D9', '8D18', '8D27', '8D42', '12D9', '12D18', '12D27', '12D42']
matrix_header_D2  = ['2D9', '2D18', '2D27', '2D42', '6D9', '6D18', '6D27', '6D42', '10D9', '10D18', '10D27', '10D42', '4D9', '4D18', '4D27', '4D42', '8D9', '8D18', '8D27', '8D42', '12D9', '12D18', '12D27', '12D42']

force_create_folder(output_dir)

# parallel profile
# PMC   parallel in monoculture and coculture
# PM    parallel in monoculture
# PC    parallel in coculture
# NPM   non-parallel and monoculture only
# NPC   non-parallel and coculture only

################################################ define output filename ################################################

pwd_QC_txt_sorted =                     '%s/SNV_QC_sorted.txt'                                      % output_dir
pwd_QC_txt_cd =                         '%s/SNV_QC_cd.txt'                                          % output_dir
pwd_QC_txt_cd_combined =                '%s/SNV_QC_cd_combined.txt'                                 % output_dir
pwd_QC_txt_ncd =                        '%s/SNV_QC_ncd.txt'                                         % output_dir
qualified_SNVs_even_flk_depth_file =    '%s/SNV_QC_ncd_even_flk_depth.txt'                          % output_dir
qualified_SNVs_diff_flk_depth_file =    '%s/SNV_QC_ncd_diff_flk_depth.txt'                          % output_dir
pwd_plot_folder_sep =                   '%s/SNV_depth_plot_sep'                                     % output_dir
pwd_plot_folder_endings =               '%s/%s_endings_%sbp'                                        % (pwd_plot_folder_sep, pwd_plot_folder, ending_len_cutoff)
pwd_plot_folder_even =                  '%s/%s_even_%s'                                             % (pwd_plot_folder_sep, pwd_plot_folder, depth_difference_cutoff)
pwd_plot_folder_diff =                  '%s/%s_diff_%s'                                             % (pwd_plot_folder_sep, pwd_plot_folder, depth_difference_cutoff)
pwd_plot_folder_unqualified =           '%s/%s_unqualified'                                         % (pwd_plot_folder_sep, pwd_plot_folder)
SNV_210_mono_uniq_file =                '%s/SNV_QC_ncd_even_flk_depth_210_monoculture_uniq.txt'     % output_dir
SNV_210_co_uniq_file =                  '%s/SNV_QC_ncd_even_flk_depth_210_coculture_uniq.txt'       % output_dir
SNV_210_concurrence_file =              '%s/SNV_QC_ncd_even_flk_depth_210_concurrence.txt'          % output_dir
SNV_210_matrix_file =                   '%s/SNV_QC_ncd_even_flk_depth_210_matrix.txt'               % output_dir
SNV_210_matrix_file_existence =         '%s/SNV_QC_ncd_even_flk_depth_210_matrix_existence.txt'     % output_dir
SNV_D2_mono_uniq_file =                 '%s/SNV_QC_ncd_even_flk_depth_D2_monoculture_uniq.txt'      % output_dir
SNV_D2_co_uniq_file =                   '%s/SNV_QC_ncd_even_flk_depth_D2_coculture_uniq.txt'        % output_dir
SNV_D2_concurrence_file =               '%s/SNV_QC_ncd_even_flk_depth_D2_concurrence.txt'           % output_dir
SNV_D2_matrix_file =                    '%s/SNV_QC_ncd_even_flk_depth_D2_matrix.txt'                % output_dir
SNV_D2_matrix_file_existence =          '%s/SNV_QC_ncd_even_flk_depth_D2_matrix_existence.txt'      % output_dir
mutation_effect_plot_folder =           '%s/SNV_QC_ncd_even_flk_depth_frequency_plot'               % output_dir
gene_level_output =                     '%s/SNV_QC_ncd_even_flk_depth_affected_genes.txt'           % output_dir


########################################## read function annotation into dict ##########################################

# store annotation results in dicts
gene_KEGG_id_dict = {}
gene_KEGG_function_dict = {}
for each_snv3 in open(annotation_file_KEGG):
    if not each_snv3.startswith('Gene_id'):
        each_snv3_split = each_snv3.strip().split('\t')
        gene_id = each_snv3_split[0]
        if len(each_snv3_split) > 1:
            ko_id = each_snv3_split[4][2:]
            ko_function = each_snv3_split[8]
            gene_KEGG_id_dict[gene_id] = ko_id
            gene_KEGG_function_dict[gene_id] = ko_function
        else:
            gene_KEGG_id_dict[gene_id] = 'NA'
            gene_KEGG_function_dict[gene_id] = 'NA'


for each_snv3 in open(annotation_file_COG_for_ko_unknown_genes):
    if not each_snv3.startswith('Query'):
        each_snv3_split = each_snv3.strip().split('\t')
        gene_id = each_snv3_split[0]
        if len(each_snv3_split) > 1:
            cog_id = each_snv3_split[1]
            cog_function = each_snv3_split[3]
            gene_KEGG_id_dict[gene_id] = cog_id
            gene_KEGG_function_dict[gene_id] = cog_function

############################################## remove continuous deletion ##############################################

# sort QC file
os.system('cat %s | sort > %s' % (SNV_quality_file, pwd_QC_txt_sorted))

# separate continuous deletion
separate_continuous_deletion(pwd_QC_txt_sorted, pwd_QC_txt_cd, pwd_QC_txt_ncd)

# combine continuous deletion
combined_continuous_deletions(pwd_QC_txt_cd, pwd_QC_txt_cd_combined)

# delete tmp file
os.system('rm %s' % pwd_QC_txt_sorted)
os.system('rm %s' % pwd_QC_txt_cd)


##################################################### Parse SNV QC #####################################################

# get reference length dict
ref_len_dict = {}
for each_ref in SeqIO.parse(combined_ref_fasta, 'fasta'):
    ref_len_dict[each_ref.id] = len(each_ref.seq)

qualified_SNVs_even_flanking_depth_file_handle = open(qualified_SNVs_even_flk_depth_file, 'w')
qualified_SNVs_diff_flanking_depth_file_handle = open(qualified_SNVs_diff_flk_depth_file, 'w')

SNV_at_endings = []
n_tst_total_unqualified = []
n_tst_each_unqualified = []
strand_bias_unqualified = []
qualified_SNVs = []
qualified_SNVs_even_flanking_depth = []
qualified_SNVs_diff_flanking_depth = []

plot_prefix_even = []
plot_prefix_diff = []
plot_prefix_unqualified = []

total_num = 0
for each_snv in open(pwd_QC_txt_ncd):

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
        plot_prefix = '%s_%s_%s_%s_%s' % (each_snv_chr, each_snv_pos, each_snv_ref, each_snv_var, each_snv_sample)

        # get SNVs located at endings
        if (each_snv_pos <= ending_len_cutoff) or (ref_len_dict[each_snv_chr] - each_snv_pos <= ending_len_cutoff):
            SNV_at_endings.append(plot_prefix)

        else:
            # get SNVs with unqualified number of reads harbouring SNV
            if each_snv_n_tst_b < min_var_reads_num:
                n_tst_total_unqualified.append(each_snv_key_with_sample_id)

                if plot_prefix not in plot_prefix_unqualified:
                    plot_prefix_unqualified.append(plot_prefix)

            # get SNVs with unqualified number of reads in each direction
            if (each_snv_n_tst_fw < min_at_each_direction) or (each_snv_n_tst_bw < min_at_each_direction):
                n_tst_each_unqualified.append(each_snv_key_with_sample_id)

                if plot_prefix not in plot_prefix_unqualified:
                    plot_prefix_unqualified.append(plot_prefix)

            # get SNVs with unqualified strand bias
            if each_snv_strand_bias > strand_bias_cutoff:
                strand_bias_unqualified.append(each_snv_key_with_sample_id)

                if plot_prefix not in plot_prefix_unqualified:
                    plot_prefix_unqualified.append(plot_prefix)

            # get SNVs with flanking depth difference higher than defined cutoff
            if (each_snv_n_tst_b >= min_var_reads_num) and (not ((each_snv_n_tst_fw < min_at_each_direction) or (each_snv_n_tst_bw < min_at_each_direction))) and (each_snv_strand_bias <= strand_bias_cutoff):
                qualified_SNVs.append(each_snv_key_with_sample_id)

            # get qualified SNV with similar flanking depth
            if (each_snv_n_tst_b >= min_var_reads_num) and (each_snv_n_tst_fw >= min_at_each_direction) and (each_snv_n_tst_bw >= min_at_each_direction) and (each_snv_strand_bias <= strand_bias_cutoff) and (each_snv_mean_depth_diff <= depth_difference_cutoff):
                qualified_SNVs_even_flanking_depth.append(each_snv_key_with_sample_id)
                plot_prefix_even.append(plot_prefix)
                qualified_SNVs_even_flanking_depth_file_handle.write(each_snv)

            # get qualified SNV with different flanking depth
            if (each_snv_n_tst_b >= min_var_reads_num) and (each_snv_n_tst_fw >= min_at_each_direction) and (each_snv_n_tst_bw >= min_at_each_direction) and (each_snv_strand_bias <= strand_bias_cutoff) and (each_snv_mean_depth_diff > depth_difference_cutoff):
                qualified_SNVs_diff_flanking_depth.append(each_snv_key_with_sample_id)
                plot_prefix_diff.append(plot_prefix)
                qualified_SNVs_diff_flanking_depth_file_handle.write(each_snv)

        total_num += 1

qualified_SNVs_even_flanking_depth_file_handle.close()
qualified_SNVs_diff_flanking_depth_file_handle.close()


# For report
print('\n########################################### Report 1 ###########################################\n')

print(datetime.now().strftime('[%Y-%m-%d %H:%M:%S]') + ' The total number of detected SNVs: %s'                                 % total_num)
print(datetime.now().strftime('[%Y-%m-%d %H:%M:%S]') + ' The number of SNVs located at ends (within %sbp): %s'                  % (ending_len_cutoff, len(SNV_at_endings)))
print(datetime.now().strftime('[%Y-%m-%d %H:%M:%S]') + ' The number of SNVs with less than %s reads harboring it: %s'           % (min_var_reads_num, len(n_tst_total_unqualified)))
print(datetime.now().strftime('[%Y-%m-%d %H:%M:%S]') + ' The number of SNVs with reads only from one direction: %s'             % len(n_tst_each_unqualified))
print(datetime.now().strftime('[%Y-%m-%d %H:%M:%S]') + ' The number of SNVs with strand bias higher than %s: %s'                % (strand_bias_cutoff, len(strand_bias_unqualified)))
print(datetime.now().strftime('[%Y-%m-%d %H:%M:%S]') + ' The number of unqualified SNVs: %s'                                    % len(plot_prefix_unqualified))
print(datetime.now().strftime('[%Y-%m-%d %H:%M:%S]') + ' The number of qualified SNVs: %s'                                      % len(qualified_SNVs))
print(datetime.now().strftime('[%Y-%m-%d %H:%M:%S]') + ' The number of qualified SNVs with flk depth (%sbp) diff >  %s%s: %s'   % (mean_depth_len, depth_difference_cutoff, '%', len(qualified_SNVs_diff_flanking_depth)))
print(datetime.now().strftime('[%Y-%m-%d %H:%M:%S]') + ' The number of qualified SNVs with flk depth (%sbp) diff <= %s%s: %s'   % (mean_depth_len, depth_difference_cutoff, '%', len(qualified_SNVs_even_flanking_depth)))


############################# separate depth plots according to provided difference cutoff #############################
separate_plot = 1
# separate plot files
if separate_plot == 1:

    # prepare folder
    if os.path.isdir(pwd_plot_folder_sep):
        shutil.rmtree(pwd_plot_folder_sep)
        if os.path.isdir(pwd_plot_folder_sep):
            shutil.rmtree(pwd_plot_folder_sep)
        os.mkdir(pwd_plot_folder_sep)
    else:
        os.mkdir(pwd_plot_folder_sep)

    for folder in [pwd_plot_folder_endings, pwd_plot_folder_even, pwd_plot_folder_diff, pwd_plot_folder_unqualified]:
        os.mkdir(folder)

    # get plots for SNV located at ending regions
    for each_even in SNV_at_endings:
        pwd_plot = '%s/%s*' % (pwd_plot_folder, each_even)
        cmd = 'cp %s %s/' % (pwd_plot, pwd_plot_folder_endings)
        os.system(cmd)

    # get plots with similar flanking depth
    for each_even in plot_prefix_even:
        pwd_plot = '%s/%s*' % (pwd_plot_folder, each_even)
        cmd = 'cp %s %s/' % (pwd_plot, pwd_plot_folder_even)
        os.system(cmd)

    # get plots with similar flanking depth
    for each_diff in plot_prefix_diff:
        pwd_plot = '%s/%s*' % (pwd_plot_folder, each_diff)
        cmd = 'cp %s %s/' % (pwd_plot, pwd_plot_folder_diff)
        os.system(cmd)

    # get plots for unqualified SNVs
    for each_unqualified in plot_prefix_unqualified:
        pwd_plot = '%s/%s*' % (pwd_plot_folder, each_unqualified)
        cmd = 'cp %s %s/' % (pwd_plot, pwd_plot_folder_unqualified)
        os.system(cmd)


####################################### compare between Monoculture and Coculture ######################################

mono_210_SNV_list = []
mono_D2_SNV_list = []
co_210_SNV_list = []
co_D2_SNV_list = []
total_SNV = 0
for each_SNV in open(qualified_SNVs_even_flk_depth_file):
    each_SNV_split = each_SNV.strip().split(',')
    treatment_id = each_SNV_split[0].split('D')[0]
    strain = each_SNV_split[1].split('_')[0]
    SNV_id = '%s|%s|%s|%s' % (each_SNV_split[1], each_SNV_split[2], each_SNV_split[3], each_SNV_split[4])

    if (strain == '2.10') and (treatment_id in ['1', '5', '9']):
        mono_210_SNV_list.append(SNV_id)
    elif (strain == '2.10') and (treatment_id in ['4', '8', '12']):
        co_210_SNV_list.append(SNV_id)
    elif (strain == 'D2') and (treatment_id in ['2', '6', '10']):
        mono_D2_SNV_list.append(SNV_id)
    elif (strain == 'D2') and (treatment_id in ['4', '8', '12']):
        co_D2_SNV_list.append(SNV_id)

    total_SNV += 1

# uniq list element
mono_210_SNV_list = unique_list_elements(mono_210_SNV_list)
co_210_SNV_list = unique_list_elements(co_210_SNV_list)
mono_D2_SNV_list = unique_list_elements(mono_D2_SNV_list)
co_D2_SNV_list = unique_list_elements(co_D2_SNV_list)

# 210
SNV_210_concurrence = set(mono_210_SNV_list).intersection(co_210_SNV_list)
SNV_210_mono_uniq = remove_l2_elements_from_l1(mono_210_SNV_list, SNV_210_concurrence)
SNV_210_co_uniq = remove_l2_elements_from_l1(co_210_SNV_list, SNV_210_concurrence)

# D2
SNV_D2_concurrence = set(mono_D2_SNV_list).intersection(co_D2_SNV_list)
SNV_D2_mono_uniq = remove_l2_elements_from_l1(mono_D2_SNV_list, SNV_D2_concurrence)
SNV_D2_co_uniq = remove_l2_elements_from_l1(co_D2_SNV_list, SNV_D2_concurrence)

# get total
total_210 = len(SNV_210_concurrence) + len(SNV_210_mono_uniq) + len(SNV_210_co_uniq)
total_D2 = len(SNV_D2_concurrence) + len(SNV_D2_mono_uniq) + len(SNV_D2_co_uniq)

# write out
write_out(SNV_210_mono_uniq,    SNV_210_mono_uniq_file)
write_out(SNV_210_co_uniq,      SNV_210_co_uniq_file)
write_out(SNV_210_concurrence,  SNV_210_concurrence_file)
write_out(SNV_D2_mono_uniq,     SNV_D2_mono_uniq_file)
write_out(SNV_D2_co_uniq,       SNV_D2_co_uniq_file)
write_out(SNV_D2_concurrence,   SNV_D2_concurrence_file)


###################################################### get matrix ######################################################

print(datetime.now().strftime('[%Y-%m-%d %H:%M:%S]') + ' Get matrix')
# get_matrix_freq(qualified_SNVs_even_flk_depth_file, SNV_210_matrix_file, SNV_D2_matrix_file)
get_matrix(qualified_SNVs_even_flk_depth_file, SNV_210_matrix_file, SNV_210_matrix_file_existence, SNV_D2_matrix_file, SNV_D2_matrix_file_existence)


#################################################### check parallel ####################################################

print(datetime.now().strftime('[%Y-%m-%d %H:%M:%S]') + ' Check parallel')

# get SNV_parallel_dict
SNV_parallel_dict_210 = check_parallel(SNV_210_matrix_file)
SNV_parallel_dict_D2  = check_parallel(SNV_D2_matrix_file)

# merge two dict
SNV_parallel_dict_combined = merge_two_dict(SNV_parallel_dict_210, SNV_parallel_dict_D2)


################################################## get_mutated_genes ###################################################

SNV_matrix_cdc_list = [SNV_210_mono_uniq_file, SNV_210_co_uniq_file, SNV_210_concurrence_file, SNV_D2_mono_uniq_file, SNV_D2_co_uniq_file, SNV_D2_concurrence_file]

SNV_to_gene_dict = {}
SNV_to_ko_id_dict = {}
SNV_to_ko_desc_dict = {}
for SNV_matrix_cdc in SNV_matrix_cdc_list:

    # output files
    SNV_matrix_cdc_path, SNV_matrix_cdc_basename, SNV_matrix_cdc_ext = sep_path_basename_ext(SNV_matrix_cdc)

    output_mutated_genes_cate =         '%s/%s_mutated_genes_category.txt'     % (output_dir, SNV_matrix_cdc_basename)
    output_mutated_genes_cate_fun =     '%s/%s_mutated_genes_cate_fun.txt'     % (output_dir, SNV_matrix_cdc_basename)
    output_summary =                    '%s/%s_summary.txt'                    % (output_dir, SNV_matrix_cdc_basename)
    effect_file =                       '%s/%s_mutation_effect.txt'            % (output_dir, SNV_matrix_cdc_basename)
    output_seq_nc =                     '%s/%s_mutated_genes_nc.fasta'         % (output_dir, SNV_matrix_cdc_basename)
    output_seq_aa =                     '%s/%s_mutated_genes_aa.fasta'         % (output_dir, SNV_matrix_cdc_basename)

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

    output_handle = open(output_summary, 'w')
    for each_snv in open(SNV_matrix_cdc):
        if not each_snv.startswith('\t'):
            each_snv_id =  each_snv.strip().split('\t')[0]
            each_snv_seq = each_snv.strip().split(',')[0].split('|')[0]
            each_snv_pos = int(each_snv.strip().split(',')[0].split('|')[1])
            each_snv_pos_wt = each_snv.strip().split(',')[0].split('|')[2]
            each_snv_pos_v = each_snv.strip().split(',')[0].split('|')[3]

            # get SNV location
            location = ''
            if each_snv_pos in coding_region_dict[each_snv_seq]:
                location = 'Coding'
            else:
                location = 'Intergenic'
                output_handle.write('%s\t%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\n' % (each_snv_id,
                                                                              SNV_parallel_dict_combined[each_snv.strip()],
                                                                              location))
                SNV_to_gene_dict[each_snv_id]    = 'Intergenic'
                SNV_to_ko_id_dict[each_snv_id]   = 'NA'
                SNV_to_ko_desc_dict[each_snv_id] = 'NA'

            # get current gene's COG ID, description, mutation type and mutation effect
            if location == 'Coding':
                for each_gene in ORF_ending_pos_dict:
                    if (each_snv_seq == ORF_seq_id_dict[each_gene]) and (ORF_ending_pos_dict[each_gene][0] <= each_snv_pos <= ORF_ending_pos_dict[each_gene][1]):

                        start_pos_raw = ORF_ending_pos_dict[each_gene][0]
                        end_pos_raw = ORF_ending_pos_dict[each_gene][1]
                        snv_pos_raw = each_snv_pos
                        start_pos_rescaled = ORF_ending_pos_dict[each_gene][0] - (ORF_ending_pos_dict[each_gene][0] - 1)
                        end_pos_rescaled = ORF_ending_pos_dict[each_gene][1] - (ORF_ending_pos_dict[each_gene][0] - 1)
                        snv_pos_rescaled = each_snv_pos - (ORF_ending_pos_dict[each_gene][0] - 1)
                        mutation_type = ''

                        # get current gene's COG ID and description
                        current_COG_id = ''
                        current_COG_function = ''
                        if each_gene in gene_KEGG_id_dict:
                            current_COG_id = gene_KEGG_id_dict[each_gene]
                            current_COG_function = gene_KEGG_function_dict[each_gene]
                        else:
                            current_COG_id = 'NA'
                            current_COG_function = 'NA'

                        # get current gene's mutation type
                        if each_snv_pos_v == '-':
                            mutation_type = 'Frameshift'
                            output_handle.write('%s\t%s\t%s\t%s\t%s\tNA\t%s\t%s\t%s\n' % (each_snv_id,
                                                                                          SNV_parallel_dict_combined[each_snv.strip()],
                                                                                          location,
                                                                                          ORF_strand_dict[each_gene],
                                                                                          each_gene,
                                                                                          mutation_type,
                                                                                          gene_KEGG_id_dict[each_gene],
                                                                                          gene_KEGG_function_dict[each_gene]))
                            SNV_to_gene_dict[each_snv_id]    = each_gene
                            SNV_to_ko_id_dict[each_snv_id]   = gene_KEGG_id_dict[each_gene]
                            SNV_to_ko_desc_dict[each_snv_id] = gene_KEGG_function_dict[each_gene]

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
                            for_write = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (each_snv_id,
                                                                                  SNV_parallel_dict_combined[each_snv.strip()],
                                                                                  location,
                                                                                  ORF_strand_dict[each_gene],
                                                                                  each_gene,
                                                                                  aa_mutation,
                                                                                  mutation_type_term,
                                                                                  gene_KEGG_id_dict[each_gene],
                                                                                  gene_KEGG_function_dict[each_gene])
                            SNV_to_gene_dict[each_snv_id]    = each_gene
                            SNV_to_ko_id_dict[each_snv_id]   = gene_KEGG_id_dict[each_gene]
                            SNV_to_ko_desc_dict[each_snv_id] = gene_KEGG_function_dict[each_gene]

                            output_handle.write(for_write)
    output_handle.close()

    # get mutation_cate_summary
    mutation_cate_dict = get_mutation_cate_summary(output_summary, output_mutated_genes_cate)
    os.system('rm %s' % output_mutated_genes_cate)


    ################################################### export sequences ###################################################

    # get the sequence of affected genes
    output_seq_nc_handle = open(output_seq_nc, 'w')
    for each_gene_nc in SeqIO.parse(combined_ref_ffn, 'fasta'):
        if each_gene_nc.id in mutation_cate_dict:
            SeqIO.write(each_gene_nc, output_seq_nc_handle, 'fasta')
    output_seq_nc_handle.close()

    output_seq_aa_handle = open(output_seq_aa, 'w')
    for each_gene_aa in SeqIO.parse(combined_ref_faa, 'fasta'):
        if each_gene_aa.id in mutation_cate_dict:
            SeqIO.write(each_gene_aa, output_seq_aa_handle, 'fasta')
    output_seq_aa_handle.close()


############################################## get gene level data matrix ##############################################

print('Get gene level data matrix')

SNV_210_matrix_file_with_gene =         '%s/SNV_QC_ncd_even_flk_depth_210_matrix_with_gene.txt'     % output_dir
SNV_D2_matrix_file_with_gene =          '%s/SNV_QC_ncd_even_flk_depth_D2_matrix_with_gene.txt'      % output_dir
SNV_210_matrix_file_gene_level =        '%s/SNV_QC_ncd_even_flk_depth_210_matrix_gene_level.txt'   % output_dir
SNV_D2_matrix_file_gene_level  =        '%s/SNV_QC_ncd_even_flk_depth_D2_matrix_gene_level.txt'    % output_dir

SNV_210_gene_level_dict = {}
SNV_210_matrix_file_with_gene_handle = open(SNV_210_matrix_file_with_gene, 'w')
for SNV_210 in open(SNV_210_matrix_file):
    if SNV_210.startswith('	1D9	1D18	1D27	1D42'):
        SNV_210_matrix_file_with_gene_handle.write('Gene\tSNV\t%s' % SNV_210)
    else:
        SNV_210_split = SNV_210.strip().split('\t')
        SNV_210_id = SNV_210_split[0]
        SNV_210_freq_list = SNV_210_split[1:]
        SNV_210_existence_list = []
        for i in SNV_210_freq_list:
            if i == '0':
                SNV_210_existence_list.append('0')
            else:
                SNV_210_existence_list.append('1')

        SNV_210_affected_gene = SNV_to_gene_dict[SNV_210_id]
        if SNV_210_affected_gene != 'Intergenic':
            SNV_210_matrix_file_with_gene_handle.write('%s\t%s' % (SNV_210_affected_gene, SNV_210))

            if SNV_210_affected_gene not in SNV_210_gene_level_dict:
                SNV_210_gene_level_dict[SNV_210_affected_gene] = SNV_210_existence_list
            else:
                SNV_210_existence_list_updated = []
                for (x, y) in zip(SNV_210_gene_level_dict[SNV_210_affected_gene], SNV_210_existence_list):
                    if (x == '1') or (y == '1'):
                        SNV_210_existence_list_updated.append('1')
                    else:
                        SNV_210_existence_list_updated.append('0')

                SNV_210_gene_level_dict[SNV_210_affected_gene] = SNV_210_existence_list_updated

SNV_210_matrix_file_with_gene_handle.close()


SNV_D2_gene_level_dict = {}
SNV_D2_matrix_file_with_gene_handle = open(SNV_D2_matrix_file_with_gene, 'w')
for SNV_D2 in open(SNV_D2_matrix_file):
    if SNV_D2.startswith('	2D9	2D18	2D27	2D42'):
        SNV_D2_matrix_file_with_gene_handle.write('Gene\tSNV\t%s' % SNV_D2)
    else:
        SNV_D2_split = SNV_D2.strip().split('\t')
        SNV_D2_id = SNV_D2_split[0]
        SNV_D2_freq_list = SNV_D2_split[1:]
        SNV_D2_existence_list = []
        for i in SNV_D2_freq_list:
            if i == '0':
                SNV_D2_existence_list.append('0')
            else:
                SNV_D2_existence_list.append('1')

        SNV_D2_affected_gene = SNV_to_gene_dict[SNV_D2_id]
        if SNV_D2_affected_gene != 'Intergenic':
            SNV_D2_matrix_file_with_gene_handle.write('%s\t%s' % (SNV_D2_affected_gene, SNV_D2))

            if SNV_D2_affected_gene not in SNV_D2_gene_level_dict:
                SNV_D2_gene_level_dict[SNV_D2_affected_gene] = SNV_D2_existence_list
            else:
                SNV_D2_existence_list_updated = []
                for (x, y) in zip(SNV_D2_gene_level_dict[SNV_D2_affected_gene], SNV_D2_existence_list):
                    if (x == '1') or (y == '1'):
                        SNV_D2_existence_list_updated.append('1')
                    else:
                        SNV_D2_existence_list_updated.append('0')

                SNV_D2_gene_level_dict[SNV_D2_affected_gene] = SNV_D2_existence_list_updated

SNV_D2_matrix_file_with_gene_handle.close()


SNV_210_matrix_file_gene_level_handle = open(SNV_210_matrix_file_gene_level, 'w')
SNV_210_matrix_file_gene_level_handle.write('	1D9	1D18	1D27	1D42	5D9	5D18	5D27	5D42	9D9	9D18	9D27	9D42	4D9	4D18	4D27	4D42	8D9	8D18	8D27	8D42	12D9	12D18	12D27	12D42\n')
for affected_gene_210 in SNV_210_gene_level_dict:
    SNV_210_matrix_file_gene_level_handle.write('%s\t%s\n' % (affected_gene_210, '\t'.join(SNV_210_gene_level_dict[affected_gene_210])))
SNV_210_matrix_file_gene_level_handle.close()


SNV_D2_matrix_file_gene_level_handle = open(SNV_D2_matrix_file_gene_level, 'w')
SNV_D2_matrix_file_gene_level_handle.write('	2D9	2D18	2D27	2D42	6D9	6D18	6D27	6D42	10D9	10D18	10D27	10D42	4D9	4D18	4D27	4D42	8D9	8D18	8D27	8D42	12D9	12D18	12D27	12D42\n')
for affected_gene_D2 in SNV_D2_gene_level_dict:
    SNV_D2_matrix_file_gene_level_handle.write('%s\t%s\n' % (affected_gene_D2, '\t'.join(SNV_D2_gene_level_dict[affected_gene_D2])))
SNV_D2_matrix_file_gene_level_handle.close()


############################################## get summary at gene level ###############################################

# get_file name list
gene_level_output_handle = open(gene_level_output, 'w')
for each_strain in ['210', 'D2']:
    affected_gene_concurrence_dict = {}
    affected_gene_uniq_list = []
    for each_parallel_cate in ['monoculture_uniq', 'coculture_uniq', 'concurrence']:

        # parse summary file
        summary_filename = '%s/SNV_QC_ncd_even_flk_depth_%s_%s_summary.txt' % (output_dir, each_strain, each_parallel_cate)
        for each_SNV in open(summary_filename):
            each_SNV_split = each_SNV.strip().split('\t')
            SNV_id = each_SNV_split[0]
            SNV_parallel_cate = each_SNV_split[1][-8:-1]
            affected_gene_id = each_SNV_split[4]

            if affected_gene_id != 'NA':

                if affected_gene_id not in affected_gene_concurrence_dict:
                    affected_gene_concurrence_dict[affected_gene_id] = [SNV_parallel_cate]
                else:
                    affected_gene_concurrence_dict[affected_gene_id].append(SNV_parallel_cate)

                # add to current strain affected_gene_uniq_list
                if affected_gene_id not in affected_gene_uniq_list:
                    affected_gene_uniq_list.append(affected_gene_id)

    for each_affected_gene in affected_gene_concurrence_dict:
        six_digit_1_sum = 0
        six_digit_2_sum = 0
        six_digit_3_sum = 0
        six_digit_4_sum = 0
        six_digit_5_sum = 0
        six_digit_6_sum = 0
        for each_six_digit in affected_gene_concurrence_dict[each_affected_gene]:
            six_digit_1_sum += int(each_six_digit[0])
            six_digit_2_sum += int(each_six_digit[1])
            six_digit_3_sum += int(each_six_digit[2])
            six_digit_4_sum += int(each_six_digit[4])
            six_digit_5_sum += int(each_six_digit[5])
            six_digit_6_sum += int(each_six_digit[6])

        six_digit_sum = '%s%s%s_%s%s%s(%s)' % (six_digit_1_sum, six_digit_2_sum, six_digit_3_sum, six_digit_4_sum, six_digit_5_sum, six_digit_6_sum, len(affected_gene_concurrence_dict[each_affected_gene]))

        # get report for each gene
        for_report = ''
        if six_digit_sum[0:3] == '000':
            for_report = '%s\t%s\tco\t%s\t%s\n' % (each_affected_gene, six_digit_sum, gene_KEGG_id_dict[each_affected_gene], gene_KEGG_function_dict[each_affected_gene])
        elif six_digit_sum[4:7] == '000':
            for_report = '%s\t%s\tmono\t%s\t%s\n' % (each_affected_gene, six_digit_sum, gene_KEGG_id_dict[each_affected_gene], gene_KEGG_function_dict[each_affected_gene])
        elif (six_digit_sum[0:3] != '000') and (six_digit_sum[4:7] != '000'):
            for_report = '%s\t%s\tboth\t%s\t%s\n' % (each_affected_gene, six_digit_sum, gene_KEGG_id_dict[each_affected_gene], gene_KEGG_function_dict[each_affected_gene])

        # write out
        gene_level_output_handle.write(for_report)

    # print the report
    print(datetime.now().strftime('[%Y-%m-%d %H:%M:%S]') + ' The total number of %s affected gene (ncd_even_depth): %s' % (each_strain, len(affected_gene_uniq_list)))

gene_level_output_handle.close()

# for report
print(datetime.now().strftime('[%Y-%m-%d %H:%M:%S]') + ' Summaries at gene level exported to %s' % gene_level_output)


################################################# plot SNV frequency ###################################################

print(datetime.now().strftime('[%Y-%m-%d %H:%M:%S]') + ' plot SNV frequency')

if plot_snv_freq is True:
    force_create_folder(mutation_effect_plot_folder)
    plot_freq(SNV_210_matrix_file, mutation_effect_plot_folder, SNV_to_gene_dict, SNV_to_ko_id_dict, SNV_to_ko_desc_dict)
    plot_freq(SNV_D2_matrix_file, mutation_effect_plot_folder, SNV_to_gene_dict, SNV_to_ko_id_dict, SNV_to_ko_desc_dict)


print(datetime.now().strftime('[%Y-%m-%d %H:%M:%S]') + ' Done!')

print('\n################################################################################################')
