import os
import shutil
import numpy as np
from scipy import stats
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


# specify wd and input files
wd = '/Users/songweizhi/Desktop/000'
deepSNV_op_cdc = 'SNV_QC_even_depth_matrix_210_frequency_cdc.txt'
deepSNV_op_cdc = 'SNV_QC_even_depth_matrix_D2_frequency_cdc.txt'
output_mutation_effect = 'SNV_mutation_effect.txt'
mutation_affect_plot_folder = 'SNV_mutation_effect_plot'
plot_mutation_effect = 1


# forward to working directory
os.chdir(wd)

# get_mutation_effect
get_mutation_effect(deepSNV_op_cdc, plot_mutation_effect)




