import os
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt


def get_affect(frequency_list, plot_filename):
    time_point = [9, 18, 27, 42]
    time_point_rescaled = []
    for each_tp in time_point:
        each_tp_rescaled = each_tp / 42
        time_point_rescaled.append(each_tp_rescaled)
    time_point_rescaled_arrary = np.array(time_point_rescaled)

    frequency_arrary = np.array(frequency_list)
    min_frequency = frequency_arrary.min()
    max_frequency = frequency_arrary.max()
    slope, intercept, r_value, p_value, std_err = stats.linregress(time_point_rescaled_arrary, frequency_arrary)

    if frequency_arrary[-1] >= 25:
        affect = 'Beneficial'
    elif (slope > 0) and (max_frequency >= 10):
        affect = 'Beneficial'
    elif (slope < 0) and (max_frequency <= 10):
        affect = 'Harmful'
    else:
        affect = 'Neutral'

    # plot
    plt.plot(time_point_rescaled_arrary, frequency_arrary, 'o')
    plt.plot(time_point_rescaled_arrary, intercept + slope*time_point_rescaled_arrary, 'r')

    # # add text
    # y_min = plt.ylim()[0]  # get the y-axes minimum value
    # y_max = plt.ylim()[1]  # get the y-axes maximum value
    #
    # # set text position
    # text_x = 0.2
    # text_y_slope = y_min + (y_max - y_min) / 5 * 4.4
    # text_y_p_value = y_min + (y_max - y_min) / 5 * 4.1
    # text_y_affect = y_min + (y_max - y_min) / 5 * 3.8
    # plt.text(text_x, text_y_slope, 'Slope: %s' % float("{0:.2f}".format(slope)))
    # plt.text(text_x, text_y_p_value, 'P_value: %s' % float("{0:.2f}".format(p_value)))
    # plt.text(text_x, text_y_affect, 'Affect: %s' % affect)
    #
    # plt.xlabel('Frequency')
    # plt.ylabel('Time point')
    # plt.savefig('%s.png' % plot_filename, dpi=300)
    # plt.close()

    return affect


# specify wd and input files
wd = '/Users/songweizhi/Desktop/000'
deepSNV_op_cdc = 'deepSNV_output_summary_all_frequency_cdc.txt'
deepSNV_mutated_genes = 'deepSNV_mutated_gene.txt'
annotation_file = 'deepSNV_mutated_gene_aa.faa.COG.arCOG.kegg.tsv'
output_combined = 'deepSNV_combined_info.txt'

# forward to working directory
os.chdir(wd)

# get co-occurrence and mutation affect
SNV_co_occurrence_dict = {}
for each_snv in open(deepSNV_op_cdc):

    if not each_snv.startswith('\t'):
        each_snv_split = each_snv.strip().split('\t')
        snv_id = each_snv_split[0]
        occur_profile = []

        if (each_snv_split[1] != '0') or (each_snv_split[2] != '0') or (each_snv_split[3] != '0') or (each_snv_split[4] != '0'):
            frequency_list = [float(each_snv_split[1]), float(each_snv_split[2]), float(each_snv_split[3]), float(each_snv_split[4])]
            png_filename = '%s_%s' % (snv_id, 'Mono210_A')
            affect = get_affect(frequency_list, png_filename)
            occur_profile.append('Mono210_A(%s)' % affect)

        if (each_snv_split[5] != '0') or (each_snv_split[6] != '0') or (each_snv_split[7] != '0') or (each_snv_split[8] != '0'):
            frequency_list = [float(each_snv_split[5]), float(each_snv_split[6]), float(each_snv_split[7]), float(each_snv_split[8])]
            png_filename = '%s_%s' % (snv_id, 'Mono210_B')
            affect = get_affect(frequency_list, png_filename)
            occur_profile.append('Mono210_B(%s)' % affect)

        if (each_snv_split[9] != '0') or (each_snv_split[10] != '0') or (each_snv_split[11] != '0') or (each_snv_split[12] != '0'):
            frequency_list = [float(each_snv_split[9]), float(each_snv_split[10]), float(each_snv_split[11]), float(each_snv_split[12])]
            png_filename = '%s_%s' % (snv_id, 'Mono210_C')
            affect = get_affect(frequency_list, png_filename)
            occur_profile.append('Mono210_C(%s)' % affect)

        if (each_snv_split[13] != '0') or (each_snv_split[14] != '0') or (each_snv_split[15] != '0') or (each_snv_split[16] != '0'):
            frequency_list = [float(each_snv_split[13]), float(each_snv_split[14]), float(each_snv_split[15]), float(each_snv_split[16])]
            png_filename = '%s_%s' % (snv_id, 'MonoD2_A')
            affect = get_affect(frequency_list, png_filename)
            occur_profile.append('MonoD2_A(%s)' % affect)

        if (each_snv_split[17] != '0') or (each_snv_split[18] != '0') or (each_snv_split[19] != '0') or (each_snv_split[20] != '0'):
            frequency_list = [float(each_snv_split[17]), float(each_snv_split[18]), float(each_snv_split[19]), float(each_snv_split[20])]
            png_filename = '%s_%s' % (snv_id, 'MonoD2_B')
            affect = get_affect(frequency_list, png_filename)
            occur_profile.append('MonoD2_B(%s)' % affect)

        if (each_snv_split[21] != '0') or (each_snv_split[22] != '0') or (each_snv_split[23] != '0') or (each_snv_split[24] != '0'):
            frequency_list = [float(each_snv_split[21]), float(each_snv_split[22]), float(each_snv_split[23]), float(each_snv_split[24])]
            png_filename = '%s_%s' % (snv_id, 'MonoD2_C')
            affect = get_affect(frequency_list, png_filename)
            occur_profile.append('MonoD2_C(%s)' % affect)

        if (each_snv_split[25] != '0') or (each_snv_split[26] != '0') or (each_snv_split[27] != '0') or (each_snv_split[28] != '0'):
            frequency_list = [float(each_snv_split[25]), float(each_snv_split[26]), float(each_snv_split[27]), float(each_snv_split[28])]
            png_filename = '%s_%s' % (snv_id, 'Coculture_A')
            affect = get_affect(frequency_list, png_filename)
            occur_profile.append('Coculture_A(%s)' % affect)

        if (each_snv_split[29] != '0') or (each_snv_split[30] != '0') or (each_snv_split[31] != '0') or (each_snv_split[32] != '0'):
            frequency_list = [float(each_snv_split[29]), float(each_snv_split[30]), float(each_snv_split[31]), float(each_snv_split[32])]
            png_filename = '%s_%s' % (snv_id, 'Coculture_B')
            affect = get_affect(frequency_list, png_filename)
            occur_profile.append('Coculture_B(%s)' % affect)

        if (each_snv_split[33] != '0') or (each_snv_split[34] != '0') or (each_snv_split[35] != '0') or (each_snv_split[36] != '0'):
            frequency_list = [float(each_snv_split[33]), float(each_snv_split[34]), float(each_snv_split[35]), float(each_snv_split[36])]
            png_filename = '%s_%s' % (snv_id, 'Coculture_C')
            affect = get_affect(frequency_list, png_filename)
            occur_profile.append('Coculture_C(%s)' % affect)

        #print('%s\t%s' % (snv_id, '|'.join(occur_profile)))
        co_occurrence = '|'.join(occur_profile)
        SNV_co_occurrence_dict[snv_id] = co_occurrence


# store annotation results in dicts
gene_KO_id_dict = {}
gene_KO_function_dict = {}
gene_COG_cat_dict = {}
gene_COG_id_dict = {}
gene_COG_function_dict = {}
for each_snv3 in open(annotation_file):
    each_snv3_split = each_snv3.strip().split('\t')
    gene_id = each_snv3_split[0]
    KO_id = each_snv3_split[6]
    if KO_id == 'No_KO':
        KO_id = 'NA'
    KO_function = each_snv3_split[7]
    COG_cat = each_snv3_split[10]
    COG_id = each_snv3_split[8]
    COG_function = each_snv3_split[9]
    gene_KO_id_dict[gene_id] = KO_id
    gene_KO_function_dict[gene_id] = KO_function
    gene_COG_cat_dict[gene_id] = COG_cat
    gene_COG_id_dict[gene_id] = COG_id
    gene_COG_function_dict[gene_id] = COG_function


# combine tables
output_combined_handle = open(output_combined, 'w')

for each_snv2 in open(deepSNV_mutated_genes):
    #print(each_snv2)
    each_snv2_split = each_snv2.strip().split('\t')
    snv_id2 = each_snv2_split[0]
    gene_id2 = each_snv2_split[3]
    current_KO_id = ''
    current_KO_function = ''
    current_COG_id = ''
    current_COG_cat = ''
    current_COG_function = ''

    if gene_id2 in gene_COG_id_dict:
        current_KO_id = gene_KO_id_dict[gene_id2]
        current_KO_function = gene_KO_function_dict[gene_id2]
        current_COG_id = gene_COG_id_dict[gene_id2]
        current_COG_cat = gene_COG_cat_dict[gene_id2]
        current_COG_function = gene_COG_function_dict[gene_id2]
    else:
        current_KO_id = 'NA'
        current_KO_function = 'NA'
        current_COG_id = 'NA'
        current_COG_cat = 'NA'
        current_COG_function = 'NA'

    for_out = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (each_snv2_split[0], SNV_co_occurrence_dict[snv_id2], each_snv2_split[1], each_snv2_split[2], each_snv2_split[3], each_snv2_split[4], each_snv2_split[5], current_KO_id, current_COG_cat, current_COG_id, current_COG_function)
    output_combined_handle.write(for_out)










