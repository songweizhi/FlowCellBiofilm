import os
import glob
import argparse


############################## specify input files ##############################

# parser = argparse.ArgumentParser()

# parser.add_argument('-in', required=True, help='deepSNV output folder')
# args = vars(parser.parse_args())

# deepSNV_output_folder = args['in']

deepSNV_output_folder = '/Users/songweizhi/Desktop/4_1_deepSNV'

#################################################################################

sample_dict = {'1': 'Mono210_A',
               '5': 'Mono210_B',
               '9': 'Mono210_C',
               '2': 'MonoD2_A',
               '6': 'MonoD2_B',
               '10': 'MonoD2_C',
               '4': 'Coculture_A',
               '8': 'Coculture_B',
               '12': 'Coculture_C'}

replicate_dict = {'1': 'A',
                  '5': 'B',
                  '9': 'C',
                  '2': 'A',
                  '6': 'B',
                  '10': 'C',
                  '4': 'A',
                  '8': 'B',
                  '12': 'C'}

sample_list = ['Mono210_A',
               'Mono210_B',
               'Mono210_C',
               'MonoD2_A',
               'MonoD2_B',
               'MonoD2_C',
               'Coculture_A',
               'Coculture_B',
               'Coculture_C']


experiment_setup_dict = {'Mono210': ['1', '5', '9'],
                         'MonoD2': ['2', '6', '10'],
                         'Coculture': ['4', '8', '12']}

sample_id_list_all = [1, 5, 9, 2, 6, 10, 4, 8, 12]
sample_id_list_210 = [1, 5, 9, 4, 8, 12]
sample_id_list_D2 =  [2, 6, 10, 4, 8, 12]

timepoint_list = ['D9', 'D18', 'D27', 'D42']
summary_file_header_all = ['1D9', '1D18', '1D27', '1D42', '5D9', '5D18', '5D27', '5D42', '9D9', '9D18', '9D27', '9D42', '2D9', '2D18', '2D27', '2D42', '6D9', '6D18', '6D27', '6D42', '10D9', '10D18', '10D27', '10D42', '4D9', '4D18', '4D27', '4D42', '8D9', '8D18', '8D27', '8D42', '12D9', '12D18', '12D27', '12D42']
summary_file_header_210 = ['1D9', '1D18', '1D27', '1D42', '5D9', '5D18', '5D27', '5D42', '9D9', '9D18', '9D27', '9D42', '4D9', '4D18', '4D27', '4D42', '8D9', '8D18', '8D27', '8D42', '12D9', '12D18', '12D27', '12D42']
summary_file_header_D2 =  ['2D9', '2D18', '2D27', '2D42', '6D9', '6D18', '6D27', '6D42', '10D9', '10D18', '10D27', '10D42', '4D9', '4D18', '4D27', '4D42', '8D9', '8D18', '8D27', '8D42', '12D9', '12D18', '12D27', '12D42']

# combined_existence_all = 'deepSNV_output_summary_all_existence.txt'
# combined_existence_210 = 'deepSNV_output_summary_210_existence.txt'
# combined_existence_D2 =  'deepSNV_output_summary_D2_existence.txt'
#
# combined_frequency_all = 'deepSNV_output_summary_all_frequency.txt'
# combined_frequency_210 = 'deepSNV_output_summary_210_frequency.txt'
# combined_frequency_D2 =  'deepSNV_output_summary_D2_frequency.txt'
#
# factor_file_all = 'deepSNV_output_summary_all_factor.txt'
# factor_file_210 = 'deepSNV_output_summary_210_factor.txt'
# factor_file_D2 =  'deepSNV_output_summary_D2_factor.txt'
#
#
# get deepSNV output file list
deepSNV_output_file_re = '%s/*.txt' % deepSNV_output_folder
deepSNV_output_file_list = [os.path.basename(file_name) for file_name in glob.glob(deepSNV_output_file_re)]

# # get total snv list and their proportions
# snv_proportion_dict = {}
# total_snv_list = []
# for each_deepSNV_output in deepSNV_output_file_list:
#     pwd_each_deepSNV_output = '%s/%s' % (deepSNV_output_folder, each_deepSNV_output)
#     treatment_id = each_deepSNV_output.split('_vs_')[0]
#     for each_snv in open(pwd_each_deepSNV_output):
#         # ignore the first line
#         if not each_snv.startswith('chr,pos,ref,var'):
#             each_snv_split = each_snv.strip().split(',')
#             snv = '|'.join(each_snv_split[0:4])
#             proportion = str(float("{0:.3f}".format(float(each_snv_split[5]) * 100)))
#             dict_key = '%s___%s' % (snv, treatment_id)
#             snv_proportion_dict[dict_key] = proportion
#             if snv not in total_snv_list:
#                 total_snv_list.append(snv)
#
# # sort total snv list
# total_snv_list = sorted(total_snv_list)
#
#
# # get summary table
# combined_existence_all_handle = open(combined_existence_all, 'w')
# combined_existence_210_handle = open(combined_existence_210, 'w')
# combined_existence_D2_handle = open(combined_existence_D2, 'w')
# combined_frequency_all_handle = open(combined_frequency_all, 'w')
# combined_frequency_210_handle = open(combined_frequency_210, 'w')
# combined_frequency_D2_handle = open(combined_frequency_D2, 'w')
#
# combined_existence_all_handle.write('\t%s\n' % '\t'.join(summary_file_header_all))
# combined_existence_210_handle.write('\t%s\n' % '\t'.join(summary_file_header_210))
# combined_existence_D2_handle.write('\t%s\n' % '\t'.join(summary_file_header_D2))
# combined_frequency_all_handle.write('\t%s\n' % '\t'.join(summary_file_header_all))
# combined_frequency_210_handle.write('\t%s\n' % '\t'.join(summary_file_header_210))
# combined_frequency_D2_handle.write('\t%s\n' % '\t'.join(summary_file_header_D2))
#
# for each_snv in total_snv_list:
#     occurrence_list_p_all = []
#     occurrence_list_p_210 = []
#     occurrence_list_p_D2 = []
#     occurrence_list_e_all = []
#     occurrence_list_e_210 = []
#     occurrence_list_e_D2 = []
#
#     # for all samples
#     for each_sample_id in sample_id_list_all:
#         for each_timepoint in timepoint_list:
#             key = '%s___%s%s' % (each_snv, each_sample_id, each_timepoint)
#             if key not in snv_proportion_dict:
#                 occurrence_list_p_all.append('0')
#                 occurrence_list_e_all.append('0')
#             if key in snv_proportion_dict:
#                 occurrence_list_p_all.append(str(snv_proportion_dict[key]))
#                 occurrence_list_e_all.append('1')
#
#     # for 210
#     for each_sample_id in sample_id_list_210:
#         for each_timepoint in timepoint_list:
#             key = '%s___%s%s' % (each_snv, each_sample_id, each_timepoint)
#             if key.startswith('2.10'):
#                 if key not in snv_proportion_dict:
#                     occurrence_list_p_210.append('0')
#                     occurrence_list_e_210.append('0')
#                 if key in snv_proportion_dict:
#                     occurrence_list_p_210.append(str(snv_proportion_dict[key]))
#                     occurrence_list_e_210.append('1')
#
#     # for D2
#     for each_sample_id in sample_id_list_D2:
#         for each_timepoint in timepoint_list:
#             key = '%s___%s%s' % (each_snv, each_sample_id, each_timepoint)
#             if key.startswith('D2'):
#                 if key not in snv_proportion_dict:
#                     occurrence_list_p_D2.append('0')
#                     occurrence_list_e_D2.append('0')
#                 if key in snv_proportion_dict:
#                     occurrence_list_p_D2.append(str(snv_proportion_dict[key]))
#                     occurrence_list_e_D2.append('1')
#
#     # write out for all
#     combined_existence_all_handle.write('%s\t%s\n' % (each_snv, '\t'.join(occurrence_list_e_all)))
#     combined_frequency_all_handle.write('%s\t%s\n' % (each_snv, '\t'.join(occurrence_list_p_all)))
#
#     # write out for 210
#     if each_snv.startswith('2.10'):
#         combined_existence_210_handle.write('%s\t%s\n' % (each_snv, '\t'.join(occurrence_list_e_210)))
#         combined_frequency_210_handle.write('%s\t%s\n' % (each_snv, '\t'.join(occurrence_list_p_210)))
#
#     # write out for D2
#     if each_snv.startswith('D2'):
#         combined_existence_D2_handle.write('%s\t%s\n' % (each_snv, '\t'.join(occurrence_list_e_D2)))
#         combined_frequency_D2_handle.write('%s\t%s\n' % (each_snv, '\t'.join(occurrence_list_p_D2)))
#
# combined_existence_all_handle.close()
# combined_existence_210_handle.close()
# combined_existence_D2_handle.close()
# combined_frequency_all_handle.close()
# combined_frequency_210_handle.close()
# combined_frequency_D2_handle.close()
#
#
# # prepare factor file for all
# factor_file_all_handle = open(factor_file_all, 'w')
# factor_file_all_handle.write('Sample\tSpecies\tReplicate\tTime\tLabel\n')
# for each_sample in summary_file_header_all:
#     sample_id = each_sample.split('D')[0]
#     timepoint = 'D%s' % each_sample.split('D')[1]
#     species = sample_dict[sample_id][:-2]
#     treatment = replicate_dict[sample_id]
#     label = sample_dict[sample_id]
#     for_write_out = '%s\t%s\t%s\t%s\t%s\n' % (each_sample, species, treatment, timepoint, label)
#     factor_file_all_handle.write(for_write_out)
# factor_file_all_handle.close()
#
# # prepare factor file for 210
# factor_file_210_handle = open(factor_file_210, 'w')
# factor_file_210_handle.write('Sample\tSpecies\tReplicate\tTime\tLabel\n')
# for each_sample in summary_file_header_210:
#     sample_id = each_sample.split('D')[0]
#     timepoint = 'D%s' % each_sample.split('D')[1]
#     species = sample_dict[sample_id][:-2]
#     treatment = replicate_dict[sample_id]
#     label = sample_dict[sample_id]
#
#     for_write_out = '%s\t%s\t%s\t%s\t%s\n' % (each_sample, species, treatment, timepoint, label)
#     factor_file_210_handle.write(for_write_out)
# factor_file_210_handle.close()
#
# # prepare factor file for D2
# factor_file_D2_handle = open(factor_file_D2, 'w')
# factor_file_D2_handle.write('Sample\tSpecies\tReplicate\tTime\tLabel\n')
# for each_sample in summary_file_header_D2:
#     sample_id = each_sample.split('D')[0]
#     timepoint = 'D%s' % each_sample.split('D')[1]
#     species = sample_dict[sample_id][:-2]
#     treatment = replicate_dict[sample_id]
#     label = sample_dict[sample_id]
#
#     for_write_out = '%s\t%s\t%s\t%s\t%s\n' % (each_sample, species, treatment, timepoint, label)
#     factor_file_D2_handle.write(for_write_out)
# factor_file_D2_handle.close()
#
