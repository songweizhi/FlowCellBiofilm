import os
import glob


os.chdir('/Users/songweizhi/Dropbox/Research/Flow_cell/deepSNV_output_summary')
deepSNV_output_folder = '/Users/songweizhi/Dropbox/Research/Flow_cell/deepSNV_output_summary/4_deepSNV'


combined_file_sorted =                       'combined_deepSNV_output_sorted.txt'
combined_file_sorted_tmp1 =                  'combined_deepSNV_output_sorted_tmp1.txt'
combined_summary_existence =                 'combined_summary_existence.txt'
combined_summary_frequency =                 'combined_summary_frequency.txt'
combined_summary_existence_sorted =          'combined_summary_existence_sorted.txt'
combined_summary_frequency_sorted =          'combined_summary_frequency_sorted.txt'
combined_summary_existence_sorted_with_num = 'deepSNV_output_summary_existence.txt'
combined_summary_frequency_sorted_with_num = 'deepSNV_output_summary_frequency.txt'


experiment_setup_dict = {'Coculture': ['4', '8', '12'], 'Mono_210': ['1', '5', '9'], 'Mono_D2': ['2', '6', '10']}
timepoint_list = ['D9', 'D18', 'D27', 'D42']


combined_summary_existence_handle = open(combined_summary_existence, 'w')
combined_summary_frequency_handle = open(combined_summary_frequency, 'w')
for each_experiment_setup in experiment_setup_dict:

    treatment_id_list = experiment_setup_dict[each_experiment_setup]
    snv_group_existence_dict = {}
    snv_group_frequency_dict = {}
    snv_group_uniq = []
    for treatment_id in treatment_id_list:

        # get deepSNV output file list
        deepSNV_output_files = '%s/%sD*.txt' % (deepSNV_output_folder, treatment_id)
        file_list = [os.path.basename(file_name) for file_name in glob.glob(deepSNV_output_files)]

        if file_list != []:
            # combine all results together
            combined_file = 'combined_deepSNV_output_treatment_%s.txt' % treatment_id
            combined_file_sorted = 'combined_deepSNV_output_treatment_%s_sorted.txt' % treatment_id

            combined_file_handle = open(combined_file, 'w')
            for each_file in file_list:
                pwd_each_file = '%s/%s' % (deepSNV_output_folder, each_file)
                each_file_basename, each_file_extension = os.path.splitext(each_file)
                day = 'D%s' % each_file_basename.split('_')[0].split('D')[1]
                for each_snv in open(pwd_each_file):
                    if not each_snv.startswith('chr,pos,ref'):
                        combined_file_handle.write('%s,%s\n' % (day, each_snv.strip()))
            combined_file_handle.close()

            # sort combined file
            os.system('cat %s | sort > %s' % (combined_file, combined_file_sorted))

            # get frequency dict
            frequency_dict = {}
            for each5 in open(combined_file_sorted):
                each5_split = each5.strip().split(',')
                key_f = '%s|%s|%s|%s>%s,%s' % (treatment_id, each5_split[1], each5_split[2], each5_split[3], each5_split[4], each5_split[0])
                key_f_value = each5_split[6]
                frequency_dict[key_f] = str(float("{0:.3f}".format(float(key_f_value) * 100)))

            # keep needed region
            snv_overall = []
            snv_uniq = []
            for each3 in open(combined_file_sorted):
                each3_split = each3.strip().split(',')
                day = each3_split[0]
                snv_no_treatment_id = '%s|%s|%s>%s' % (each3_split[1], each3_split[2], each3_split[3], each3_split[4])

                if snv_no_treatment_id not in snv_group_uniq:
                    snv_group_uniq.append(snv_no_treatment_id)

                snv = '%s|%s|%s|%s>%s' % (treatment_id, each3_split[1], each3_split[2], each3_split[3], each3_split[4])
                needed = '%s,%s' % (snv, day)
                snv_overall.append(needed)
                if snv not in snv_uniq:
                    snv_uniq.append(snv)

            # get table
            for each_uniq_snv in snv_uniq:
                each_uniq_snv_existence_list = []
                each_uniq_snv_frequency_list = []
                for time_point in timepoint_list:
                    key = '%s,%s' % (each_uniq_snv, time_point)
                    if key in snv_overall:
                        existence = '1'
                        frequency = frequency_dict[key]
                    else:
                        existence = '0'
                        frequency = '0'
                    each_uniq_snv_existence_list.append(existence)
                    each_uniq_snv_frequency_list.append(frequency)
                snv_group_existence_dict[each_uniq_snv] = ' '.join(each_uniq_snv_existence_list)
                snv_group_frequency_dict[each_uniq_snv] = ' '.join(each_uniq_snv_frequency_list)

    # combine three replicates together
    for each_snv_group_uniq in snv_group_uniq:
        parallel_existence_list = []
        parallel_frequency_list = []

        for treatment_id2 in treatment_id_list:
            key_group = '%s|%s' % (treatment_id2, each_snv_group_uniq)
            if key_group in snv_group_existence_dict:
                parallel_existence_list.append(snv_group_existence_dict[key_group])
                parallel_frequency_list.append(snv_group_frequency_dict[key_group])
            else:
                parallel_existence_list.append('0 0 0 0')
                parallel_frequency_list.append('0 0 0 0')
        combined_summary_existence_handle.write('%s\t\t%s\t\t%s\n' % (each_experiment_setup, each_snv_group_uniq, '\t\t'.join(parallel_existence_list)))
        combined_summary_frequency_handle.write('%s\t\t%s\t\t%s\n' % (each_experiment_setup, each_snv_group_uniq, '\t\t'.join(parallel_frequency_list)))
combined_summary_existence_handle.close()
combined_summary_frequency_handle.close()


# sort summary
os.system('cat %s | sort > %s' % (combined_summary_existence, combined_summary_existence_sorted))
os.system('cat %s | sort > %s' % (combined_summary_frequency, combined_summary_frequency_sorted))


# add detected number
combined_summary_sorted_with_num_handle_existence = open(combined_summary_existence_sorted_with_num, 'w')
#combined_summary_sorted_with_num_handle_existence.write('Treatment\t\tSNV\t\tD9A D18A D27A D42A\t\tD9B D18B D27B D42B\t\tD9C D18C D27C D42C\t\tNumber\n')
for each4 in open(combined_summary_existence_sorted):
    each4_split = each4.strip().split('\t\t')
    detected_number = 3 - each4_split.count('0 0 0 0')
    combined_summary_sorted_with_num_handle_existence.write('%s\t\t%s\n' % (each4.strip(), detected_number))
combined_summary_sorted_with_num_handle_existence.close()

combined_summary_sorted_with_num_handle_frequency = open(combined_summary_frequency_sorted_with_num, 'w')
#combined_summary_sorted_with_num_handle_frequency.write('Treatment\t\tSNV\t\tD9A D18A D27A D42A\t\tD9B D18B D27B D42B\t\tD9C D18C D27C D42C\t\tNumber\n')
for each6 in open(combined_summary_frequency_sorted):
    each6_split = each6.strip().split('\t\t')
    detected_number = 3 - each6_split.count('0 0 0 0')
    combined_summary_sorted_with_num_handle_frequency.write('%s\t\t%s\n' % (each6.strip(), detected_number))
combined_summary_sorted_with_num_handle_frequency.close()

# remove tmp files
os.system('rm combined*.txt')
