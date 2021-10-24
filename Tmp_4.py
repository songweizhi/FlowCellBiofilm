import os
import glob


# file in/out
hmmer_output_folder = '/Users/songweizhi/Desktop/HMMER_wd'
file_extension      = 'tbl'
file_out            = '/Users/songweizhi/Desktop/hmmer_output_df.txt'


# get genome list in provided folder
hmmer_output_file_re = '%s/*.%s' % (hmmer_output_folder, file_extension)
hmmer_output_file_list = glob.glob(hmmer_output_file_re)

hmm_id_to_description_dict = {}
genome_list = []
identified_hmm_list = []
hmm_count_dict_of_dict = {}
for hmmer_output_file in hmmer_output_file_list:
    hmmer_output_file_name = hmmer_output_file.split('/')[-1]
    current_genome_hmm_count_dict = {}
    for each_line in open(hmmer_output_file):
        if not each_line.startswith('#'):
            each_line_split = each_line.strip().split(' ')
            each_line_split_no_empty_element = [i for i in each_line_split if i != '']
            gene_id = each_line_split_no_empty_element[0]
            hmm_id = each_line_split_no_empty_element[4]
            description = '_'.join(each_line_split_no_empty_element[22:])

            if hmm_id not in hmm_id_to_description_dict:
                hmm_id_to_description_dict[hmm_id] = [description]
            else:
                if description not in hmm_id_to_description_dict[hmm_id]:
                    hmm_id_to_description_dict[hmm_id].append(description)

            if hmm_id not in identified_hmm_list:
                identified_hmm_list.append(hmm_id)

            if hmm_id not in current_genome_hmm_count_dict:
                current_genome_hmm_count_dict[hmm_id] = 1
            else:
                current_genome_hmm_count_dict[hmm_id] += 1
    hmm_count_dict_of_dict[hmmer_output_file_name] = current_genome_hmm_count_dict
    genome_list.append(hmmer_output_file_name)


identified_hmm_list_sorted = sorted(identified_hmm_list)


file_out_handle = open(file_out, 'w')
file_out_handle.write('Genome\t%s\n' % '\t'.join(identified_hmm_list_sorted))
for genome in sorted(genome_list):
    current_genome_hmm_stats_dict = hmm_count_dict_of_dict[genome]
    current_genome_hmm_stats_list = []
    for hmm in identified_hmm_list_sorted:
        if hmm in current_genome_hmm_stats_dict:
            current_genome_hmm_stats_list.append(current_genome_hmm_stats_dict[hmm])
        else:
            current_genome_hmm_stats_list.append(0)

    current_genome_hmm_stats_list_str = [str(i) for i in current_genome_hmm_stats_list]
    for_write_out = '%s\t%s\n' % (genome, '\t'.join(current_genome_hmm_stats_list_str))
    file_out_handle.write(for_write_out)
file_out_handle.close()


print(hmm_id_to_description_dict)