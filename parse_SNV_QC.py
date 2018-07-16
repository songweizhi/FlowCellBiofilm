import os


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


def get_summary_matrix(SNV_file_in):
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

    SNV_file_out_existence_210 = '%s_summary_210_existence.txt' % SNV_file_in_basename
    SNV_file_out_frequency_210 = '%s_summary_210_frequency.txt' % SNV_file_in_basename
    SNV_file_out_existence_D2 = '%s_summary_D2_existence.txt' % SNV_file_in_basename
    SNV_file_out_frequency_D2 = '%s_summary_D2_frequency.txt' % SNV_file_in_basename

    SNV_file_out_existence_210_cdc = '%s_summary_210_existence_cdc.txt' % SNV_file_in_basename
    SNV_file_out_frequency_210_cdc = '%s_summary_210_frequency_cdc.txt' % SNV_file_in_basename
    SNV_file_out_existence_D2_cdc = '%s_summary_D2_existence_cdc.txt' % SNV_file_in_basename
    SNV_file_out_frequency_D2_cdc = '%s_summary_D2_frequency_cdc.txt' % SNV_file_in_basename

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


########################################## specify input files and parameters ##########################################

wd = '/Users/songweizhi/Desktop/666666'
SNV_quality_file = 'deepSNV_output_combined_QC.txt'
min_var_reads_num = 10
min_at_each_direction = 1
strand_bias_cutoff = 20
depth_difference_cutoff = 15
mean_depth_len = 2000
os.chdir(wd)


############################################ define the name of output files ###########################################

matrix_header_210 = ['1D9', '1D18', '1D27', '1D42', '5D9', '5D18', '5D27', '5D42', '9D9', '9D18', '9D27', '9D42', '4D9', '4D18', '4D27', '4D42', '8D9', '8D18', '8D27', '8D42', '12D9', '12D18', '12D27', '12D42']
matrix_header_D2 = ['2D9', '2D18', '2D27', '2D42', '6D9', '6D18', '6D27', '6D42', '10D9', '10D18', '10D27', '10D42', '4D9', '4D18', '4D27', '4D42', '8D9', '8D18', '8D27', '8D42', '12D9', '12D18', '12D27', '12D42']

qualified_SNVs_even_flanking_depth_file = 'deepSNV_output_qualified_even_flanking_depth.txt'
qualified_SNVs_diff_flanking_depth_file = 'deepSNV_output_qualified_diff_flanking_depth.txt'


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


############################## get matrix for SNVs with similar/different flanking depth ###############################

get_summary_matrix(qualified_SNVs_even_flanking_depth_file)
get_summary_matrix(qualified_SNVs_diff_flanking_depth_file)








