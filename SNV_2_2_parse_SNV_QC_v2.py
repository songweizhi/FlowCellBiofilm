import os


############################################### input file and parameters ##############################################

# wd
wd = '/Users/songweizhi/Dropbox/Research/Flow_cell_datasets/parse_SNV_QC_v2'
os.chdir(wd)

# SNV QC output (nonsubsampled: /srv/scratch/z5039045/Flow_cell_biofilm/4_2_SNV_QC/output_f50000bp_1000mer_dl2000bp/SNV_QC.txt)
SNV_QC_nonsubsampled = 'SNV_QC.txt'

# The minimum number of reads harboring SNV
min_var_reads_num = 10

# The minimum number of reads harboring SNV at each direction
min_at_each_direction = 1

# strand_bias cutoff
strand_bias_cutoff = 20

# depth difference cutoff
depth_difference_cutoff = 15

# flanking length for mean depth calculation
mean_depth_len = 2000


######################################################### Main #########################################################

qualified_SNVs_even_flanking_depth_file = 'Qualified_SNVs_even_flanking_depth.txt'
qualified_SNVs_diff_flanking_depth_file = 'Qualified_SNVs_diff_flanking_depth.txt'

qualified_SNVs_even_flanking_depth_file_handle = open(qualified_SNVs_even_flanking_depth_file, 'w')
qualified_SNVs_diff_flanking_depth_file_handle = open(qualified_SNVs_diff_flanking_depth_file, 'w')

n_tst_total_unqualified = []
n_tst_each_unqualified = []
strand_bias_unqualified = []
qualified_SNVs = []
qualified_SNVs_even_flanking_depth = []
qualified_SNVs_diff_flanking_depth = []
total_num = 0
for each_snv in open(SNV_QC_nonsubsampled):

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

        # get SNVs with flanking depth difference higher than defined cutoff
        if (each_snv_n_tst_b >= min_var_reads_num) and (not ((each_snv_n_tst_fw < min_at_each_direction) or (each_snv_n_tst_bw < min_at_each_direction))) and (each_snv_strand_bias <= strand_bias_cutoff):
            qualified_SNVs.append(each_snv_key_with_sample_id)

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


# For report
print('The total number of detected SNVs: %s' % total_num)
print('The number of SNVs with less than %s reads harboring it: %s' % (min_var_reads_num, len(n_tst_total_unqualified)))
print('The number of SNVs with reads only from one direction: %s' % len(n_tst_each_unqualified))
print('The number of SNVs with strand bias higher than %s: %s' % (strand_bias_cutoff, len(strand_bias_unqualified)))
print('The number of qualified SNVs: %s' % len(qualified_SNVs))
print('The number of qualified SNVs with flanking depth (%sbp) difference higher than %s: %s' % (mean_depth_len, depth_difference_cutoff, len(qualified_SNVs_diff_flanking_depth)))
print('The number of qualified SNVs with flanking depth (%sbp) difference not higher than %s: %s' % (mean_depth_len, depth_difference_cutoff, len(qualified_SNVs_even_flanking_depth)))
print('Qualified SNVs with different mean flanking depth exported to: %s' % qualified_SNVs_diff_flanking_depth_file)
print('Qualified SNVs with similar mean flanking depth exported to: %s' % qualified_SNVs_even_flanking_depth_file)
