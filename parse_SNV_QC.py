import os


########################################## specify input files and parameters ##########################################

wd = '/Users/songweizhi/Desktop/666666'
SNV_quality_file = 'deepSNV_output_combined_QC.txt'
min_var_reads_num = 10
min_at_each_direction = 1
strand_bias_cutoff = 20
depth_difference_cutoff = 15
mean_depth_len = 2000
os.chdir(wd)


############################################ parse quality of detected SNVs ############################################

n_tst_total_unqualified = []
n_tst_each_unqualified = []
strand_bias_unqualified = []
depth_diff_unqualified = []
qualified_SNVs = []
total_num = 0
for each_snv in open(SNV_quality_file):
    each_snv_split = each_snv.strip().split(',')
    if not each_snv.startswith('sample'):
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

        # get qualified SNV list
        if (each_snv_n_tst_b >= min_var_reads_num) and (each_snv_n_tst_fw >= min_at_each_direction) and (each_snv_n_tst_bw >= min_at_each_direction) and (each_snv_strand_bias <= strand_bias_cutoff) and (each_snv_mean_depth_diff <= depth_difference_cutoff):
            qualified_SNVs.append(each_snv_key_with_sample_id)

        total_num += 1

# report
print('The total number of detected SNVs: %s' % total_num)
print('The number of SNVs with less than %s reads harboring it: %s' % (min_var_reads_num, len(n_tst_total_unqualified)))
print('The number of SNVs with reads only from one direction: %s' % len(n_tst_each_unqualified))
print('The number of SNVs with strand bias higher than %s%s: %s' % (strand_bias_cutoff, '%',len(strand_bias_unqualified)))
print('The number of SNVs with flanking depth (%sbp) difference higher than %s%s: %s' % (mean_depth_len, depth_difference_cutoff, '%', len(depth_diff_unqualified)))
print('The number of qualified SNVs: %s' % len(qualified_SNVs))


################################### get matrix for SNVs with similar flanking depth ####################################





################################## get matrix for SNVs with different flanking depth ###################################














