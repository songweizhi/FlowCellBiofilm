ref_length_dict = {}
pwd_deepSNV_output_combined = ''
depth_file_folder = ''
depth_flanking_length = 10
print('Plotting depth for detected SNVs')
pwd_SNV_depth_plot = ''
def plot_sam_depth():
    pass



for each_snv in open(pwd_deepSNV_output_combined):
    each_snv_split = each_snv.strip().split(',')
    each_snv_sample = each_snv_split[0]
    each_snv_chr = each_snv_split[1]
    each_snv_pos = int(each_snv_split[2])
    each_snv_ref = each_snv_split[3]
    each_snv_var = each_snv_split[4]
    each_snv_p_val = float("{0:.3f}".format(float(each_snv_split[5])))
    each_snv_freq = float("{0:.3f}".format(float(each_snv_split[6])))
    n_tst_fw = int(each_snv_split[8])
    cov_tst_fw = int(each_snv_split[9])
    n_tst_bw = int(each_snv_split[10])
    cov_tst_bw = int(each_snv_split[11])

    each_snv_key_no_sample_id = '%s|%s|%s|%s' % (each_snv_chr, each_snv_pos, each_snv_ref, each_snv_var)
    each_snv_key_with_sample_id = '%s|%s|%s|%s|%s' % (each_snv_sample, each_snv_chr, each_snv_pos, each_snv_ref, each_snv_var)

    # plot depth
    pwd_depth_file = '%s/%s.depth' % (depth_file_folder, each_snv_sample)

    # get start position to plot
    plot_start = each_snv_pos - depth_flanking_length
    if plot_start <= 0:
        plot_start = 1

    # get end position to plot
    plot_end = each_snv_pos + depth_flanking_length
    if plot_end >= ref_length_dict[each_snv_chr]:
        plot_end = ref_length_dict[each_snv_chr]

    # png file name
    each_snv_id_modified = '_'.join(each_snv_key_no_sample_id.split('|'))
    plot_file_name = '%s_%s_f%sbp_%smer_md%sbp_%s' % (each_snv_id_modified, each_snv_sample, depth_flanking_length, depth_kmer, flanking_length_for_mean_depth, snv_flanking_depth_diff_dict[each_snv_key_with_sample_id])
    pwd_plot_file_name = '%s/%s' % (pwd_SNV_depth_plot, plot_file_name)

    # plot depth
    plot_sam_depth(pwd_depth_file, each_snv_chr, plot_start, plot_end, depth_kmer, str(each_snv_pos), pwd_plot_file_name)
