import os
import glob
import shutil
import argparse


############################################### input file and parameters ##############################################

parser = argparse.ArgumentParser(description='', add_help=False)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')

optional.add_argument('-h', action='help', help='Show this help message and exit')
required.add_argument('-wd', dest='SNV_QC_WD', nargs='?', required=True,  type=str, help='deepSNV QC wd')
required.add_argument('-min_both', dest='MIN_BOTH', nargs='?', required=True, type=int, help='The minimum number of reads harboring SNV')
required.add_argument('-min_each', dest='MIN_EACH', nargs='?', required=True, type=int, help='The minimum number of reads harboring SNV at each direction')
required.add_argument('-strand_bias', dest='STRAND_BIAS', nargs='?', required=True, type=int, help='strand_bias cutoff')
required.add_argument('-depth_diff', dest='DEPTH_DIFF', nargs='?', required=True,  type=int, help='depth difference cutoff')
required.add_argument('-deplen', dest='DEPLEN', nargs='?', required=True, type=int, help='flanking length for mean depth calculation')
required.add_argument('-sep_plot', dest='SEP_PLOT', nargs='?', required=False, type=int, help='separate depth plot with provide cutoff')

args = vars(parser.parse_args())
SNV_qc_wd = args['SNV_QC_WD']
min_var_reads_num = args['MIN_BOTH']
min_at_each_direction = args['MIN_EACH']
strand_bias_cutoff = args['STRAND_BIAS']
depth_difference_cutoff = args['DEPTH_DIFF']
mean_depth_len = args['DEPLEN']
separate_plot = args['SEP_PLOT']


# input QC txt
QC_txt = 'SNV_QC.txt'
pwd_QC_txt = '%s/%s' % (SNV_qc_wd, QC_txt)


######################################################### Main #########################################################

qualified_SNVs_even_flanking_depth_file = '%s/Qualified_SNVs_even_flanking_depth.txt' % SNV_qc_wd
qualified_SNVs_diff_flanking_depth_file = '%s/Qualified_SNVs_diff_flanking_depth.txt' % SNV_qc_wd

qualified_SNVs_even_flanking_depth_file_handle = open(qualified_SNVs_even_flanking_depth_file, 'w')
qualified_SNVs_diff_flanking_depth_file_handle = open(qualified_SNVs_diff_flanking_depth_file, 'w')

n_tst_total_unqualified = []
n_tst_each_unqualified = []
strand_bias_unqualified = []
qualified_SNVs = []
qualified_SNVs_even_flanking_depth = []
qualified_SNVs_diff_flanking_depth = []

plot_prefix_even = []
plot_prefix_diff = []

total_num = 0
for each_snv in open(pwd_QC_txt):

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
print('The total number of detected SNVs: %s' % total_num)
print('The number of SNVs with less than %s reads harboring it: %s' % (min_var_reads_num, len(n_tst_total_unqualified)))
print('The number of SNVs with reads only from one direction: %s' % len(n_tst_each_unqualified))
print('The number of SNVs with strand bias higher than %s: %s' % (strand_bias_cutoff, len(strand_bias_unqualified)))
print('The number of qualified SNVs: %s' % len(qualified_SNVs))
print('The number of qualified SNVs with flanking depth (%sbp) difference higher than %s: %s' % (mean_depth_len, depth_difference_cutoff, len(qualified_SNVs_diff_flanking_depth)))
print('The number of qualified SNVs with flanking depth (%sbp) difference not higher than %s: %s' % (mean_depth_len, depth_difference_cutoff, len(qualified_SNVs_even_flanking_depth)))
print('Qualified SNVs with different mean flanking depth exported to: %s' % qualified_SNVs_diff_flanking_depth_file)
print('Qualified SNVs with similar mean flanking depth exported to: %s' % qualified_SNVs_even_flanking_depth_file)


################################ separate plots according to provided difference cutoff ################################

# separate plot files
if separate_plot ==1:

    plot_folder = 'SNV_depth_plot'
    plot_folder_even = '%s_even_%s' % (plot_folder, depth_difference_cutoff)
    plot_folder_diff = '%s_diff_%s' % (plot_folder, depth_difference_cutoff)
    plot_folder_ignored = '%s_ignored_%s' % (plot_folder, depth_difference_cutoff)

    pwd_plot_folder = '%s/%s' % (SNV_qc_wd, plot_folder)
    pwd_plot_folder_even = '%s/%s' % (SNV_qc_wd, plot_folder_even)
    pwd_plot_folder_diff = '%s/%s' % (SNV_qc_wd, plot_folder_diff)
    pwd_plot_folder_ignored = '%s/%s' % (SNV_qc_wd, plot_folder_ignored)

    # prepare folder
    if os.path.isdir(pwd_plot_folder_even):
        shutil.rmtree(pwd_plot_folder_even)
        os.mkdir(pwd_plot_folder_even)
    else:
        os.mkdir(pwd_plot_folder_even)

    if os.path.isdir(pwd_plot_folder_diff):
        shutil.rmtree(pwd_plot_folder_diff)
        os.mkdir(pwd_plot_folder_diff)
    else:
        os.mkdir(pwd_plot_folder_diff)

    if os.path.isdir(pwd_plot_folder_ignored):
        shutil.rmtree(pwd_plot_folder_ignored)
        os.mkdir(pwd_plot_folder_ignored)
    else:
        os.mkdir(pwd_plot_folder_ignored)


    # get plots with similar flanking deoth
    for each_even in plot_prefix_even:
        pwd_plot = '%s/%s*' % (pwd_plot_folder, each_even)
        cmd = 'cp %s %s/' % (pwd_plot, pwd_plot_folder_even)
        os.system(cmd)

    # get plots with similar flanking deoth
    for each_diff in plot_prefix_diff:
        pwd_plot = '%s/%s*' % (pwd_plot_folder, each_diff)
        cmd = 'cp %s %s/' % (pwd_plot, pwd_plot_folder_diff)
        os.system(cmd)

