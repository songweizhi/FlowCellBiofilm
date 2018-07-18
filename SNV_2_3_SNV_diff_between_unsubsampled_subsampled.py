import os


def get_uniq_snv_list(input_QC_file):
    uniq_SNV_list = []
    for each_snv_unsubsampled in open(input_QC_file):
        if not each_snv_unsubsampled.startswith('sample'):
            each_snv_split = each_snv_unsubsampled.strip().split(',')
            each_snv_chr = each_snv_split[1]
            each_snv_pos = int(each_snv_split[2])
            each_snv_ref = each_snv_split[3]
            each_snv_var = each_snv_split[4]
            each_snv_key_no_sample_id = '%s|%s|%s|%s' % (each_snv_chr, each_snv_pos, each_snv_ref, each_snv_var)
            if each_snv_key_no_sample_id not in uniq_SNV_list:
                uniq_SNV_list.append(each_snv_key_no_sample_id)

    return uniq_SNV_list


# input files
wd = '/Users/songweizhi/Dropbox/Research/Flow_cell/SNV_2_3_SNV_diff_between_unsubsampled_subsampled'
deepSNV_output_unsubsampled = 'deepSNV_output_combined_QC.txt'
deepSNV_output_subsampled = 'deepSNV_output_combined_QC_subsampled.txt'

os.chdir(wd)

# get SNV list
uniq_SNV_list_unsubsampled = get_uniq_snv_list(deepSNV_output_unsubsampled)
uniq_SNV_list_subsampled = get_uniq_snv_list(deepSNV_output_subsampled)

# get shared and uniq SNVs
shared_snv_list = []
unsubsampled_uniq_snv_list = []
for each_snv in uniq_SNV_list_unsubsampled:
    if each_snv in uniq_SNV_list_subsampled:
        shared_snv_list.append(each_snv)
    else:
        unsubsampled_uniq_snv_list.append(each_snv)

subsampled_uniq_snv_list = []
for each_snv_2 in uniq_SNV_list_subsampled:
    if each_snv_2 not in uniq_SNV_list_unsubsampled:
        subsampled_uniq_snv_list.append(each_snv_2)

# report
print('The number of SNVs detected from unsubsampled samples: %s' % len(uniq_SNV_list_unsubsampled))
print('The number of SNVs detected from subsampled samples: %s' % len(uniq_SNV_list_subsampled))
print('The number of shared SNVs: %s' % len(shared_snv_list))
print('The number of SNVs uniq to unsubsampled samples: %s' % len(unsubsampled_uniq_snv_list))
print('The number of SNVs uniq to subsampled samples: %s' % len(subsampled_uniq_snv_list))

