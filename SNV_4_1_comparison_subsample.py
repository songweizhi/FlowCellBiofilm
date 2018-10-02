
def unique_list_elements(list_input):

    list_output = []
    for each_element in list_input:
        if each_element not in list_output:
            list_output.append(each_element)

    return list_output


##################################################### input files ######################################################

# deepSNV output (nonsubsampled: /srv/scratch/z5039045/Flow_cell_biofilm/4_2_SNV_QC/output_f50000bp_1000mer_dl2000bp/SNV_combined.txt)
SNVs_no_subsample = '/Users/songweizhi/Dropbox/Research/Flow_cell_datasets/Effect_of_subsample/SNV_combined_nonsubsample.txt'

# deepSNV output (subsampled: /srv/scratch/z5039045/Flow_cell_biofilm/4_2_SNV_QC_subsampled/output_f50000bp_1000mer_dl2000bp/SNV_combined.txt)
SNVs_subsampled = '/Users/songweizhi/Dropbox/Research/Flow_cell_datasets/Effect_of_subsample/SNV_combined_subsample.txt'


# get nonsubsampled SNV list
SNV_list_no_subsample = []
for each_SNV_1 in open(SNVs_no_subsample):
    each_SNV_1_split = each_SNV_1.strip().split(',')
    each_SNV_1_id = '%s|%s|%s|%s|%s' % (each_SNV_1_split[0], each_SNV_1_split[1], each_SNV_1_split[2], each_SNV_1_split[3], each_SNV_1_split[4])
    SNV_list_no_subsample.append(each_SNV_1_id)


# get subsampled SNV list
SNV_list_subsampled = []
for each_SNV_2 in open(SNVs_subsampled):
    each_SNV_2_split = each_SNV_2.strip().split(',')
    each_SNV_2_id = '%s|%s|%s|%s|%s' % (each_SNV_2_split[0], each_SNV_2_split[1], each_SNV_2_split[2], each_SNV_2_split[3], each_SNV_2_split[4])
    SNV_list_subsampled.append(each_SNV_2_id)


# get shared SNV
SNVs_shared = sorted(set(SNV_list_no_subsample).intersection(SNV_list_subsampled))


# get SNVs unique to nonsubsampled
SNVs_unique_to_unsubsampled = []
for each_1 in SNV_list_no_subsample:
    if each_1 not in SNVs_shared:
        SNVs_unique_to_unsubsampled.append(each_1)


# get SNVs unique to subsampled
SNVs_unique_to_subsampled = []
for each_2 in SNV_list_subsampled:
    if each_2 not in SNVs_shared:
        SNVs_unique_to_subsampled.append(each_2)


# for report
print('SNV_list_nonsubsampled: %s' % len(SNV_list_no_subsample))
print('SNV_list_subsampled: %s' % len(SNV_list_subsampled))
print('SNVs_shared: %s' % len(SNVs_shared))
print('SNVs_unique_to_nonsubsampled: %s' % len(SNVs_unique_to_unsubsampled))
print('SNVs_unique_to_subsampled: %s' % len(SNVs_unique_to_subsampled))

