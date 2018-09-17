
def unique_list_elements(list_input):

    list_output = []
    for each_element in list_input:
        if each_element not in list_output:
            list_output.append(each_element)

    return list_output


def remove_l2_elements_from_l1(l1, l2):

    l1_new = []
    for each in l1:
        if each not in l2:
            l1_new.append(each)
    return l1_new


##################################################### input files ######################################################

# deepSNV output (nonsubsampled: /srv/scratch/z5039045/Flow_cell_biofilm/4_2_SNV_QC/output_f50000bp_1000mer_dl2000bp/SNV_QC.txt)
#SNVs_QC_nonsubsampled = '/Users/songweizhi/Dropbox/Research/Flow_cell_datasets/Effect_of_coculture/SNV_QC.txt'
SNVs_QC_nonsubsampled = '/Users/songweizhi/Dropbox/Research/Flow_cell_datasets/parse_SNV_QC_v2/Qualified_SNVs_even_flanking_depth.txt'
#SNVs_QC_nonsubsampled = '/Users/songweizhi/Dropbox/Research/Flow_cell_datasets/parse_SNV_QC_v2/Qualified_SNVs_diff_flanking_depth.txt'


mono_210_SNV_list = []
co_210_SNV_list = []
mono_D2_SNV_list = []
co_D2_SNV_list = []
total_SNV = 0
for each_SNV in open(SNVs_QC_nonsubsampled):
    each_SNV_split = each_SNV.strip().split(',')
    treatment_id = each_SNV_split[0].split('D')[0]
    strain = each_SNV_split[1].split('_')[0]
    SNV_id = '%s|%s|%s|%s' % (each_SNV_split[1], each_SNV_split[2], each_SNV_split[3], each_SNV_split[4])

    if (strain == '2.10') and (treatment_id in ['1', '5', '9']):
        mono_210_SNV_list.append(SNV_id)

    elif (strain == '2.10') and (treatment_id in ['4', '8', '12']):
        co_210_SNV_list.append(SNV_id)

    elif (strain == 'D2') and (treatment_id in ['2', '6', '10']):
        mono_D2_SNV_list.append(SNV_id)

    elif (strain == 'D2') and (treatment_id in ['4', '8', '12']):
        co_D2_SNV_list.append(SNV_id)

    total_SNV += 1


# uniq list element
mono_210_SNV_list = unique_list_elements(mono_210_SNV_list)
co_210_SNV_list = unique_list_elements(co_210_SNV_list)
mono_D2_SNV_list = unique_list_elements(mono_D2_SNV_list)
co_D2_SNV_list = unique_list_elements(co_D2_SNV_list)


# 210
SNV_210_concurrence = set(mono_210_SNV_list).intersection(co_210_SNV_list)
SNV_210_mono_uniq = remove_l2_elements_from_l1(mono_210_SNV_list, SNV_210_concurrence)
SNV_210_co_uniq = remove_l2_elements_from_l1(co_210_SNV_list, SNV_210_concurrence)

# D2
SNV_D2_concurrence = set(mono_D2_SNV_list).intersection(co_D2_SNV_list)
SNV_D2_mono_uniq = remove_l2_elements_from_l1(mono_D2_SNV_list, SNV_D2_concurrence)
SNV_D2_co_uniq = remove_l2_elements_from_l1(co_D2_SNV_list, SNV_D2_concurrence)

# get total
total_210 = len(SNV_210_concurrence) + len(SNV_210_mono_uniq) + len(SNV_210_co_uniq)
total_D2 = len(SNV_D2_concurrence) + len(SNV_D2_mono_uniq) + len(SNV_D2_co_uniq)

# For report
print('Qualified SNVs in total: %s' % total_SNV)
print('Qualified 2.10 SNVs: %s' % total_210)
print('Qualified D2 SNVs: %s' % total_D2)
print()
print('Qualified 2.10 SNVs uniq to monoculture: %s (%s' % (len(SNV_210_mono_uniq), float("{0:.1f}".format(len(SNV_210_mono_uniq)/total_210 * 100))) + '%)')
print('Qualified 2.10 SNVs uniq to coculture: %s (%s' % (len(SNV_210_co_uniq), float("{0:.1f}".format(len(SNV_210_co_uniq)/total_210 * 100))) + '%)')
print('Qualified 2.10 SNVs concurrent in both: %s (%s' % (len(SNV_210_concurrence), float("{0:.1f}".format(len(SNV_210_concurrence)/total_210 * 100))) + '%)')
print()
print('Qualified D2 SNVs uniq to monoculture: %s (%s' % (len(SNV_D2_mono_uniq), float("{0:.1f}".format(len(SNV_D2_mono_uniq)/total_D2 * 100))) + '%)')
print('Qualified D2 SNVs uniq to coculture: %s (%s' % (len(SNV_D2_co_uniq), float("{0:.1f}".format(len(SNV_D2_co_uniq)/total_D2 * 100))) + '%)')
print('Qualified D2 SNVs concurrent in both: %s (%s' % (len(SNV_D2_concurrence), float("{0:.1f}".format(len(SNV_D2_concurrence)/total_D2 * 100))) + '%)')


