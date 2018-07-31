
file_in = '/Users/songweizhi/Desktop/FC/SNV_QC_even_depth_matrix_D2_frequency_cdc_summary.txt'

for each_snv in open(file_in):
    each_snv_split = each_snv.strip().split('\t')
    effects = each_snv_split[6]
    effect_list = effects.split('|')

    if len(effect_list) == 1:
        pass
        #print(effect_list)

    else:
        #pass
        print(effect_list)



































