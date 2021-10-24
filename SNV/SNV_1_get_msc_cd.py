import numpy as np


def get_msc(mutant_freq_start, mutant_freq_end, generation_num):
    malthusian_selection_coefficients = np.log((mutant_freq_end / mutant_freq_start) / ((1 - mutant_freq_end) / (1 - mutant_freq_start))) / generation_num
    malthusian_selection_coefficients = float("{0:.2f}".format(malthusian_selection_coefficients))
    return malthusian_selection_coefficients


def get_msc_list(freq_list, iRep_list):

    if freq_list[0] == 0:
        freq_list_msc_1 = '-'
    else:
        freq_list_msc_1 = get_msc((1 / 662), freq_list[0], (iRep_list[0] + iRep_list[1]) / 2)

    if freq_list[1] == 0:
        freq_list_msc_2 = '-'
    else:
        if freq_list[0] == 0:
            freq_list_msc_2 = get_msc((1 / 662), freq_list[1], (iRep_list[1] + iRep_list[2]) / 2)
        else:
            freq_list_msc_2 = get_msc(freq_list[0], freq_list[1], (iRep_list[1] + iRep_list[2]) / 2)

    if freq_list[2] == 0:
        freq_list_msc_3 = '-'
    else:
        if freq_list[1] == 0:
            freq_list_msc_3 = get_msc((1 / 662), freq_list[2], (iRep_list[2] + iRep_list[3]) / 2)
        else:
            freq_list_msc_3 = get_msc(freq_list[1], freq_list[2], (iRep_list[2] + iRep_list[3]) / 2)

    if freq_list[3] == 0:
        freq_list_msc_4 = '-'
    else:
        if freq_list[2] == 0:
            freq_list_msc_4 = get_msc((1 / 662), freq_list[3], ((iRep_list[3] + iRep_list[4]) / 2)*(15/9))
        else:
            freq_list_msc_4 = get_msc(freq_list[2], freq_list[3], ((iRep_list[3] + iRep_list[4]) / 2)*(15/9))

    msc_list = [freq_list_msc_1, freq_list_msc_2, freq_list_msc_3, freq_list_msc_4]

    return msc_list


wd = '/Users/songweizhi/FC_Biofilm/OneStep_MinBoth_10_MinEach_1_StrandBias_10_DepthDiff_15'
strain       = 'D2'  # 210 or D2
iRep_values  = '/Users/songweizhi/FC_Biofilm/iRep_values.txt'
matrix       = '%s/SNV_QC_ncd_even_flk_depth_%s_matrix.txt'     % (wd, strain)
msc_file_out = '%s/SNV_QC_ncd_even_flk_depth_%s_matrix_msc.txt' % (wd, strain)


iRep_value_dict = {}
for iRep_value in open(iRep_values):
    iRep_value_split = iRep_value.strip().split('\t')
    iRep_value_dict[iRep_value_split[0]] = float(iRep_value_split[1])


msc_file_out_handle = open(msc_file_out, 'w')
header_line = []
for snv_freq in open(matrix):
    snv_freq_split = snv_freq.strip().split('\t')
    if snv_freq.startswith('	'):
        header_line = snv_freq_split
        msc_file_out_handle.write(snv_freq)
    else:
        snv_id = snv_freq_split[0]
        mono_A = [float(i) for i in snv_freq_split[1:5]]
        mono_B = [float(i) for i in snv_freq_split[5:9]]
        mono_C = [float(i) for i in snv_freq_split[9:13]]
        co_A   = [float(i) for i in snv_freq_split[13:17]]
        co_B   = [float(i) for i in snv_freq_split[17:21]]
        co_C   = [float(i) for i in snv_freq_split[21:25]]

        if snv_id.startswith('2.10'):
            mono_A_iRep = [iRep_value_dict['210D0'], iRep_value_dict['Mono_210_D9_A'], iRep_value_dict['Mono_210_D18_A'], iRep_value_dict['Mono_210_D27_A'], iRep_value_dict['Mono_210_D42_A']]
            mono_B_iRep = [iRep_value_dict['210D0'], iRep_value_dict['Mono_210_D9_B'], iRep_value_dict['Mono_210_D18_B'], iRep_value_dict['Mono_210_D27_B'], iRep_value_dict['Mono_210_D42_B']]
            mono_C_iRep = [iRep_value_dict['210D0'], iRep_value_dict['Mono_210_D9_C'], iRep_value_dict['Mono_210_D18_C'], iRep_value_dict['Mono_210_D27_C'], iRep_value_dict['Mono_210_D42_C']]
            co_A_iRep   = [iRep_value_dict['210D0'], iRep_value_dict['Co_210_D9_A'], iRep_value_dict['Co_210_D18_A'], iRep_value_dict['Co_210_D27_A'], iRep_value_dict['Co_210_D42_A']]
            co_B_iRep   = [iRep_value_dict['210D0'], iRep_value_dict['Co_210_D9_B'], iRep_value_dict['Co_210_D18_B'], iRep_value_dict['Co_210_D27_B'], iRep_value_dict['Co_210_D42_B']]
            co_C_iRep   = [iRep_value_dict['210D0'], iRep_value_dict['Co_210_D9_C'], iRep_value_dict['Co_210_D18_C'], iRep_value_dict['Co_210_D27_C'], iRep_value_dict['Co_210_D42_C']]
        else:
            mono_A_iRep = [iRep_value_dict['D2D0'], iRep_value_dict['Mono_D2_D9_A'], iRep_value_dict['Mono_D2_D18_A'], iRep_value_dict['Mono_D2_D27_A'], iRep_value_dict['Mono_D2_D42_A']]
            mono_B_iRep = [iRep_value_dict['D2D0'], iRep_value_dict['Mono_D2_D9_B'], iRep_value_dict['Mono_D2_D18_B'], iRep_value_dict['Mono_D2_D27_B'], iRep_value_dict['Mono_D2_D42_B']]
            mono_C_iRep = [iRep_value_dict['D2D0'], iRep_value_dict['Mono_D2_D9_C'], iRep_value_dict['Mono_D2_D18_C'], iRep_value_dict['Mono_D2_D27_C'], iRep_value_dict['Mono_D2_D42_C']]
            co_A_iRep   = [iRep_value_dict['D2D0'], iRep_value_dict['Co_D2_D9_A'], iRep_value_dict['Co_D2_D18_A'], iRep_value_dict['Co_D2_D27_A'], iRep_value_dict['Co_D2_D42_A']]
            co_B_iRep   = [iRep_value_dict['D2D0'], iRep_value_dict['Co_D2_D9_B'], iRep_value_dict['Co_D2_D18_B'], iRep_value_dict['Co_D2_D27_B'], iRep_value_dict['Co_D2_D42_B']]
            co_C_iRep   = [iRep_value_dict['D2D0'], iRep_value_dict['Co_D2_D9_C'], iRep_value_dict['Co_D2_D18_C'], iRep_value_dict['Co_D2_D27_C'], iRep_value_dict['Co_D2_D42_C']]

        mono_A_msc = get_msc_list(mono_A, mono_A_iRep)
        mono_B_msc = get_msc_list(mono_B, mono_B_iRep)
        mono_C_msc = get_msc_list(mono_C, mono_C_iRep)
        co_A_msc   = get_msc_list(co_A, co_A_iRep)
        co_B_msc   = get_msc_list(co_B, co_B_iRep)
        co_C_msc   = get_msc_list(co_C, co_C_iRep)

        # print('mono_A\t%s\t%s\t%s' % (mono_A, mono_A_iRep, mono_A_msc))
        # print('mono_B\t%s\t%s\t%s' % (mono_B, mono_B_iRep, mono_B_msc))
        # print('mono_C\t%s\t%s\t%s' % (mono_C, mono_C_iRep, mono_C_msc))
        # print('co_A\t%s\t%s\t%s'   % (co_A, co_A_iRep, co_A_msc))
        # print('co_B\t%s\t%s\t%s'   % (co_B, co_B_iRep, co_B_msc))
        # print('co_C\t%s\t%s\t%s'   % (co_C, co_C_iRep, co_C_msc))

        snv_msc = [snv_id] + mono_A_msc + mono_B_msc + mono_C_msc + co_A_msc + co_B_msc + co_C_msc
        snv_msc_str = [str(i) for i in snv_msc]
        msc_file_out_handle.write('%s\n' % '\t'.join(snv_msc_str))

msc_file_out_handle.close()

