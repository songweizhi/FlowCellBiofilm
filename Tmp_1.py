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

