
dict_a = {'A':'aa', 'B':'bb'}

dict_b = {'C':'cc', 'D':'dd'}






def merge_two_dict(dict_a, dict_b):

    dict_c = dict_a.copy()
    dict_c.update(dict_b)

    return dict_c


print(merge_two_dict(dict_a, dict_b))


