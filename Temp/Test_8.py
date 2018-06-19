import os


wd = '/Users/songweizhi/Desktop'
deepSNV_output = 'deepSNV_output_summary_all_existence.txt'
deepSNV_output_cdc = 'deepSNV_output_summary_all_existence_cdc.txt'

os.chdir(wd)


num_list = []
for each in open(deepSNV_output_cdc):

    m = each.strip().split('\t')[0].split('|')[1]

    if '-' in m:
        print(m)
        m_l = int(m.split('-')[0])
        m_r = int(m.split('-')[1])
        m_list = list(range(m_l, m_r + 1))
        print(m_list)
        for each_m in m_list:
            num_list.append(int(each_m))


    else:
        num_list.append(int(m))




num_list2 = []
for each2 in open(deepSNV_output):
    if not each2.startswith('\t'):
        m2 = each2.strip().split('\t')[0].split('|')[1]
        num_list2.append(int(m2))


print(num_list)
print(num_list2)


print(len(num_list))
print(len(num_list2))

