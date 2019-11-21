
total = 0
for each in open('/Users/songweizhi/Desktop/4_1_deepSNV_subsampled.txt'):
    if not each.startswith('chr'):
        each_split = each.strip().split(',')
        p_value = float(each_split[4])

        total += 1
print(total)


total_st_005 = 0
total_lt_005 = 0

for each2 in open('/Users/songweizhi/Desktop/4_1_deepSNV_subsampled_sig0.9999999.txt'):
    if not each2.startswith('chr'):
        each_split2 = each2.strip().split(',')
        p_value2 = float(each_split2[4])

        if p_value2 > 0.05:
            total_lt_005 += 1
        else:
            total_st_005 += 1

print(total_st_005)
print(total_lt_005)








