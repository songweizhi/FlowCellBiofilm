import os

dict = {}

sample_list = []
for each in open('/Users/songweizhi/Desktop/aaa.txt'):
    print(each)
    dict[each.strip().split('\t')[0]] = each.strip().split('\t')[1]
    sample_list.append(each.strip().split('\t')[1])
print(dict)

print(sample_list)