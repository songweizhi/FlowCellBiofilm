

gene_list_1 = []
for each_1 in open('/Users/songweizhi/Desktop/FC/deepSNV_mutated_gene_copy.txt'):
    gene = each_1.strip().split('\t')[3]
    if gene not in gene_list_1:
        gene_list_1.append(gene)

print(gene_list_1)

print(len(gene_list_1))


gene_list_2 = []
for each_2 in open('/Users/songweizhi/Desktop/FC/deepSNV_mutated_gene_sorted.txt'):
    gene2 = each_2.strip().split('\t')[3]
    if gene2 not in gene_list_2:
        gene_list_2.append(gene2)

print(gene_list_2)
print(len(gene_list_2))

print(gene_list_1 == gene_list_2)



print(len(set(gene_list_1).intersection(gene_list_2)))
