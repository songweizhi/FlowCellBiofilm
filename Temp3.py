

file_in = '/Users/songweizhi/Desktop/FC_000/SNV_QC_even_depth_matrix_D2_frequency_cdc_mutated_genes.txt'
file_out = '/Users/songweizhi/Desktop/FC_000/SNV_QC_even_depth_matrix_D2_frequency_cdc_mutated_genes_category.txt'

def get_mutation_cate_summary(file_in, file_out):

    mutation_cate_dict = {}
    for each_snv in open(file_in):
        each_snv_split = each_snv.strip().split('\t')
        mutated_gene = each_snv_split[3]
        mutation_cate = each_snv_split[5]

        if mutated_gene != 'NA':
            if mutated_gene not in mutation_cate_dict:
                mutation_cate_dict[mutated_gene] = [mutation_cate]
            else:
                mutation_cate_dict[mutated_gene].append(mutation_cate)

    mutation_cate_list = ['Missense', 'Nonsense', 'Silent', 'Fragment_deletion', 'Frameshift']

    file_out_handle = open(file_out, 'w')
    file_out_handle.write('\tMis\tNon\tSilen\tFD\tFS\n')
    for each_gene in mutation_cate_dict:

        mutation_cate_occurence = []
        for each_cate in mutation_cate_list:
            mutation_cate_occurence.append(str(mutation_cate_dict[each_gene].count(each_cate)))
        file_out_handle.write('%s\t%s\n' % (each_gene, '\t'.join(mutation_cate_occurence)))

    file_out_handle.close()
