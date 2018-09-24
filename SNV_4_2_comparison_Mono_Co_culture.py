import os


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


def write_out(op_file, list):

    op_file_handle = open(op_file, 'w')

    for each in list:
        op_file_handle.write('%s\n' % each)

    op_file_handle.close()


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


# export SNVs to file
output_file_path = '/Users/songweizhi/Dropbox/Research/Flow_cell_datasets/Monoculture_VS_Coculture'

SNV_210_mono_uniq_file =   '%s/SNV_210_mono_uniq.txt'   % output_file_path
SNV_210_co_uniq_file =     '%s/SNV_210_co_uniq.txt'     % output_file_path
SNV_210_concurrence_file = '%s/SNV_210_concurrence.txt' % output_file_path
SNV_D2_mono_uniq_file =    '%s/SNV_D2_mono_uniq.txt'    % output_file_path
SNV_D2_co_uniq_file =      '%s/SNV_D2_co_uniq.txt'      % output_file_path
SNV_D2_concurrence_file =  '%s/SNV_D2_concurrence.txt'  % output_file_path


# write out
write_out(SNV_210_mono_uniq_file, SNV_210_mono_uniq)
write_out(SNV_210_co_uniq_file, SNV_210_co_uniq)
write_out(SNV_210_concurrence_file, SNV_210_concurrence)
write_out(SNV_D2_mono_uniq_file, SNV_D2_mono_uniq)
write_out(SNV_D2_co_uniq_file, SNV_D2_co_uniq)
write_out(SNV_D2_concurrence_file, SNV_D2_concurrence)


# run SNV_3_get_mutated_genes.py
get_mutated_genes_script = '/Users/songweizhi/PycharmProjects/FlowCellBiofilm/SNV_3_get_mutated_genes_v2.py'
get_mutated_genes_cmd_210_mono_uniq = 'python3 %s -fasta combined_ref.fasta -gff combined_ref.gff -ffn combined_ref.ffn -faa combined_ref.faa -SNV_matrix_cdc %s' % (get_mutated_genes_script, SNV_210_mono_uniq_file)
get_mutated_genes_cmd_210_co_uniq = 'python3 %s -fasta combined_ref.fasta -gff combined_ref.gff -ffn combined_ref.ffn -faa combined_ref.faa -SNV_matrix_cdc %s' % (get_mutated_genes_script, SNV_210_co_uniq_file)
get_mutated_genes_cmd_210_concurrence = 'python3 %s -fasta combined_ref.fasta -gff combined_ref.gff -ffn combined_ref.ffn -faa combined_ref.faa -SNV_matrix_cdc %s' % (get_mutated_genes_script, SNV_210_concurrence_file)
get_mutated_genes_cmd_D2_mono_uniq = 'python3 %s -fasta combined_ref.fasta -gff combined_ref.gff -ffn combined_ref.ffn -faa combined_ref.faa -SNV_matrix_cdc %s' % (get_mutated_genes_script, SNV_D2_mono_uniq_file)
get_mutated_genes_cmd_D2_co_uniq = 'python3 %s -fasta combined_ref.fasta -gff combined_ref.gff -ffn combined_ref.ffn -faa combined_ref.faa -SNV_matrix_cdc %s' % (get_mutated_genes_script, SNV_D2_co_uniq_file)
get_mutated_genes_cmd_D2_concurrence = 'python3 %s -fasta combined_ref.fasta -gff combined_ref.gff -ffn combined_ref.ffn -faa combined_ref.faa -SNV_matrix_cdc %s' % (get_mutated_genes_script, SNV_D2_concurrence_file)
print()

os.chdir(output_file_path)


# store annotation results in dicts
annotation_file = 'combined_ref_aa.faa.COG.arCOG.kegg'

gene_COG_id_dict = {}
gene_COG_function_dict = {}
for each_snv3 in open(annotation_file):
    each_snv3_split = each_snv3.strip().split('\t')
    gene_id = each_snv3_split[0]
    COG_cat = each_snv3_split[10]
    COG_id = each_snv3_split[8]
    COG_function = each_snv3_split[9]
    gene_COG_id_dict[gene_id] = COG_id
    gene_COG_function_dict[gene_id] = COG_function



# print('get_mutated_genes_cmd_210_mono_uniq')
# os.system(get_mutated_genes_cmd_210_mono_uniq)
# print('get_mutated_genes_cmd_210_co_uniq')
# os.system(get_mutated_genes_cmd_210_co_uniq)
# print('get_mutated_genes_cmd_210_concurrence')
# os.system(get_mutated_genes_cmd_210_concurrence)
# print('get_mutated_genes_cmd_D2_mono_uniq')
# os.system(get_mutated_genes_cmd_D2_mono_uniq)
# print('get_mutated_genes_cmd_D2_co_uniq')
# os.system(get_mutated_genes_cmd_D2_co_uniq)
# print('get_mutated_genes_cmd_D2_concurrence')
# os.system(get_mutated_genes_cmd_D2_concurrence)

all_affected_gene_list = []
gene_occurence_dict = {}
gene_occurence_count_dict = {}

strain_list = ['210', 'D2']

for each_gene in open('SNV_D2_co_uniq_summary.txt'):
    gene_id = each_gene.split('\t')[3]
    if gene_id != 'NA':

        if gene_id not in all_affected_gene_list:
            all_affected_gene_list.append(gene_id)

        if gene_id not in gene_occurence_count_dict:
            gene_occurence_count_dict[gene_id] = 1
        else:
            gene_occurence_count_dict[gene_id] += 1
        if gene_id not in gene_occurence_dict:
            gene_occurence_dict[gene_id] = ['co']

for each_gene in open('SNV_D2_mono_uniq_summary.txt'):
    gene_id = each_gene.split('\t')[3]
    if gene_id != 'NA':

        if gene_id not in all_affected_gene_list:
            all_affected_gene_list.append(gene_id)

        if gene_id not in gene_occurence_count_dict:
            gene_occurence_count_dict[gene_id] = 1
        else:
            gene_occurence_count_dict[gene_id] += 1
        if gene_id not in gene_occurence_dict:
            gene_occurence_dict[gene_id] = ['mono']
        else:
            if 'mono' not in gene_occurence_dict[gene_id]:
                gene_occurence_dict[gene_id].append('mono')

for each_gene in open('SNV_D2_concurrence_summary.txt'):
    gene_id = each_gene.split('\t')[3]
    if gene_id != 'NA':

        if gene_id not in all_affected_gene_list:
            all_affected_gene_list.append(gene_id)

        if gene_id not in gene_occurence_count_dict:
            gene_occurence_count_dict[gene_id] = 1
        else:
            gene_occurence_count_dict[gene_id] += 1
        if gene_id not in gene_occurence_dict:
            gene_occurence_dict[gene_id] = ['both']
        else:
            if 'both' not in gene_occurence_dict[gene_id]:
                gene_occurence_dict[gene_id].append('both')


for each_affected_gene in all_affected_gene_list:

    existence = ''
    existence_profile = gene_occurence_dict[each_affected_gene]
    #print(existence_profile)
    if len(existence_profile) == 1:
        if existence_profile == ['mono']:
            existence = 'mono'
        elif existence_profile == ['co']:
            existence = 'co'
        elif existence_profile == ['both']:
            existence = 'both'
    elif len(existence_profile) > 1:
        existence = 'both'

    current_COG_id2 = ''
    current_COG_function2 = ''
    if each_affected_gene in gene_COG_id_dict:
        current_COG_id2 = gene_COG_id_dict[each_affected_gene]
        current_COG_function2 = gene_COG_function_dict[each_affected_gene]
    else:
        current_COG_id2 = 'NA'
        current_COG_function2 = 'NA'

    for_report = '%s\t%s\t%s\t%s\t%s' % (each_affected_gene, existence, gene_occurence_count_dict[each_affected_gene], current_COG_id2, current_COG_function2)

    print(for_report)
















