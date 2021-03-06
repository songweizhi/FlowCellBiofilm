import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_ext = os.path.splitext(file_name)

    return file_path, file_basename, file_ext


def get_rf_pos_list(gene_len):

    rf_pos_list = []
    n = 0
    while n < (gene_len / 3):
        start_pos = 1 + 3 * n
        rf_pos_list.append([start_pos, start_pos + 1, start_pos + 2])
        n += 1

    return rf_pos_list


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
    file_out_handle.write('Gene\tMis\tNon\tSilen\tFD\tFS\n')
    for each_gene in mutation_cate_dict:

        mutation_cate_occurence = []
        for each_cate in mutation_cate_list:
            mutation_cate_occurence.append(str(mutation_cate_dict[each_gene].count(each_cate)))
        file_out_handle.write('%s\t%s\n' % (each_gene, '\t'.join(mutation_cate_occurence)))

    file_out_handle.close()

    return mutation_cate_dict


def check_existence(num_list_1, num_list_2):

    existence_list = []
    for each_num in num_list_1:
        if each_num in num_list_2:
            existence_list.append(1)
        else:
            existence_list.append(0)

    if 1 in existence_list:
        return True
    else:
        return False


def check_parallel(effects):

    effect_list = effects.split('|')

    treatment_id_list = []
    for each_effect in effect_list:
        each_effect_treatment_id = int(each_effect.split('_')[0])
        treatment_id_list.append(each_effect_treatment_id)

    SNV_type = ''

    if (len(treatment_id_list) == 1) and ((check_existence(treatment_id_list, [1, 5, 9]) is True) or (check_existence(treatment_id_list, [2, 6, 10]) is True)):
        SNV_type = 'SM'

    elif (len(treatment_id_list) == 1) and (treatment_id_list[0] in [4, 8, 12]):
        SNV_type = 'SC'

    elif (len(treatment_id_list) > 1) and ((check_existence(treatment_id_list, [1, 5, 9]) is True) or (check_existence(treatment_id_list, [2, 6, 10]) is True)) and (check_existence(treatment_id_list, [4, 8, 12]) is False):
        SNV_type = 'PM'

    elif (len(treatment_id_list) > 1) and (check_existence(treatment_id_list, [1, 5, 9, 2, 6, 10]) is False) and (check_existence(treatment_id_list, [4, 8, 12]) is True):
        SNV_type = 'PC'

    elif (len(treatment_id_list) > 1) and ((check_existence(treatment_id_list, [1, 5, 9]) is True) or (check_existence(treatment_id_list, [2, 6, 10]) is True)) and (check_existence(treatment_id_list, [4, 8, 12]) is True):
        SNV_type = 'PMC'

    return SNV_type


parser = argparse.ArgumentParser()
parser.add_argument('-SNV_matrix_cdc', required=True, help='SNV_matrix_cdc file')
parser.add_argument('-fasta', required=False, help='fasta file')
parser.add_argument('-gff', required=False, help='gff file')
parser.add_argument('-ffn', required=False, help='ffn file')
parser.add_argument('-faa', required=False, help='faa file')

args = vars(parser.parse_args())
SNV_matrix_cdc = args['SNV_matrix_cdc']
combined_ref_fasta = args['fasta']
combined_ref_gff = args['gff']
combined_ref_ffn = args['ffn']
combined_ref_faa = args['faa']

annotation_file = 'combined_ref_aa.faa.COG.arCOG.kegg'
#effect_file = 'SNV_QC_diff_depth_matrix_210_frequency_cdc_mutation_effect.txt'
transl_table = 11


# output files
SNV_matrix_cdc_path, SNV_matrix_cdc_basename, SNV_matrix_cdc_ext = sep_path_basename_ext(SNV_matrix_cdc)
output_mutated_genes = '%s_mutated_genes.txt' % SNV_matrix_cdc_basename
output_mutated_genes_cate = '%s_mutated_genes_category.txt' % SNV_matrix_cdc_basename
output_mutated_genes_cate_fun = '%s_mutated_genes_cate_fun.txt' % SNV_matrix_cdc_basename

output_summary = '%s_summary.txt' % SNV_matrix_cdc_basename
effect_file = '%s_mutation_effect.txt' % SNV_matrix_cdc_basename

output_seq_nc = '%s_mutated_genes_nc.fasta' % SNV_matrix_cdc_basename
output_seq_aa = '%s_mutated_genes_aa.fasta' % SNV_matrix_cdc_basename


# store sequences in dict
ref_seq_dict = {}
for each_seq in SeqIO.parse(combined_ref_fasta, 'fasta'):
    ref_seq_dict[each_seq.id] = str(each_seq.seq)


# get seq id list
seq_id_list = []
for cds in open(combined_ref_gff):
    if (cds.startswith('2.10')) or (cds.startswith('D2')):
        seq_ID = cds.strip().split('\t')[0]
        if seq_ID not in seq_id_list:
            seq_id_list.append(seq_ID)


# get ending positions for all genes
ORF_ending_pos_dict = {}
ORF_seq_id_dict = {}
ORF_strand_dict = {}
for cds2 in open(combined_ref_gff):
    if (cds2.startswith('2.10')) or (cds2.startswith('D2')):
        cds2_split = cds2.strip().split('\t')
        seq_id2 = cds2_split[0]
        gene_id2 = cds2_split[8].split(';')[0].split('=')[1]
        start_pos2 = int(cds2_split[3])
        end_pos2 = int(cds2_split[4])
        strand = cds2_split[6]
        ORF_ending_pos_dict[gene_id2] = [start_pos2, end_pos2]
        ORF_seq_id_dict[gene_id2] = seq_id2
        ORF_strand_dict[gene_id2] = strand


# initialize coding region dict
coding_region_dict = {}
for each_id in seq_id_list:
    coding_region_dict[each_id] = []

# add coding pos to coding_region_dict
for each_cds in open(combined_ref_gff):
    if (each_cds.startswith('2.10')) or (each_cds.startswith('D2')):
        each_cds_split = each_cds.strip().split('\t')
        seq_id = each_cds_split[0]
        start_pos = int(each_cds_split[3])
        end_pos = int(each_cds_split[4])
        pos_list = list(range(start_pos, end_pos + 1))
        for each_pos in pos_list:
            coding_region_dict[seq_id].append(each_pos)


output_handle = open(output_mutated_genes, 'w')
for each_snv in open(SNV_matrix_cdc):
    if not each_snv.startswith('\t'):
#        print(each_snv)
        each_snv_seq = each_snv.strip().split(',')[0].split('|')[0]
        each_snv_pos = each_snv.strip().split(',')[0].split('|')[1]
        each_snv_pos_wt = each_snv.strip().split(',')[0].split('|')[2]
        each_snv_pos_v = each_snv.strip().split(',')[0].split('|')[3]

        # get all affected genes for
        each_snv_pos = int(each_snv_pos)

        # get SNV location
        location = ''
        if each_snv_pos in coding_region_dict[each_snv_seq]:
            location = 'Coding'
        else:
            location = 'Intergenic'
            #print('%s\t%s' % (each_snv.strip().split('\t')[0], location))
            output_handle.write('%s\t%s\tNA\tNA\tNA\tNA\n' % (each_snv.strip().split('\t')[0], location))

        # get mutation type
        if location == 'Coding':
            for each_gene in ORF_ending_pos_dict:
                if (each_snv_seq == ORF_seq_id_dict[each_gene]) and (ORF_ending_pos_dict[each_gene][0] <= each_snv_pos <= ORF_ending_pos_dict[each_gene][1]):

                    start_pos_raw = ORF_ending_pos_dict[each_gene][0]
                    end_pos_raw = ORF_ending_pos_dict[each_gene][1]
                    snv_pos_raw = each_snv_pos
                    start_pos_rescaled = ORF_ending_pos_dict[each_gene][0] - (ORF_ending_pos_dict[each_gene][0] - 1)
                    end_pos_rescaled = ORF_ending_pos_dict[each_gene][1] - (ORF_ending_pos_dict[each_gene][0] - 1)
                    snv_pos_rescaled = each_snv_pos - (ORF_ending_pos_dict[each_gene][0] - 1)
                    mutation_type = ''

                    if each_snv_pos_v == '-':
                        mutation_type = 'Frameshift'
                        output_handle.write('%s\t%s\t%s\t%s\tNA\t%s\n' % (each_snv.strip().split('\t')[0], location, ORF_strand_dict[each_gene], each_gene, mutation_type))

                    elif each_snv_pos_v in ['A', 'T', 'C', 'G']:

                        # get all reading frame
                        rf_pos_list = get_rf_pos_list(snv_pos_rescaled)

                        # get the mutated reading frame with rescaled position
                        mutated_rf_rescaled = []
                        for each_rf in rf_pos_list:
                            if snv_pos_rescaled in each_rf:
                                mutated_rf_rescaled = each_rf

                        # get the mutated reading frame with raw position
                        mutated_rf_raw = []
                        for each_rescaled_pos in mutated_rf_rescaled:
                            raw_pos = each_rescaled_pos + ORF_ending_pos_dict[each_gene][0] - 1
                            mutated_rf_raw.append(raw_pos)

                        # get sequence of raw and mutated reading frame
                        rf_seq_raw = ref_seq_dict[each_snv_seq][(mutated_rf_raw[0] - 1):(mutated_rf_raw[0] + 2)]
                        rf_seq_mutated = ''
                        for each_bp in mutated_rf_raw:
                            current_bp = ''
                            if each_bp == each_snv_pos:
                                current_bp = each_snv_pos_v
                            else:
                                current_bp = ref_seq_dict[each_snv_seq][each_bp - 1]
                            rf_seq_mutated += current_bp

                        if ORF_strand_dict[each_gene] == '-':
                            rf_seq_raw = str(Seq(rf_seq_raw, generic_dna).reverse_complement())
                            rf_seq_mutated = str(Seq(rf_seq_mutated, generic_dna).reverse_complement())

                        rf_seq_raw_aa = str(SeqRecord(Seq(rf_seq_raw)).seq.translate(table=transl_table))
                        rf_seq_mutated_aa = str(SeqRecord(Seq(rf_seq_mutated)).seq.translate(table=transl_table))

                        mutation_type_term = ''
                        if rf_seq_raw_aa == rf_seq_mutated_aa:
                            mutation_type_term = 'Silent'
                        elif rf_seq_mutated_aa == '*':
                            mutation_type_term = 'Nonsense'
                        else:
                            mutation_type_term = 'Missense'
                        aa_mutation = '%s(%s)->%s(%s)' % (rf_seq_raw, rf_seq_raw_aa, rf_seq_mutated, rf_seq_mutated_aa)

                        # print out
                        for_write = '%s\t%s\t%s\t%s\t%s\t%s\n' % (each_snv.strip().split('\t')[0], location, ORF_strand_dict[each_gene], each_gene, aa_mutation, mutation_type_term)
                        output_handle.write(for_write)
output_handle.close()


# get_mutation_cate_summary
mutation_cate_dict = get_mutation_cate_summary(output_mutated_genes, output_mutated_genes_cate)


# report
print('The number of affected genes: %s' % len(mutation_cate_dict))


################################################### export sequences ###################################################

# get the sequence of affected genes
output_seq_nc_handle = open(output_seq_nc, 'w')
for each_gene_nc in SeqIO.parse(combined_ref_ffn, 'fasta'):
    if each_gene_nc.id in mutation_cate_dict:
        SeqIO.write(each_gene_nc, output_seq_nc_handle, 'fasta')
output_seq_nc_handle.close()

output_seq_aa_handle = open(output_seq_aa, 'w')
for each_gene_aa in SeqIO.parse(combined_ref_faa, 'fasta'):
    if each_gene_aa.id in mutation_cate_dict:
        SeqIO.write(each_gene_aa, output_seq_aa_handle, 'fasta')
output_seq_aa_handle.close()


########################################## read function annotation into dict ##########################################

# store annotation results in dicts
gene_COG_cat_dict = {}
gene_COG_id_dict = {}
gene_COG_function_dict = {}
for each_snv3 in open(annotation_file):
    each_snv3_split = each_snv3.strip().split('\t')
    gene_id = each_snv3_split[0]
    COG_cat = each_snv3_split[10]
    COG_id = each_snv3_split[8]
    COG_function = each_snv3_split[9]
    gene_COG_cat_dict[gene_id] = COG_cat
    gene_COG_id_dict[gene_id] = COG_id
    gene_COG_function_dict[gene_id] = COG_function


############################################### combine mutation effect ################################################


output_summary_handle = open(output_summary, 'w')
for each_snv4 in open(output_mutated_genes):
    snv_id4 = each_snv4.strip().split('\t')[0]
    gene_id4 = each_snv4.strip().split('\t')[3]

    current_COG_id2 = ''
    current_COG_function2 = ''
    if gene_id4 in gene_COG_id_dict:
        current_COG_id2 = gene_COG_id_dict[gene_id4]
        current_COG_function2 = gene_COG_function_dict[gene_id4]
    else:
        current_COG_id2 = 'NA'
        current_COG_function2 = 'NA'

    for_write2 = '%s\t%s\t%s\n' % (each_snv4.strip(), current_COG_id2, current_COG_function2)
    output_summary_handle.write(for_write2)
output_summary_handle.close()


################################################### remove tmp file ####################################################

os.system('rm %s' % output_mutated_genes_cate)
os.system('rm %s' % output_mutated_genes)
os.system('rm %s' % output_seq_nc)
os.system('rm %s' % output_seq_aa)

