import os
import shutil
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna


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


############################################### input file and parameters ##############################################

parser = argparse.ArgumentParser(description='', add_help=False)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')

optional.add_argument('-h', action='help', help='Show this help message and exit')
required.add_argument('-min_both', dest='MIN_BOTH', nargs='?', required=True, type=int, help='The minimum number of reads harboring SNV')
required.add_argument('-min_each', dest='MIN_EACH', nargs='?', required=True, type=int, help='The minimum number of reads harboring SNV at each direction')
required.add_argument('-strand_bias', dest='STRAND_BIAS', nargs='?', required=True, type=int, help='strand_bias cutoff')
required.add_argument('-depth_diff', dest='DEPTH_DIFF', nargs='?', required=True,  type=int, help='depth difference cutoff')
required.add_argument('-deplen', dest='DEPLEN', nargs='?', required=True, type=int, help='flanking length for mean depth calculation')
required.add_argument('-sep_plot', dest='SEP_PLOT', nargs='?', required=False, type=int, help='separate depth plot with provide cutoff')

args = vars(parser.parse_args())
min_var_reads_num = args['MIN_BOTH']
min_at_each_direction = args['MIN_EACH']
strand_bias_cutoff = args['STRAND_BIAS']
depth_difference_cutoff = args['DEPTH_DIFF']
mean_depth_len = args['DEPLEN']
separate_plot = args['SEP_PLOT']

combined_ref_fasta = '/Users/songweizhi/Dropbox/Research/Flow_cell_datasets/reference_files/combined_ref.fasta'
combined_ref_gff = '/Users/songweizhi/Dropbox/Research/Flow_cell_datasets/reference_files/combined_ref.gff'
combined_ref_ffn = '/Users/songweizhi/Dropbox/Research/Flow_cell_datasets/reference_files/combined_ref.ffn'
combined_ref_faa = '/Users/songweizhi/Dropbox/Research/Flow_cell_datasets/reference_files/combined_ref.faa'
annotation_file = '/Users/songweizhi/Dropbox/Research/Flow_cell_datasets/reference_files/combined_ref_aa.faa.COG.arCOG.kegg'
transl_table = 11
pwd_QC_txt = 'SNV_QC.txt'


######################################################### Main #########################################################

qualified_SNVs_even_flanking_depth_file = 'SNV_QC_even_depth.txt'
qualified_SNVs_diff_flanking_depth_file = 'SNV_QC_diff_depth.txt'

qualified_SNVs_even_flanking_depth_file_handle = open(qualified_SNVs_even_flanking_depth_file, 'w')
qualified_SNVs_diff_flanking_depth_file_handle = open(qualified_SNVs_diff_flanking_depth_file, 'w')

n_tst_total_unqualified = []
n_tst_each_unqualified = []
strand_bias_unqualified = []
qualified_SNVs = []
qualified_SNVs_even_flanking_depth = []
qualified_SNVs_diff_flanking_depth = []

plot_prefix_even = []
plot_prefix_diff = []
plot_prefix_unqualified = []

total_num = 0
for each_snv in open(pwd_QC_txt):

    if not each_snv.startswith('sample'):
        each_snv_split = each_snv.strip().split(',')
        each_snv_sample = each_snv_split[0]
        each_snv_chr = each_snv_split[1]
        each_snv_pos = int(each_snv_split[2])
        each_snv_ref = each_snv_split[3]
        each_snv_var = each_snv_split[4]
        each_snv_p_val = float(each_snv_split[5])
        each_snv_freq = float(each_snv_split[6])
        each_snv_n_tst_b = int(each_snv_split[7])
        each_snv_n_tst_fw = int(each_snv_split[8])
        each_snv_n_tst_bw = int(each_snv_split[9])
        each_snv_strand_bias = float(each_snv_split[10])
        each_snv_mean_depth_diff = float(each_snv_split[14])

        each_snv_key_no_sample_id = '%s|%s|%s|%s' % (each_snv_chr, each_snv_pos, each_snv_ref, each_snv_var)
        each_snv_key_with_sample_id = '%s|%s|%s|%s|%s' % (each_snv_sample, each_snv_chr, each_snv_pos, each_snv_ref, each_snv_var)
        plot_prefix = '%s_%s_%s_%s_%s' % (each_snv_chr, each_snv_pos, each_snv_ref, each_snv_var, each_snv_sample)

        # get SNVs with unqualified number of reads harbouring SNV
        if each_snv_n_tst_b < min_var_reads_num:
            n_tst_total_unqualified.append(each_snv_key_with_sample_id)

            if plot_prefix not in plot_prefix_unqualified:
                plot_prefix_unqualified.append(plot_prefix)

        # get SNVs with unqualified number of reads in each direction
        if (each_snv_n_tst_fw < min_at_each_direction) or (each_snv_n_tst_bw < min_at_each_direction):
            n_tst_each_unqualified.append(each_snv_key_with_sample_id)

            if plot_prefix not in plot_prefix_unqualified:
                plot_prefix_unqualified.append(plot_prefix)

        # get SNVs with unqualified strand bias
        if each_snv_strand_bias > strand_bias_cutoff:
            strand_bias_unqualified.append(each_snv_key_with_sample_id)

            if plot_prefix not in plot_prefix_unqualified:
                plot_prefix_unqualified.append(plot_prefix)

        # get SNVs with flanking depth difference higher than defined cutoff
        if (each_snv_n_tst_b >= min_var_reads_num) and (not ((each_snv_n_tst_fw < min_at_each_direction) or (each_snv_n_tst_bw < min_at_each_direction))) and (each_snv_strand_bias <= strand_bias_cutoff):
            qualified_SNVs.append(each_snv_key_with_sample_id)

        # get qualified SNV with similar flanking depth
        if (each_snv_n_tst_b >= min_var_reads_num) and (each_snv_n_tst_fw >= min_at_each_direction) and (each_snv_n_tst_bw >= min_at_each_direction) and (each_snv_strand_bias <= strand_bias_cutoff) and (each_snv_mean_depth_diff <= depth_difference_cutoff):
            qualified_SNVs_even_flanking_depth.append(each_snv_key_with_sample_id)
            plot_prefix_even.append(plot_prefix)
            qualified_SNVs_even_flanking_depth_file_handle.write(each_snv)

        # get qualified SNV with different flanking depth
        if (each_snv_n_tst_b >= min_var_reads_num) and (each_snv_n_tst_fw >= min_at_each_direction) and (each_snv_n_tst_bw >= min_at_each_direction) and (each_snv_strand_bias <= strand_bias_cutoff) and (each_snv_mean_depth_diff > depth_difference_cutoff):
            qualified_SNVs_diff_flanking_depth.append(each_snv_key_with_sample_id)
            plot_prefix_diff.append(plot_prefix)
            qualified_SNVs_diff_flanking_depth_file_handle.write(each_snv)

        total_num += 1

qualified_SNVs_even_flanking_depth_file_handle.close()
qualified_SNVs_diff_flanking_depth_file_handle.close()


# For report
print('\n############################## report ##############################')
print('The total number of detected SNVs: %s' % total_num)
print('The number of SNVs with less than %s reads harboring it: %s' % (min_var_reads_num, len(n_tst_total_unqualified)))
print('The number of SNVs with reads only from one direction: %s' % len(n_tst_each_unqualified))
print('The number of SNVs with strand bias higher than %s: %s' % (strand_bias_cutoff, len(strand_bias_unqualified)))
print('The number of unqualified SNVs: %s' % len(plot_prefix_unqualified))
print('The number of qualified SNVs: %s' % len(qualified_SNVs))
print('The number of qualified SNVs with flanking depth (%sbp) difference higher than %s: %s' % (mean_depth_len, depth_difference_cutoff, len(qualified_SNVs_diff_flanking_depth)))
print('The number of qualified SNVs with flanking depth (%sbp) difference not higher than %s: %s' % (mean_depth_len, depth_difference_cutoff, len(qualified_SNVs_even_flanking_depth)))
print('Details exported to %s and %s.' % (qualified_SNVs_even_flanking_depth_file, qualified_SNVs_diff_flanking_depth_file))


################################ separate plots according to provided difference cutoff ################################

# separate plot files
if separate_plot == 1:

    pwd_plot_folder = 'SNV_depth_plot'
    pwd_plot_folder_even = '%s_even_%s' % (pwd_plot_folder, depth_difference_cutoff)
    pwd_plot_folder_diff = '%s_diff_%s' % (pwd_plot_folder, depth_difference_cutoff)
    pwd_plot_folder_unqualified = '%s_unqualified_%s' % (pwd_plot_folder, depth_difference_cutoff)

    # prepare folder
    if os.path.isdir(pwd_plot_folder_even):
        shutil.rmtree(pwd_plot_folder_even)
        os.mkdir(pwd_plot_folder_even)
    else:
        os.mkdir(pwd_plot_folder_even)

    if os.path.isdir(pwd_plot_folder_diff):
        shutil.rmtree(pwd_plot_folder_diff)
        os.mkdir(pwd_plot_folder_diff)
    else:
        os.mkdir(pwd_plot_folder_diff)

    if os.path.isdir(pwd_plot_folder_unqualified):
        shutil.rmtree(pwd_plot_folder_unqualified)
        os.mkdir(pwd_plot_folder_unqualified)
    else:
        os.mkdir(pwd_plot_folder_unqualified)


    # get plots with similar flanking deoth
    for each_even in plot_prefix_even:
        pwd_plot = '%s/%s*' % (pwd_plot_folder, each_even)
        cmd = 'cp %s %s/' % (pwd_plot, pwd_plot_folder_even)
        os.system(cmd)

    # get plots with similar flanking deoth
    for each_diff in plot_prefix_diff:
        pwd_plot = '%s/%s*' % (pwd_plot_folder, each_diff)
        cmd = 'cp %s %s/' % (pwd_plot, pwd_plot_folder_diff)
        os.system(cmd)

    # get plots for unqualified SNVs
    for each_unqualified in plot_prefix_unqualified:
        pwd_plot = '%s/%s*' % (pwd_plot_folder, each_unqualified)
        cmd = 'cp %s %s/' % (pwd_plot, pwd_plot_folder_unqualified)
        os.system(cmd)


####################################### compare between Monoculture and Coculture ######################################

mono_210_SNV_list = []
co_210_SNV_list = []
mono_D2_SNV_list = []
co_D2_SNV_list = []
total_SNV = 0
for each_SNV in open(qualified_SNVs_even_flanking_depth_file):
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
print('\n############################## report ##############################')
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
#output_file_path = '/Users/songweizhi/Dropbox/Research/Flow_cell_datasets/Monoculture_VS_Coculture'

SNV_210_mono_uniq_file =   'SNV_210_mono_uniq.txt'
SNV_210_co_uniq_file =     'SNV_210_co_uniq.txt'
SNV_210_concurrence_file = 'SNV_210_concurrence.txt'
SNV_D2_mono_uniq_file =    'SNV_D2_mono_uniq.txt'
SNV_D2_co_uniq_file =      'SNV_D2_co_uniq.txt'
SNV_D2_concurrence_file =  'SNV_D2_concurrence.txt'


# write out
write_out(SNV_210_mono_uniq_file, SNV_210_mono_uniq)
write_out(SNV_210_co_uniq_file, SNV_210_co_uniq)
write_out(SNV_210_concurrence_file, SNV_210_concurrence)
write_out(SNV_D2_mono_uniq_file, SNV_D2_mono_uniq)
write_out(SNV_D2_co_uniq_file, SNV_D2_co_uniq)
write_out(SNV_D2_concurrence_file, SNV_D2_concurrence)


# store annotation results in dicts
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


################################################## get_mutated_genes ###################################################

SNV_matrix_cdc_list = [SNV_210_mono_uniq_file, SNV_210_co_uniq_file, SNV_210_concurrence_file, SNV_D2_mono_uniq_file, SNV_D2_co_uniq_file, SNV_D2_concurrence_file]
print(SNV_matrix_cdc_list)

print('\n############################## report ##############################')

for SNV_matrix_cdc in SNV_matrix_cdc_list:
    print(SNV_matrix_cdc)

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

    # get mutation_cate_summary
    mutation_cate_dict = get_mutation_cate_summary(output_mutated_genes, output_mutated_genes_cate)

    # report
    print('The number of affected genes for %s: %s' % (SNV_matrix_cdc, len(mutation_cate_dict)))


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

    # remove tmp file
    os.system('rm %s' % output_mutated_genes_cate)
    os.system('rm %s' % output_mutated_genes)
    os.system('rm %s' % output_seq_nc)
    os.system('rm %s' % output_seq_aa)


########################################################################################################################

print('\n############################## report ##############################')
strain_list = ['210', 'D2']
for strain in strain_list:

    # output file
    output_mutated_genes = 'mutated_genes_%s.txt' % strain
    output_mutated_genes_handle = open(output_mutated_genes, 'w')

    all_affected_gene_list = []
    gene_occurence_dict = {}
    gene_occurence_count_dict = {}

    for each_gene in open('SNV_%s_co_uniq_summary.txt' % strain):
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

    for each_gene in open('SNV_%s_mono_uniq_summary.txt' % strain):
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

    for each_gene in open('SNV_%s_concurrence_summary.txt' % strain):
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

    mutated_gene_nomo_uniq = 0
    mutated_gene_co_uniq = 0
    mutated_gene_concurrence = 0
    for each_mutated_gene in gene_occurence_dict:
        if gene_occurence_dict[each_mutated_gene] == ['mono']:
            mutated_gene_nomo_uniq += 1
        elif gene_occurence_dict[each_mutated_gene] == ['co']:
            mutated_gene_co_uniq += 1
        else:
            mutated_gene_concurrence += 1

    # for report:
    print('The number of mutated %s gene uniq to monoculture: %s' % (strain, mutated_gene_nomo_uniq))
    print('The number of mutated %s gene uniq to coculture: %s' % (strain, mutated_gene_co_uniq))
    print('The number of mutated %s gene with concurrence: %s' % (strain, mutated_gene_concurrence))
    print('Details exported to %s\n' % output_mutated_genes)

    # write out
    for each_affected_gene in all_affected_gene_list:

        existence = ''
        existence_profile = gene_occurence_dict[each_affected_gene]
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

        for_report = '%s\t%s\t%s\t%s\t%s\n' % (each_affected_gene, existence, gene_occurence_count_dict[each_affected_gene], current_COG_id2, current_COG_function2)
        output_mutated_genes_handle.write(for_report)

    output_mutated_genes_handle.close()

