import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna


def get_rf_pos_list(gene_len):

    rf_pos_list = []
    n = 0
    while n < (gene_len / 3):
        start_pos = 1 + 3 * n
        rf_pos_list.append([start_pos, start_pos + 1, start_pos + 2])
        n += 1

    return rf_pos_list


# specify wd and input files
wd = '/Users/songweizhi/Desktop/FC'

# output files
file_in = 'deepSNV_mutated_gene.txt'

# forward to working directory
os.chdir(wd)

gene_snv_dict = {}
gene_mutation_type_dict = {}

for each_snv in open(file_in):
    print(each_snv.strip())
    each_snv_split = each_snv.strip().split('\t')

    affected_gene = each_snv_split[3]

    if affected_gene != 'NA':

        affect = each_snv_split[5]

        print(affected_gene)

