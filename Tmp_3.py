import os
import shutil
import argparse
import numpy as np
from datetime import datetime
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from scipy import stats
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import os
import shutil


output_dir =                                '/Users/songweizhi/Desktop/gene_level_summary'
gene_level_output_tmp =                     '%s/SNV_QC_ncd_even_flk_depth_affected_genes_tmp.txt'           % output_dir
gene_level_output =                         '%s/SNV_QC_ncd_even_flk_depth_affected_genes.txt'           % output_dir
annotation_file_KEGG =                      '/Users/songweizhi/Research/Flow_cell_datasets/reference_files/combined_ref_ko_assignment_ABCD.txt'
annotation_file_COG_for_ko_unknown_genes =  '/Users/songweizhi/Research/Flow_cell_datasets/reference_files/function_unknown_genes_query_to_cog.txt'

########################################## read function annotation into dict ##########################################

# store annotation results in dicts
gene_KEGG_id_dict = {}
gene_KEGG_function_dict = {}
for each_snv3 in open(annotation_file_KEGG):
    if not each_snv3.startswith('Gene_id'):
        each_snv3_split = each_snv3.strip().split('\t')
        gene_id = each_snv3_split[0]
        if len(each_snv3_split) > 1:
            ko_id = each_snv3_split[4][2:]
            ko_function = each_snv3_split[8]
            gene_KEGG_id_dict[gene_id] = ko_id
            gene_KEGG_function_dict[gene_id] = ko_function
        else:
            gene_KEGG_id_dict[gene_id] = 'NA'
            gene_KEGG_function_dict[gene_id] = 'NA'


for each_snv3 in open(annotation_file_COG_for_ko_unknown_genes):
    if not each_snv3.startswith('Query'):
        each_snv3_split = each_snv3.strip().split('\t')
        gene_id = each_snv3_split[0]
        if len(each_snv3_split) > 1:
            cog_id = each_snv3_split[1]
            cog_function = each_snv3_split[3]
            gene_KEGG_id_dict[gene_id] = cog_id
            gene_KEGG_function_dict[gene_id] = cog_function


############################################## get summary at gene level ###############################################

# get file name list
gene_level_output_handle_tmp = open(gene_level_output_tmp, 'w')
for each_strain in ['210', 'D2']:
    affected_gene_uniq_list = set()
    affected_gene_to_SNV_dict = {}
    for each_parallel_cate in ['monoculture_uniq', 'coculture_uniq', 'concurrence']:

        # parse summary file
        summary_filename = '%s/SNV_QC_ncd_even_flk_depth_%s_%s_summary.txt' % (output_dir, each_strain, each_parallel_cate)
        for each_SNV in open(summary_filename):
            each_SNV_split = each_SNV.strip().split('\t')
            SNV_id = each_SNV_split[0]
            SNV_parallel_cate = each_SNV_split[1][-8:-1]
            #SNV_parallel_cate = each_SNV_split[1]

            affected_gene_id = each_SNV_split[4]
            affected_gene_ko = each_SNV_split[7]
            affected_gene_ko_desc = each_SNV_split[8]

            #print(each_SNV_split)
            if affected_gene_id != 'NA':
                gene_level_output_handle_tmp.write('%s\t%s\t%s\t%s\t%s\n' % (affected_gene_id, SNV_id, SNV_parallel_cate, affected_gene_ko, affected_gene_ko_desc))



            #print(SNV_parallel_cate)

            # if affected_gene_id != 'NA':
            #
            #     if affected_gene_id not in affected_gene_concurrence_dict:
            #         affected_gene_concurrence_dict[affected_gene_id] = [SNV_parallel_cate]
            #     else:
            #         affected_gene_concurrence_dict[affected_gene_id].append(SNV_parallel_cate)
            #
            #     # add to current strain affected_gene_uniq_list
            #     if affected_gene_id not in affected_gene_uniq_list:
            #         affected_gene_uniq_list.append(affected_gene_id)


    # print the report
    print(datetime.now().strftime('[%Y-%m-%d %H:%M:%S]') + ' The total number of %s affected gene (ncd_even_depth): %s' % (each_strain, len(affected_gene_uniq_list)))

gene_level_output_handle_tmp.close()

# sort gene_level_output
os.system('cat %s | sort > %s' % (gene_level_output_tmp, gene_level_output))
os.system('rm %s' % gene_level_output_tmp)

# for report
print(datetime.now().strftime('[%Y-%m-%d %H:%M:%S]') + ' Summaries at gene level exported to %s' % gene_level_output)
