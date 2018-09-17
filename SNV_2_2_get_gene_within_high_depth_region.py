import os


############################################### input file and parameters ##############################################

# wd
wd = '/Users/songweizhi/Dropbox/Research/Flow_cell_datasets/high_depth_region'

gff_file = '2.10.gff'
annotation_file = 'combined_ref_aa.faa.COG.arCOG.kegg'
# starting_bp = 866221
# ending_bp = 941248
sequence_id = '2.10_chromosome'
starting_bp = 2350500
ending_bp = 2355408

########################################################################################################################

os.chdir(wd)
output_file_name = '%s/%s_%s_%s_gene_annotation.txt' % ('/Users/songweizhi/Desktop', sequence_id, starting_bp, ending_bp)


gene_list = []
for each_cds in open(gff_file):
    if not each_cds.startswith('#'):
        each_cds_split = each_cds.strip().split('\t')
        seq_id = each_cds_split[0]
        start_bp = int(each_cds_split[3])
        end_bp = int(each_cds_split[4])
        gene_id = each_cds_split[8].split(';')[0].split('=')[-1]
        if (seq_id == sequence_id) and (start_bp >= starting_bp) and (end_bp <= ending_bp):
            gene_list.append(gene_id)

print(gene_list)
print(len(gene_list))

output_file_handle = open(output_file_name, 'w')
for each_snv3 in open(annotation_file):

    each_snv3_split = each_snv3.strip().split('\t')
    g_id = each_snv3_split[0]
    COG_id = each_snv3_split[8]
    COG_function = each_snv3_split[9]

    if g_id in gene_list:

        for_out = '%s\t%s\t%s\n' % (g_id, COG_id, COG_function)

        print(for_out)
        output_file_handle.write(for_out)
output_file_handle.close()




