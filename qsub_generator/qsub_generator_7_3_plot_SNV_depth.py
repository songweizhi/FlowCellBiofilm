import os
import glob
import shutil

###################################### CONFIGURATION ######################################

nodes_number = 1
ppn_number = 1
memory = 10
walltime_needed = '11:59:00'
email = 'wythe1987@163.com'
modules_needed = ['python/3.5.2']

wd = '/Users/songweizhi/Desktop'
outputs_folder = 'qsub_get_SNV_depth_plot'
wd_on_katana = '/srv/scratch/z5039045/Flow_cell_biofilm/4_3_SNV_depth'

plot_sam_depth_script = '/srv/scratch/z5039045/Scripts/plot_sam_depth.py'
snv_cdc_file = 'deepSNV_output_summary_all_frequency_cdc.txt'
flanking_length = 50000
kmer = 1000


summary_file_header_all = ['1D9', '1D18', '1D27', '1D42', '5D9', '5D18', '5D27', '5D42', '9D9', '9D18', '9D27', '9D42', '2D9', '2D18', '2D27', '2D42', '6D9', '6D18', '6D27', '6D42', '10D9', '10D18', '10D27', '10D42', '4D9', '4D18', '4D27', '4D42', '8D9', '8D18', '8D27', '8D42', '12D9', '12D18', '12D27', '12D42']
reference_dict = {'1': '2.10wt_illumina.fasta',
                  '5': '2.10wt_illumina.fasta',
                  '9': '2.10wt_illumina.fasta',
                  '2': 'D2_pacbio.fasta',
                  '6': 'D2_pacbio.fasta',
                  '10': 'D2_pacbio.fasta',
                  '4': 'combined_references.fasta',
                  '8': 'combined_references.fasta',
                  '12': 'combined_references.fasta'}

ref_length_dict = {'2.10_chromosome': 3758219,
                   '2.10_plasmid1': 237747,
                   '2.10_plasmid2': 94490,
                   '2.10_plasmid3': 70353,
                   'D2_c': 4010148,
                   'D2_p': 1010395}

###########################################################################################

os.chdir(wd)

# create outputs folder
if not os.path.exists(outputs_folder):
    os.makedirs(outputs_folder)
else:
    shutil.rmtree(outputs_folder)
    os.makedirs(outputs_folder)

# Prepare header
line_1 = '#!/bin/bash\n'
line_2 = '#PBS -l nodes=' + str(nodes_number) + ':ppn=' + str(ppn_number) + '\n'
line_3 = '#PBS -l vmem=' + str(memory) + 'gb\n'
line_4 = '#PBS -l walltime=' + walltime_needed + '\n'
line_5 = '#PBS -j oe\n'
line_6 = '#PBS -M ' + email + '\n'
line_7 = '#PBS -m ae\n'
header = line_1 + line_2 + line_3 + line_4 + line_5 + line_6 + line_7

# Prepare module lines
module_lines = ''
for module in modules_needed:
    module_lines += 'module load ' + module + '\n'

################################################################

for each_snv in open(snv_cdc_file):
    if not each_snv.startswith('	'):
        each_snv_id = each_snv.strip().split('\t')[0]
        each_snv_id_split = each_snv_id.split('|')
        each_snv_pos = each_snv_id_split[1]
        ref_seq_id = each_snv_id_split[0]
        time_points = each_snv.strip().split('\t')[1:]

        n = 0
        for each_time_point in time_points:

            if each_time_point != '0':
                current_time_point = summary_file_header_all[n]
                depth_file = '/srv/scratch/z5039045/Flow_cell_biofilm/3_novoalign/%s.depth' % current_time_point
                reference_file = '/srv/scratch/z5039045/Flow_cell_biofilm/0_References/%s' % reference_dict[current_time_point.split('D')[0]]
                reference_seq = ref_seq_id
                each_snv_id_modified = '_'.join(each_snv_id.split('|'))
                plot_name = '%s_%s_f%sbp_%smer' % (each_snv_id_modified, current_time_point, flanking_length, kmer)
                cmd = ''

                if '-' in each_snv_pos:
                    each_snv_pos_start = int(each_snv_pos.split('-')[0])
                    each_snv_pos_end = int(each_snv_pos.split('-')[1])

                    plot_start = each_snv_pos_start - flanking_length
                    if plot_start <= 0:
                        plot_start = 1

                    plot_end = each_snv_pos_end + flanking_length
                    if plot_end >= ref_length_dict[reference_seq]:
                        plot_end = ref_length_dict[reference_seq]

                    marker_lines = '%s,%s' % (each_snv_pos_start, each_snv_pos_end)
                    cmd = 'python3 %s -r %s -d %s -i %s -s %s -e %s -k %s -l %s -o %s' % (plot_sam_depth_script, reference_file, depth_file, reference_seq, plot_start, plot_end, kmer, marker_lines, plot_name)
                else:

                    plot_start = int(each_snv_pos) - flanking_length
                    if plot_start <= 0:
                        plot_start = 1

                    plot_end = int(each_snv_pos) + flanking_length
                    if plot_end >= ref_length_dict[reference_seq]:
                        plot_end = ref_length_dict[reference_seq]

                    marker_lines = each_snv_pos
                    cmd = 'python3 %s -r %s -d %s -i %s -s %s -e %s -k %s -l %s -o %s' % (plot_sam_depth_script, reference_file, depth_file, reference_seq, plot_start, plot_end, kmer, marker_lines, plot_name)

                output_handle = open('%s/qsub_get_SNV_reads_num_%s.sh' % (outputs_folder, plot_name), 'w')
                output_handle.write(header)
                output_handle.write(module_lines)
                output_handle.write('cd %s\n' % wd_on_katana)
                output_handle.write('%s\n' % cmd)
                output_handle.close()

            n += 1
