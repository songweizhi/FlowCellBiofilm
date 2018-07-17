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

wd = '/Users/songweizhi/Desktop/666666'

strain = '210'
strain = 'D2'
flanking_length = 100000

outputs_folder = 'qsub_plot_diff_flanking_depth_SNVs_%s_%sbp' % (strain, flanking_length)
wd_on_katana = '/srv/scratch/z5039045/Flow_cell_biofilm/4_4_plot_diff_flanking_depth_SNVs_%sbp' % (flanking_length)

plot_sam_depth_script = '/srv/scratch/z5039045/Scripts/plot_sam_depth.py'
snv_cdc_file = 'deepSNV_output_qualified_diff_flanking_depth_summary_%s_frequency_cdc.txt' % strain

kmer = 1000

reference_dict = {'1': '2.10wt_illumina.fasta',
                  '5': '2.10wt_illumina.fasta',
                  '9': '2.10wt_illumina.fasta',
                  '2': 'D2_pacbio.fasta',
                  '6': 'D2_pacbio.fasta',
                  '10': 'D2_pacbio.fasta',
                  '4': 'combined_ref.fasta',
                  '8': 'combined_ref.fasta',
                  '12': 'combined_ref.fasta'}

os.chdir(wd)

###########################################################################################

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


current_matrix_header = []
for each_snv in open(snv_cdc_file):

    if each_snv.startswith('	'):
        current_matrix_header = each_snv.strip().split('\t')

    if not each_snv.startswith('	'):
        each_snv_id = each_snv.strip().split('\t')[0]
        each_snv_id_renamed = '_'.join(each_snv_id.split('|'))
        each_snv_id_split = each_snv_id.split('|')
        ref_seq_id = each_snv_id_split[0]
        each_snv_pos = each_snv_id_split[1]
        time_points = each_snv.strip().split('\t')[1:]

        qsub_file_name = 'qsub_plot_depth_%s.sh' % each_snv_id_renamed
        output_handle = open('%s/%s' % (outputs_folder, qsub_file_name), 'w')
        output_handle.write(header)
        output_handle.write(module_lines)
        output_handle.write('cd %s\n' % wd_on_katana)

        n = 0
        for each_time_point in time_points:
            
            if each_time_point != '0':
                current_time_point = current_matrix_header[n]
                depth_file = '%s.depth' % current_time_point
                reference_file = reference_dict[current_time_point.split('D')[0]]
                reference_seq = ref_seq_id
                plot_name = '%s_%s_f%sbp_%smer' % (each_snv_id_renamed, current_time_point, flanking_length, kmer)

                cmd = ''
                if '-' in each_snv_pos:
                    each_snv_pos_start = int(each_snv_pos.split('-')[0])
                    each_snv_pos_end = int(each_snv_pos.split('-')[1])
                    plot_start = each_snv_pos_start - flanking_length
                    plot_end = each_snv_pos_end + flanking_length
                    cmd = 'python3 %s -r ../0_References/%s -d ../3_novoalign/%s -i %s -s %s -e %s -k %s -l %s -o %s\n' % (plot_sam_depth_script, reference_file, depth_file, reference_seq, plot_start, plot_end, kmer, each_snv_pos, plot_name)

                else:
                    plot_start = int(each_snv_pos) - flanking_length
                    plot_end = int(each_snv_pos) + flanking_length
                    cmd = 'python3 %s -r ../0_References/%s -d ../3_novoalign/%s -i %s -s %s -e %s -k %s -l %s -o %s\n' % (plot_sam_depth_script, reference_file, depth_file, reference_seq, plot_start, plot_end, kmer, each_snv_pos, plot_name)

                output_handle.write(cmd)

            n += 1

        output_handle.close()

