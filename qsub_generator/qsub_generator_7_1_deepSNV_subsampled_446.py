import os
import glob
import shutil

###################################### CONFIGURATION ######################################

nodes_number = 1
ppn_number = 1
memory = 10
walltime_needed = '11:59:00'
email = 'wythe1987@163.com'
modules_needed = ['python/3.5.2', 'R/3.4.2']

wd = '/Users/songweizhi/Desktop/'
outputs_folder = 'qsub_deepSNV_subsampled_446'
wd_on_katana = '/srv/scratch/z5039045/Flow_cell_biofilm/4_1_deepSNV_subsampled_446'
sample_prefix_file = '/Users/songweizhi/Dropbox/Research/Flow_cell_datasets/sample_group.txt'

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

sample_id_to_group_dict = {}
for each in open(sample_prefix_file):
    each_split = each.strip().split('\t')
    sample_id_to_group_dict[each_split[0]] = each_split[1]


for sample in open(sample_prefix_file):
    sample_id = sample.strip().split('\t')[0]
    sample_group = sample_id_to_group_dict[sample_id]

    deepSNV_cmd = ''
    if sample_group == 'mono_210':
        deepSNV_cmd = 'python3 /srv/scratch/z5039045/Scripts/deep_SNV_runner.py -r ../0_References/2.10wt_illumina.fasta -q ../3_novoalign_subsampled_446/%s.bam -c ../3_novoalign_subsampled_446/210WTD0.bam' % (sample_id)

    if sample_group == 'mono_D2':
        deepSNV_cmd = 'python3 /srv/scratch/z5039045/Scripts/deep_SNV_runner.py -r ../0_References/D2_pacbio.fasta -q ../3_novoalign_subsampled_446/%s.bam -c ../3_novoalign_subsampled_446/D2D0.bam' % (sample_id)

    if sample_group == 'coculture':
        deepSNV_cmd = 'python3 /srv/scratch/z5039045/Scripts/deep_SNV_runner.py -r ../0_References/combined_references.fasta -q ../3_novoalign_subsampled_446/%s.bam -c ../3_novoalign_subsampled_446/coculture_D0.bam' % (sample_id)

    # print(deepSNV_cmd)

    output_handle = open('%s/qsub_deepSNV_subsampled_%s.sh' % (outputs_folder, sample_id), 'w')
    output_handle.write(header)
    output_handle.write(module_lines)
    output_handle.write('cd %s\n' % wd_on_katana)
    output_handle.write('%s\n' % deepSNV_cmd)
    output_handle.close()
