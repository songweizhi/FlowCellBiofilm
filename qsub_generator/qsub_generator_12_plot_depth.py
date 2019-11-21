import os
import glob
import shutil

###################################### CONFIGURATION ######################################

nodes_number = 1
ppn_number = 1
memory = 30
walltime_needed = '11:59:00'
email = 'wythe1987@163.com'
modules_needed = ['python/3.5.2']

wd = '/Users/songweizhi/Desktop'
outputs_folder = 'qsub_plot_sam_depth'
wd_on_katana = '/srv/scratch/z5039045/Flow_cell_biofilm/3_novoalign_subsampled_662'

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
sample_prefix_file = '/Users/songweizhi/Dropbox/Research/Flow_cell_datasets/sample_prefix_combined.txt'

for sample in open(sample_prefix_file):
    sample = sample.strip()

    ref_file_name = ''
    if sample.split('D')[0] in ['1', '5', '9']:
        ref_file_name = '2.10wt_illumina.fasta'

    elif sample.split('D')[0] in ['2', '6', '10']:
        ref_file_name = 'D2_pacbio.fasta'

    elif sample.split('D')[0] in ['4', '8', '12']:
        ref_file_name = 'combined_references.fasta'

    cmd = 'python3 /srv/scratch/z5039045/Scripts/plot_sam_depth.py -r ../0_References/%s -d %s.depth -k 1000\n' % (ref_file_name, sample)


    output_handle = open('%s/qsub_plot_depth_%s.sh' % (outputs_folder, sample), 'w')
    output_handle.write('%s\n' % header)
    output_handle.write(module_lines)
    output_handle.write('\ncd %s\n' % wd_on_katana)
    output_handle.write(cmd)
    output_handle.close()

