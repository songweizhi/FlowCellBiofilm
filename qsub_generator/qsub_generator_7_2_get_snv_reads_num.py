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
#outputs_folder = 'qsub_get_SNV_reads_num'
outputs_folder = 'qsub_get_SNV_reads_num_subsampled'

#wd_on_katana = '/srv/scratch/z5039045/Flow_cell_biofilm/4_2_SNV_reads_num'
wd_on_katana = '/srv/scratch/z5039045/Flow_cell_biofilm/4_2_SNV_reads_num_subsampled'

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
sample_prefix_file = '/Users/songweizhi/Dropbox/Research/Flow_cell/sample_prefix_combined.txt'
#bam_file_folder = '/srv/scratch/z5039045/Flow_cell_biofilm/3_novoalign'
bam_file_folder = '/srv/scratch/z5039045/Flow_cell_biofilm/3_novoalign_subsampled'

#deepSNV_output_folder = '/srv/scratch/z5039045/Flow_cell_biofilm/4_1_deepSNV'
deepSNV_output_folder = '/srv/scratch/z5039045/Flow_cell_biofilm/4_1_deepSNV_subsampled'

pwd_python_script = '/srv/scratch/z5039045/Scripts/SNV_2_1_get_SNV_reads_num.py'

for sample in open(sample_prefix_file):
    sample_id = sample.strip()
    pwd_bam_file = '%s/%s.bam' % (bam_file_folder, sample_id)
    cmd = 'python3 %s -bf %s -df %s -i %s' % (pwd_python_script, bam_file_folder, deepSNV_output_folder, sample_id)
    print(cmd)
    output_handle = open('%s/qsub_get_SNV_reads_num_%s.sh' % (outputs_folder, sample_id), 'w')
    output_handle.write(header)
    output_handle.write(module_lines)
    output_handle.write('cd %s\n' % wd_on_katana)
    output_handle.write('%s\n' % cmd)
    output_handle.close()
