import os
import glob
import shutil

###################################### CONFIGURATION ######################################

nodes_number = 1
ppn_number = 1
memory = 10
walltime_needed = '02:59:00'
email = 'wythe1987@163.com'
modules_needed = []

wd = '/Users/songweizhi/Desktop/Flow_cell_biofilm'
outputs_folder = 'qsub_combine_reads'
wd_on_katana = '/srv/scratch/z5039045/Flow_cell_biofilm/raw_reads'

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


#print(header)
#print(module_lines)

################################################################

for each in open('flow_cell_biofilm_sample_prefix.txt'):
    each = each.strip()
    output_handle = open('%s/qsub_combine_reads_%s.sh' % (outputs_folder, each), 'w')
    run_1_reads_R1 = 'run_1_Q30/%s_R1_001_Q30_P.fastq' % each
    run_1_reads_R2 = 'run_1_Q30/%s_R2_001_Q30_P.fastq' % each
    run_2_reads_R1 = 'run_2_Q30/%s_R1_001_Q30_P.fastq' % each
    run_2_reads_R2 = 'run_2_Q30/%s_R2_001_Q30_P.fastq' % each
    run_3_reads_R1 = 'run_3_Q30/%s_R1_001_Q30_P.fastq' % each
    run_3_reads_R2 = 'run_3_Q30/%s_R2_001_Q30_P.fastq' % each
    output_handle.write(header)
    output_handle.write(module_lines)
    output_handle.write('cd %s\n' % wd_on_katana)
    output_handle.write('cat %s %s %s > combined_reads/%s_R1_Q30_P.fastq\n' % (run_1_reads_R1, run_2_reads_R1, run_3_reads_R1, each.split('_')[0]))
    output_handle.write('cat %s %s %s > combined_reads/%s_R2_Q30_P.fastq\n' % (run_1_reads_R2, run_2_reads_R2, run_3_reads_R2, each.split('_')[0]))
    output_handle.close()
