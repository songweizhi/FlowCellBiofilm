import os
import glob
import shutil

###################################### CONFIGURATION ######################################

nodes_number = 1
ppn_number = 1
memory = 10
walltime_needed = '02:59:00'
email = 'wythe1987@163.com'
modules_needed = ['java/8u91', 'fastqc/0.11.6']

wd = '/Users/songweizhi/Desktop/Flow_cell_biofilm'
outputs_folder = 'qsub_QC_run_3'
wd_on_katana = '/srv/scratch/z5039045/Flow_cell_biofilm/raw_reads/run_3'

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

for each in open('flow_cell_biofilm_sample_prefix.txt'):
    each = each.strip()
    output_handle = open('%s/qsub_QC_%s.sh' % (outputs_folder, each), 'w')
    reads_R1 = '%s_R1_001.fastq.gz' % each
    reads_R2 = '%s_R2_001.fastq.gz' % each
    output_handle.write(header)
    output_handle.write(module_lines)
    output_handle.write('cd %s\n' % wd_on_katana)
    output_handle.write('java -jar /share/apps/trimmomatic/0.33/trimmomatic-0.33.jar PE %s_R1_001.fastq %s_R2_001.fastq %s_R1_001_Q30_P.fastq %s_R1_001_Q30_UP.fastq %s_R2_001_Q30_P.fastq %s_R2_001_Q30_UP.fastq ILLUMINACLIP:/share/apps/trimmomatic/0.33/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:30 TRAILING:30 CROP:240 HEADCROP:10 SLIDINGWINDOW:6:30 MINLEN:36\n' % (each, each, each, each, each, each))
    output_handle.write('fastqc %s_R1_001_Q30_P.fastq %s_R2_001_Q30_P.fastq\n' % (each, each))
    output_handle.write('rm %s_R1_001_Q30_UP.fastq\n' % each)
    output_handle.write('rm %s_R2_001_Q30_UP.fastq\n' % each)
    output_handle.write('rm %s_R1_001_Q30_P_fastqc.zip\n' % each)
    output_handle.write('rm %s_R2_001_Q30_P_fastqc.zip\n' % each)
    output_handle.close()

