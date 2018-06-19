import os
import glob
import shutil

###################################### CONFIGURATION ######################################

nodes_number = 1
ppn_number = 6
memory = 30
walltime_needed = '11:59:00'
email = 'wythe1987@163.com'
modules_needed = ['bowtie/2.3.2']

wd = '/Users/songweizhi/Desktop'
outputs_folder = 'qsub_10_mapHGT_mapping'
wd_on_katana = '/srv/scratch/z5039045/Flow_cell_biofilm/6_mapHGT'

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
sample_prefix_file = 'qsub_generator_8_uniq_ids.txt'

for sample in open(sample_prefix_file):
    sample = sample.strip()
    r1 = '%s_R1_Q30_P.fasta' % sample
    r2 = '%s_R2_Q30_P.fasta' % sample
    pwd_r1 = '/srv/scratch/z5039045/Flow_cell_biofilm/5_Daisy/%s' % r1
    pwd_r2 = '/srv/scratch/z5039045/Flow_cell_biofilm/5_Daisy/%s' % r2
    cmd = 'bowtie2 -x combined_ref -1 %s -2 %s -S %s.sam -p 6 -f\n' % (pwd_r1, pwd_r2, sample)

    output_handle = open('%s/qsub_mapping_%s.sh' % (outputs_folder, sample), 'w')
    output_handle.write('%s\n' % header)
    output_handle.write(module_lines)
    output_handle.write('\ncd %s\n' % wd_on_katana)
    output_handle.write(cmd)
    output_handle.close()

