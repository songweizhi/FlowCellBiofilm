import os
import glob
import shutil

###################################### CONFIGURATION ######################################

nodes_number = 1
ppn_number = 1
memory = 10
walltime_needed = '11:59:00'
email = 'wythe1987@163.com'
modules_needed = ['idba/1.1.3']

wd = '/Users/songweizhi/Desktop'
outputs_folder = 'qsub_fq2fa'
wd_on_katana = '/srv/scratch/z5039045/Flow_cell_biofilm/2_combined_reads'

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
sample_prefix_file = 'uniq_ids.txt'

for sample in open(sample_prefix_file):
    sample = sample.strip()
    r1 = '%s_R1_Q30_P.fastq' % sample.strip()
    r1_fa = '%s_R1_Q30_P.fasta' % sample.strip()
    r1_gz = '%s_R1_Q30_P.fastq.gz' % sample.strip()
    r2 = '%s_R2_Q30_P.fastq' % sample.strip()
    r2_fa = '%s_R2_Q30_P.fasta' % sample.strip()
    r2_gz = '%s_R2_Q30_P.fastq.gz' % sample.strip()

    cmd_1 = 'gunzip %s' % r1_gz
    cmd_2 = 'gunzip %s' % r2_gz

    cmd_3 = 'fq2fa %s %s' % (r1, r1_fa)
    cmd_4 = 'fq2fa %s %s' % (r2, r2_fa)

    cmd_5 = 'gzip %s' % r1
    cmd_6 = 'gzip %s' % r2

    cmd_7 = 'mv %s ../5_Daisy/' % r1_fa
    cmd_8 = 'mv %s ../5_Daisy/' % r2_fa

    output_handle = open('%s/qsub_fq2fa_%s.sh' % (outputs_folder, sample), 'w')
    output_handle.write(header)
    output_handle.write(module_lines)
    output_handle.write('cd %s\n' % wd_on_katana)
    output_handle.write('%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n' % (cmd_1, cmd_2, cmd_3, cmd_4, cmd_5, cmd_6, cmd_7, cmd_8))


    output_handle.close()
