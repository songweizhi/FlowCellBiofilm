import os
import glob
import shutil

###################################### CONFIGURATION ######################################

nodes_number = 1
ppn_number = 5
memory = 60
walltime_needed = '47:59:00'
email = 'wythe1987@163.com'
modules_needed = ['intel/17.0.4.196', 'gcc/6.2.0', 'clang/3.5.0', 'cmake/3.7.2', 'python/2.7.12', 'java/8u91', 'samtools/0.1.19', 'bedtools/2.27.1']

wd = '/Users/songweizhi/Desktop'
outputs_folder = 'qsub_Daisy'
wd_on_katana = '/srv/scratch/z5039045/Flow_cell_biofilm/5_Daisy'

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
    pwd_daisy = '/srv/scratch/z5039045/Softwares/daisy/daisy.py'
    for each in ['2.10_chromosome', '2.10_plasmid1', '2.10_plasmid2', '2.10_plasmid3', 'D2_c', 'D2_p']:

        ar = ''
        dr = ''
        d = ''
        if each.startswith('2.10_'):
            ar = '2.10wt_illumina.fasta'
            dr = 'D2_pacbio.fasta'
            d = 'D2_seq_ids.txt'
        if each.startswith('D2_'):
            ar = 'D2_pacbio.fasta'
            dr = '2.10wt_illumina.fasta'
            d = '210_seq_ids.txt'

        od = '%s_%s' % (sample, each)
        cmd_1 = 'mkdir %s' % od
        cmd_2 = 'python %s -r1 %s_R1_Q30_P.fasta -r2 %s_R2_Q30_P.fasta -ar %s -dr %s -a "%s" -d %s -new -od %s/ -t %s' % (pwd_daisy, sample, sample, ar, dr, each, d, od, od)

        output_handle = open('%s/qsub_daisy_%s.sh' % (outputs_folder, od), 'w')
        output_handle.write('%s\n' % header)
        output_handle.write('module unload intel/11.1.080\n')
        output_handle.write(module_lines)
        output_handle.write('cd ~\n')
        output_handle.write('. mypythonenv/bin/activate\n')
        output_handle.write('\ncd %s\n' % wd_on_katana)
        output_handle.write('%s\n' % cmd_1)
        output_handle.write('%s\n' % cmd_2)
        output_handle.close()
