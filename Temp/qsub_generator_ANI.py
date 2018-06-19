import os
import glob
import shutil

###################################### CONFIGURATION ######################################

nodes_number = 1
ppn_number = 1
memory = 60
walltime_needed = '11:59:00'
email = 'wythe1987@163.com'
modules_needed = ['java/8u162']

wd = '/Users/songweizhi/Desktop'
outputs_folder = 'qsub_ANI'
wd_on_katana = '/srv/scratch/z5039045/ANI_calculator_wd/download_genome'

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
sample_prefix_file = 'ANI.txt'

n = 1
for each_line in open(sample_prefix_file):

    output_handle = open('%s/qsub_ANI_%s.sh' % (outputs_folder, n), 'w')
    output_handle.write('%s\n' % header)
    output_handle.write(module_lines)
    output_handle.write('\ncd %s\n' % '/srv/scratch/z5039045/ANI_calculator_wd/download_genome')
    output_handle.write(each_line)
    output_handle.close()
    n += 1

