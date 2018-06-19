import os
import glob
import shutil

###################################### CONFIGURATION ######################################

nodes_number = 1
ppn_number = 1
memory = 60
walltime_needed = '05:59:00'
email = 'wythe1987@163.com'
modules_needed = []

wd = '/Users/songweizhi/Desktop/phylosift'
outputs_folder = 'qsub_phylosift_human_gut'
#outputs_folder = 'qsub_phylosift_NorthSea'
wd_on_katana = '/srv/scratch/z5039045/MetaCHIP/phylosift/phylosift_out_human_gut'
#wd_on_katana = '/srv/scratch/z5039045/MetaCHIP/phylosift/phylosift_out_NorthSea'

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

for each in open('human_gut_bins.txt'):
#for each in open('northsea_bins.txt'):
    each = each.strip()
    output_handle = open('%s/qsub_phylosift_%s.sh' % (outputs_folder, each), 'w')
    output_handle.write(header)
    output_handle.write(module_lines)
    output_handle.write('cd %s\n' % wd_on_katana)
    output_handle.write('/home/z5039045/phylosift_v1.0.1/bin/phylosift all /srv/scratch/z5039045/MetaCHIP/phylosift/human_gut_138bins/%s.fasta --out=/srv/scratch/z5039045/MetaCHIP/phylosift/phylosift_out_human_gut/%s\n' % (each, each))
    #output_handle.write('/home/z5039045/phylosift_v1.0.1/bin/phylosift all /srv/scratch/z5039045/MetaCHIP/phylosift/NorthSea_47bins/%s.fasta --out=/srv/scratch/z5039045/MetaCHIP/phylosift/phylosift_out_NorthSea/%s\n' % (each, each))

    output_handle.close()
