import os
import glob
import shutil

###################################### CONFIGURATION ######################################

nodes_number = 1
ppn_number = 6
memory = 60
walltime_needed = '23:59:00'
email = 'wythe1987@163.com'
modules_needed = ['novocraft/3.04.06']

wd = '/Users/songweizhi/Dropbox/Research/Flow_cell'
outputs_folder = 'qsub_novoalign_mono_210'
#outputs_folder = 'qsub_novoalign_mono_D2'
#outputs_folder = 'qsub_novoalign_coculture'

wd_on_katana = '/srv/scratch/z5039045/Flow_cell_biofilm/3_novoalign'

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
sample_prefix_file = '/Users/songweizhi/Dropbox/Research/Flow_cell/sample_prefix_mono_210.txt'
#sample_prefix_file = '/Users/songweizhi/Dropbox/Research/Flow_cell/sample_prefix_mono_D2.txt'
#sample_prefix_file = '/Users/songweizhi/Dropbox/Research/Flow_cell/sample_prefix_coculture.txt'

for each in open(sample_prefix_file):
    each = each.strip()
    output_handle = open('%s/qsub_gunzip_%s.sh' % (outputs_folder, each), 'w')
    reads_R1 = '/srv/scratch/z5039045/Flow_cell_biofilm/2_combined_reads/%s_R1_Q30_P.fastq' % each
    reads_R2 = '/srv/scratch/z5039045/Flow_cell_biofilm/2_combined_reads/%s_R2_Q30_P.fastq' % each
    output_handle.write(header)
    output_handle.write(module_lines)
    output_handle.write('cd %s\n' % wd_on_katana)
    output_handle.write('novoalign -d /srv/scratch/z5039045/Flow_cell_biofilm/0_References/2.10wt_illumina.fasta.index -f %s %s -g 35 -x 10 -oSAM -oFullNW > %s.sam\n' % (reads_R1, reads_R2, each))
    print('novoalign -d /srv/scratch/z5039045/Flow_cell_biofilm/0_References/2.10wt_illumina.fasta.index -f %s %s -g 35 -x 10 -oSAM -oFullNW > %s.sam' % (reads_R1, reads_R2, each))

    #output_handle.write('novoalign -d /srv/scratch/z5039045/Flow_cell_biofilm/0_References/D2_pacbio.fasta.index -f %s %s -g 35 -x 10 -oSAM -oFullNW > %s.sam\n' % (reads_R1, reads_R2, each))
    #print('novoalign -d /srv/scratch/z5039045/Flow_cell_biofilm/0_References/D2_pacbio.fasta.index -f %s %s -g 35 -x 10 -oSAM -oFullNW > %s.sam' % (reads_R1, reads_R2, each))

    #output_handle.write('novoalign -d /srv/scratch/z5039045/Flow_cell_biofilm/0_References/combined_references.fasta.index -f %s %s -g 35 -x 10 -oSAM -oFullNW > %s.sam\n' % (reads_R1, reads_R2, each))
    #print('novoalign -d /srv/scratch/z5039045/Flow_cell_biofilm/0_References/combined_references.fasta.index -f %s %s -g 35 -x 10 -oSAM -oFullNW > %s.sam' % (reads_R1, reads_R2, each))

    output_handle.close()
