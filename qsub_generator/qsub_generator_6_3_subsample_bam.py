import os
import glob
import shutil

###################################### calculate average depth from bam depth file #####################################

total_length_dict = {'1': 4090456, '5': 4090456, '9': 4090456, '2': 5020543, '6': 5020543, '10': 5020543, '4': 9110999, '8': 9110999, '12': 9110999, '210WT': 4090456, 'D2': 5020543, 'coculture_': 9110999}
depth_file_folder = '/Users/songweizhi/Desktop/depth_file'
depth_file_re = '%s/*.depth' % depth_file_folder
depth_file_list = [os.path.basename(file_name) for file_name in glob.glob(depth_file_re)]

# for each_file in depth_file_list:
#
#     if each_file == 'D2D0.depth':
#         depth_file_treatment_id = 'D2'
#     else:
#         depth_file_treatment_id = each_file.split('D')[0]
#
#     pwd_depth_file = '%s/%s' % (depth_file_folder, each_file)
#
#     total_depth = 0
#     for each_bp in open(pwd_depth_file):
#         each_bp_split = each_bp.strip().split('\t')
#         depth_value = int(each_bp_split[2])
#         total_depth += depth_value
#
#     average_depth = total_depth/(total_length_dict[depth_file_treatment_id])
#     average_depth = float("{0:.0f}".format(average_depth))
#     print('%s\t%s' % (each_file, average_depth))


############################################### get subsample percentage ###############################################

depth_summary_file = '/Users/songweizhi/Dropbox/Research/Flow_cell/subsample_wd/bam_depth.txt'

depth_list = []
depth_dict = {}
for each_depth in open(depth_summary_file):
    each_depth_split = each_depth.strip().split('\t')
    depth_list.append(int(each_depth_split[1]))
    depth_dict[each_depth_split[0]] = int(each_depth_split[1])

min_depth = min(depth_list)

subsample_percent_file = '/Users/songweizhi/Dropbox/Research/Flow_cell/subsample_wd/bam_depth_subsample_percent.txt'
subsample_percent_file_handle = open(subsample_percent_file, 'w')
subsample_percent_file_handle.write('Sample\tDepth\tPercent\tFinal_depth\n')
subsample_percent_dict = {}
for each_key in depth_dict:
    percent = min_depth/depth_dict[each_key]
    percent = float("{0:.2f}".format(percent))
    subsample_percent_file_handle.write('%s\t%s\t%s\t%s\n' % (each_key, depth_dict[each_key], percent, min_depth))
    subsample_percent_dict[each_key] = percent
subsample_percent_file_handle.close()


################################################### prepare qsub file ##################################################

############ CONFIGURATION ############

nodes_number = 1
ppn_number = 1
memory = 10
walltime_needed = '11:59:00'
email = 'wythe1987@163.com'
modules_needed = ['samtools/1.7']

wd = '/Users/songweizhi/Desktop'
outputs_folder = 'qsub_subsample'
wd_on_katana = '/srv/scratch/z5039045/Flow_cell_biofilm/3_novoalign'

#######################################

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

#######################################

for each_percent in subsample_percent_dict:

    cmd = 'samtools view -s %s -b %s.bam > subsampled/%s.bam' % (subsample_percent_dict[each_percent], each_percent, each_percent)

    output_handle = open('%s/qsub_subsample_%s.sh' % (outputs_folder, each_percent), 'w')
    output_handle.write(header)
    output_handle.write(module_lines)
    output_handle.write('cd %s\n' % wd_on_katana)
    output_handle.write('%s\n' % cmd)
    output_handle.write('cd subsampled\n')
    output_handle.write('samtools index %s.bam\n' % each_percent)
    output_handle.close()
