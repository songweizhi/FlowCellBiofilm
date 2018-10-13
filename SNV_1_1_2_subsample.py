import os
import glob
import shutil


# define input
depth_summary = '/Users/songweizhi/Desktop/depth_summary.txt'
#depth_cutoff = 446.1
depth_cutoff = 662.4

nonsubsampled_bam_folder = '/srv/scratch/z5039045/Flow_cell_biofilm/3_novoalign_nonsubsampled'
#subsampled_bam_folder = '/srv/scratch/z5039045/Flow_cell_biofilm/3_novoalign_subsampled_446'
subsampled_bam_folder = '/srv/scratch/z5039045/Flow_cell_biofilm/3_novoalign_subsampled_662'



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
#wd_on_katana = '/srv/scratch/z5039045/Flow_cell_biofilm/3_novoalign_subsampled_446'
wd_on_katana = '/srv/scratch/z5039045/Flow_cell_biofilm/3_novoalign_subsampled_662'

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

# get subsample_ratio_dict
sample_id_list = []
subsample_ratio_dict = {}
for each_depth in open(depth_summary):

    each_depth_split = each_depth.strip().split('\t')
    each_depth_sample = each_depth_split[1][:-6]
    each_depth_ref = each_depth_split[2]
    each_depth_value = float(each_depth_split[3])
    each_depth_id = '%s_%s' % (each_depth_sample, each_depth_ref)

    # get sample_id_list
    if each_depth_sample not in sample_id_list:
        sample_id_list.append(each_depth_sample)

    # get subsample_ratio_dict
    subsample_ratio = 0
    if each_depth_value <= depth_cutoff:
        subsample_ratio = 1
    else:
        subsample_ratio = float("{0:.3f}".format(depth_cutoff/each_depth_value))

    subsample_ratio_dict[each_depth_id] = subsample_ratio



sample_list_mono_210 = ['1D18', '1D27', '1D42', '1D9', '5D18', '5D27', '5D42', '5D9', '9D18', '9D27', '9D42', '9D9', '210WTD0']
sample_list_mono_D2 = ['2D18', '2D27', '2D42', '2D9', '6D18', '6D27', '6D42', '6D9', '10D18', '10D27', '10D42', '10D9', 'D2D0']
sample_list_co = ['4D18', '4D27', '4D42', '4D9', '8D18', '8D27', '8D42', '8D9', '12D18', '12D27', '12D42', '12D9', 'coculture_D0']


for sample_mono_210 in sample_list_mono_210:

    # subsample_cmd
    subsample_cmd_mono_210 = 'samtools view -s %s -b %s/%s.bam > %s.bam\n' % (subsample_ratio_dict['%s_210' % sample_mono_210], nonsubsampled_bam_folder, sample_mono_210, sample_mono_210)

    qsub_file_mono_210 = '%s/qsub_subsample_%s.sh' % (outputs_folder, sample_mono_210)
    qsub_file_mono_210_handle = open(qsub_file_mono_210, 'w')
    qsub_file_mono_210_handle.write(header)
    qsub_file_mono_210_handle.write(module_lines)
    qsub_file_mono_210_handle.write('cd %s\n' % wd_on_katana)
    qsub_file_mono_210_handle.write(subsample_cmd_mono_210)
    qsub_file_mono_210_handle.close()



for sample_mono_D2 in sample_list_mono_D2:

    # subsample_cmd
    subsample_cmd_mono_D2 = 'samtools view -s %s -b %s/%s.bam > %s.bam' % (subsample_ratio_dict['%s_D2' % sample_mono_D2], nonsubsampled_bam_folder, sample_mono_D2, sample_mono_D2)

    qsub_file_mono_D2 = '%s/qsub_subsample_%s.sh' % (outputs_folder, sample_mono_D2)
    qsub_file_mono_D2_handle = open(qsub_file_mono_D2, 'w')
    qsub_file_mono_D2_handle.write(header)
    qsub_file_mono_D2_handle.write(module_lines)
    qsub_file_mono_D2_handle.write('cd %s\n' % wd_on_katana)
    qsub_file_mono_D2_handle.write(subsample_cmd_mono_D2)
    qsub_file_mono_D2_handle.close()


for sample_co in sample_list_co:

    # split sam cmd
    split_cmd_210c = 'samtools view -b %s/%s.bam 2.10_chromosome > %s_2.10_chromosome.bam\n' % (nonsubsampled_bam_folder, sample_co, sample_co)
    split_cmd_210p1 = 'samtools view -b %s/%s.bam 2.10_plasmid1 > %s_2.10_plasmid1.bam\n' % (nonsubsampled_bam_folder, sample_co, sample_co)
    split_cmd_210p2 = 'samtools view -b %s/%s.bam 2.10_plasmid2 > %s_2.10_plasmid2.bam\n' % (nonsubsampled_bam_folder, sample_co, sample_co)
    split_cmd_210p3 = 'samtools view -b %s/%s.bam 2.10_plasmid3 > %s_2.10_plasmid3.bam\n' % (nonsubsampled_bam_folder, sample_co, sample_co)
    split_cmd_D2c = 'samtools view -b %s/%s.bam D2_c > %s_D2_c.bam\n' % (nonsubsampled_bam_folder, sample_co, sample_co)
    split_cmd_D2p = 'samtools view -b %s/%s.bam D2_p > %s_D2_p.bam\n' % (nonsubsampled_bam_folder, sample_co, sample_co)

    # combine sam cmd
    merge_cmd_210 = 'samtools merge %s_2.10.bam %s_2.10_chromosome.bam %s_2.10_plasmid1.bam %s_2.10_plasmid2.bam %s_2.10_plasmid3.bam\n' % (sample_co, sample_co, sample_co, sample_co, sample_co)
    merge_cmd_D2 = 'samtools merge %s_D2.bam %s_D2_c.bam %s_D2_p.bam\n' % (sample_co, sample_co, sample_co)

    # subsample sam cmd
    subsample_cmd_co_210 = 'samtools view -s %s -b %s_2.10.bam > subsampled_%s_2.10.bam\n' % (subsample_ratio_dict['%s_210' % sample_co], sample_co, sample_co)
    subsample_cmd_co_D2 = 'samtools view -s %s -b %s_D2.bam > subsampled_%s_D2.bam\n' % (subsample_ratio_dict['%s_D2' % sample_co], sample_co, sample_co)

    # combine subsampled sam cmd
    merge_cmd_subsampled = 'samtools merge %s.bam subsampled_%s_2.10.bam subsampled_%s_D2.bam\n' % (sample_co, sample_co, sample_co)

    # index sam file
    index_cmd_subsampled = 'samtools index %s.bam\n' % (sample_co)



    qsub_file_co = '%s/qsub_subsample_%s.sh' % (outputs_folder, sample_co)
    qsub_file_co_handle = open(qsub_file_co, 'w')
    qsub_file_co_handle.write(header)
    qsub_file_co_handle.write(module_lines)
    qsub_file_co_handle.write('cd %s\n' % wd_on_katana)

    # split sam cmd
    qsub_file_co_handle.write(split_cmd_210c)
    qsub_file_co_handle.write(split_cmd_210p1)
    qsub_file_co_handle.write(split_cmd_210p2)
    qsub_file_co_handle.write(split_cmd_210p3)
    qsub_file_co_handle.write(split_cmd_D2c)
    qsub_file_co_handle.write(split_cmd_D2p)
    qsub_file_co_handle.write('\n')

    # combine sam cmd
    qsub_file_co_handle.write(merge_cmd_210)
    qsub_file_co_handle.write(merge_cmd_D2)
    qsub_file_co_handle.write('\n')

    # subsample sam cmd
    qsub_file_co_handle.write(subsample_cmd_co_210)
    qsub_file_co_handle.write(subsample_cmd_co_D2)
    qsub_file_co_handle.write('\n')

    # combine subsampled sam cmd
    qsub_file_co_handle.write(merge_cmd_subsampled)
    qsub_file_co_handle.write('\n')

    # index sam file
    qsub_file_co_handle.write(index_cmd_subsampled)
    qsub_file_co_handle.write('\n')

    # remove tmp files
    qsub_file_co_handle.write('rm %s_2.10_chromosome.bam\n' % sample_co)
    qsub_file_co_handle.write('rm %s_2.10_plasmid1.bam\n' % sample_co)
    qsub_file_co_handle.write('rm %s_2.10_plasmid2.bam\n' % sample_co)
    qsub_file_co_handle.write('rm %s_2.10_plasmid3.bam\n' % sample_co)
    qsub_file_co_handle.write('rm %s_D2_c.bam\n' % sample_co)
    qsub_file_co_handle.write('rm %s_D2_p.bam\n' % sample_co)
    qsub_file_co_handle.write('rm %s_2.10.bam\n' % sample_co)
    qsub_file_co_handle.write('rm %s_D2.bam\n' % sample_co)
    qsub_file_co_handle.write('rm subsampled_%s_2.10.bam\n' % sample_co)
    qsub_file_co_handle.write('rm subsampled_%s_D2.bam\n' % sample_co)

    qsub_file_co_handle.close()

