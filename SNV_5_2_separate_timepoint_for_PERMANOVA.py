import os


#wd = '/Users/songweizhi/Dropbox/Research/Flow_cell/deepSNV_output_summary'
wd = '/Users/songweizhi/Desktop/aaa'

timepoint_list = ['D9', 'D18', 'D27', 'D42']

infile = 'deepSNV_output_summary_210_existence_cdc.txt'
#infile = 'deepSNV_output_summary_210_frequency_cdc.txt'
#infile = 'deepSNV_output_summary_D2_existence_cdc.txt'
#infile = 'deepSNV_output_summary_D2_frequency_cdc.txt'


os.chdir(wd)

infile_basename, ext = os.path.splitext(infile)

for timepoint in timepoint_list:

    out_file_name = '%s_%s%s' % (infile_basename, timepoint, ext)
    out_file_handle = open(out_file_name, 'w')

    for each in open(infile):
        each_split = each.strip().split('\t')

        needed = ''
        if timepoint == 'D9':
            if each.startswith('	'):
                needed = '%s\t%s\t%s\t%s\t%s\t%s' % (each_split[0], each_split[4], each_split[8], each_split[12], each_split[16], each_split[20])
            else:
                needed = '%s\t%s\t%s\t%s\t%s\t%s\t%s' % (each_split[0], each_split[1], each_split[5], each_split[9], each_split[13], each_split[17], each_split[21])

        if timepoint == 'D18':
            if each.startswith('	'):
                needed = '%s\t%s\t%s\t%s\t%s\t%s' % (each_split[1], each_split[5], each_split[9], each_split[13], each_split[17], each_split[21])
            else:
                needed = '%s\t%s\t%s\t%s\t%s\t%s\t%s' % (each_split[0], each_split[2], each_split[6], each_split[10], each_split[14], each_split[18], each_split[22])

        if timepoint == 'D27':
            if each.startswith('	'):
                needed = '%s\t%s\t%s\t%s\t%s\t%s' % (each_split[2], each_split[6], each_split[10], each_split[14], each_split[18], each_split[22])
            else:
                needed = '%s\t%s\t%s\t%s\t%s\t%s\t%s' % (each_split[0], each_split[3], each_split[7], each_split[11], each_split[15], each_split[19], each_split[23])

        if timepoint == 'D42':
            print(each_split)
            if each.startswith('	'):
                needed = '%s\t%s\t%s\t%s\t%s\t%s' % (each_split[3], each_split[7], each_split[11], each_split[15], each_split[19], each_split[23])
            else:
                needed = '%s\t%s\t%s\t%s\t%s\t%s\t%s' % (each_split[0], each_split[4], each_split[8], each_split[12], each_split[16], each_split[20], each_split[24])

        needed_split = needed.split('\t')
        if not needed_split[1:] == ['0', '0', '0', '0', '0', '0']:
            out_file_handle.write('%s\n' % needed)

    out_file_handle.close()


