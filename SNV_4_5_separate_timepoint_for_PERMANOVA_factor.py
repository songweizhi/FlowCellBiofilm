import os


#wd = '/Users/songweizhi/Dropbox/Research/Flow_cell/deepSNV_output_summary'
wd = '/Users/songweizhi/Desktop/Permanova'

timepoint_list = ['D9', 'D18', 'D27', 'D42']

infile = 'deepSNV_output_summary_210_factor.txt'
#infile = 'deepSNV_output_summary_D2_factor.txt'

os.chdir(wd)

infile_basename, ext = os.path.splitext(infile)

for timepoint in timepoint_list:

    out_file_name = '%s_%s%s' % (infile_basename, timepoint, ext)
    out_file_handle = open(out_file_name, 'w')
    out_file_handle.write('Sample\tSpecies\tReplicate\tTime\tLabel\n')

    for each in open(infile):
        each_split = each.strip().split('\t')

        if each_split[0].split('D')[-1] == timepoint[1:]:
            out_file_handle.write(each)

    out_file_handle.close()
