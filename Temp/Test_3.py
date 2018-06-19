import os
import glob

os.chdir('/Users/songweizhi/Desktop')
infile_re = 'in/*.txt'
infile_list = [os.path.basename(file_name) for file_name in glob.glob(infile_re)]


for each_file in infile_list:

    pwd_each_file_in = 'in/%s' % each_file
    pwd_each_file_out = 'out/%s' % each_file

    each_file_out_handle = open(pwd_each_file_out, 'w')
    for each in open(pwd_each_file_in):
        each_split = each.strip().split('\t')
        if not each_split[1:] == ['0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0']:
            each_file_out_handle.write(each)
    each_file_out_handle.close()
