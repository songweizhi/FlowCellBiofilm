import os
import glob


# define input
depth_file_folder = '/Users/songweizhi/Desktop/depth_files_446'
depth_file_ext = 'depth'
output_depth_summary = '/Users/songweizhi/Desktop/depth_summary_446.txt'


# get file list
depth_file_re = '%s/*.%s' % (depth_file_folder, depth_file_ext)
depth_file_list = [os.path.basename(file_name) for file_name in glob.glob(depth_file_re)]


# output depth summary
output_depth_summary_handle = open(output_depth_summary, 'w')

for depth_file in depth_file_list:

    pwd_depth_file = '%s/%s' % (depth_file_folder, depth_file)

    # get total length and depth
    total_len_210 = 0
    total_dep_210 = 0
    total_len_D2 = 0
    total_dep_D2 = 0
    for each_pos in open(pwd_depth_file):
        each_pos_split = each_pos.strip().split('\t')
        each_pos_ref = each_pos_split[0]
        each_pos_loc = each_pos_split[1]
        each_pos_dep = int(each_pos_split[2])

        # get total length and depth for 2.10
        if each_pos_ref.startswith('2.10'):
            total_len_210 += 1
            total_dep_210 += each_pos_dep

        # get total length and depth for D2
        if each_pos_ref.startswith('D2'):
            total_len_D2 += 1
            total_dep_D2 += each_pos_dep

    # report
    if depth_file in ['1D18.depth', '1D27.depth', '1D42.depth', '1D9.depth',
                      '5D18.depth', '5D27.depth', '5D42.depth', '5D9.depth',
                      '9D18.depth', '9D27.depth', '9D42.depth', '9D9.depth', '210WTD0.depth']:

        average_depth_210 = float("{0:.1f}".format(total_dep_210 / total_len_210))
        output_depth_summary_handle.write('Monoculture\t%s\t210\t%s\n' % (depth_file, average_depth_210))

    if depth_file in ['2D18.depth', '2D27.depth', '2D42.depth', '2D9.depth',
                      '6D18.depth', '6D27.depth', '6D42.depth', '6D9.depth',
                      '10D18.depth', '10D27.depth', '10D42.depth', '10D9.depth', 'D2D0.depth']:

        average_depth_D2 = float("{0:.1f}".format(total_dep_D2 / total_len_D2))
        output_depth_summary_handle.write('Monoculture\t%s\tD2\t%s\n' % (depth_file, average_depth_D2))

    if depth_file in ['4D18.depth', '4D27.depth', '4D42.depth', '4D9.depth',
                      '8D18.depth', '8D27.depth', '8D42.depth', '8D9.depth',
                      '12D18.depth', '12D27.depth', '12D42.depth', '12D9.depth', 'coculture_D0.depth']:

        average_depth_210 = float("{0:.1f}".format(total_dep_210 / total_len_210))
        average_depth_D2 = float("{0:.1f}".format(total_dep_D2 / total_len_D2))
        output_depth_summary_handle.write('Coculture\t%s\t210\t%s\n' % (depth_file, average_depth_210))
        output_depth_summary_handle.write('Coculture\t%s\tD2\t%s\n' % (depth_file, average_depth_D2))

output_depth_summary_handle.close()

