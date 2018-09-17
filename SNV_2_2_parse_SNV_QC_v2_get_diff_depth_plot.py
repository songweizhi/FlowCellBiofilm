import os
import glob

############################################### input file and parameters ##############################################

# wd
wd = '/Users/songweizhi/Dropbox/Research/Flow_cell_datasets/parse_SNV_QC_v2'
os.chdir(wd)

diff_depth_cutoff = 15
depth_plots_folder_all = 'SNV_depth_plot'
depth_plots_folder_diff = '%s/SNV_depth_plot_diff_%s' % (wd, diff_depth_cutoff)
depth_plots_folder_even = '%s/SNV_depth_plot_even_%s' % (wd, diff_depth_cutoff)

depth_plot_files = '%s/%s/*.png' % (wd, depth_plots_folder_all)
depth_plot_file_list = [os.path.basename(file_name) for file_name in glob.glob(depth_plot_files)]

for each_plot in depth_plot_file_list:
    pwd_plot = '%s/%s/%s' % (wd, depth_plots_folder_all, each_plot)
    depth_difference = float('.'.join(each_plot.split('_')[-1].split('.')[:-1]))

    if depth_difference > diff_depth_cutoff:
        cmd = 'cp %s %s/' % (pwd_plot, depth_plots_folder_diff)
    else:
        cmd = 'cp %s %s/' % (pwd_plot, depth_plots_folder_even)

    os.system(cmd)


