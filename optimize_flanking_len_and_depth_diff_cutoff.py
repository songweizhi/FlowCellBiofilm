import os
import glob


wd = '/Users/songweizhi/Desktop/4_2_SNV_QC'
QC_result_folders = ['output_f50000bp_1000mer_dl500bp', 'output_f50000bp_1000mer_dl1000bp', 'output_f50000bp_1000mer_dl2000bp', 'output_f50000bp_1000mer_dl5000bp']
cutoff = 15

os.chdir(wd)


for QC_results in QC_result_folders:

    folder_lt_20 = '%s_higher_%s' % (QC_results, cutoff)
    folder_st_20 = '%s_lower_%s' % (QC_results, cutoff)

    os.mkdir(folder_lt_20)
    os.mkdir(folder_st_20)

    png_files_re = '%s/SNV_depth_plot/*.png' % (QC_results)
    png_file_list = [os.path.basename(file_name) for file_name in glob.glob(png_files_re)]

    for each_png in png_file_list:

        pwd_each_png = '%s/SNV_depth_plot/%s' % (QC_results, each_png)
        each_png_depth_diff = float(each_png.split('_')[-1][:-4])
        each_png_newname = '%s_%s' % (each_png_depth_diff, each_png)

        if each_png_depth_diff > cutoff:
            os.system('cp %s %s/%s' % (pwd_each_png, folder_lt_20, each_png_newname))

        if each_png_depth_diff < cutoff:
            os.system('cp %s %s/%s' % (pwd_each_png, folder_st_20, each_png_newname))


