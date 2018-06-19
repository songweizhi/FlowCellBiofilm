import os


wd = '/Users/songweizhi/Desktop/FC'
deepSNV_output = 'deepSNV_output_summary_all_existence.txt'
deepSNV_output_cdc = 'deepSNV_output_summary_all_existence_cdc.txt'

os.chdir(wd)


deepSNV_output_cdc_handle = open(deepSNV_output_cdc, 'w')
current_seq = ''
current_pos_start = 0
current_profile_start = ''
current_pos = 0
deleted_ncs = ''
for each_snv in open(deepSNV_output):
    if not each_snv.startswith('\t'):
        each_snv_seq = each_snv.strip().split('\t')[0].split('|')[0]
        each_snv_pos = int(each_snv.strip().split('\t')[0].split('|')[1])
        each_snv_pos_wt = each_snv.strip().split('\t')[0].split('|')[2]
        each_snv_pos_v = each_snv.strip().split('\t')[0].split('|')[3]

        if (each_snv_pos_v == '-') and (current_seq == '') and (current_pos == 0) and (deleted_ncs == ''):
            current_seq = each_snv_seq
            current_pos_start = each_snv_pos
            current_profile_start = '\t'.join(each_snv.strip().split('\t')[1:])
            current_pos = each_snv_pos
            deleted_ncs = each_snv_pos_wt
        elif (each_snv_seq == current_seq) and (each_snv_pos == (current_pos + 1)) and (each_snv_pos_v != '-'):
            for_out = '%s|%s-%s|%s|-\t%s\n' % (current_seq, current_pos_start, current_pos, deleted_ncs, current_profile_start)
            deepSNV_output_cdc_handle.write(for_out)
            deepSNV_output_cdc_handle.write(each_snv)
            current_seq = ''
            current_pos_start = 0
            current_profile_start = ''
            current_pos = 0
            deleted_ncs = ''

        elif (each_snv_seq == current_seq) and (each_snv_pos == (current_pos + 1)) and (each_snv_pos_v == '-'):
            current_pos = each_snv_pos
            deleted_ncs += each_snv_pos_wt
        elif each_snv_pos != (current_pos + 1):
            if (current_seq != '') and (current_pos != 0) and (deleted_ncs != ''):
                if len(deleted_ncs) == 1:
                    for_out = '%s|%s|%s|-\t%s\n' % (current_seq, current_pos, deleted_ncs, current_profile_start)
                elif len(deleted_ncs) > 1:
                    for_out = '%s|%s-%s|%s|-\t%s\n' % (current_seq, current_pos_start, current_pos, deleted_ncs, current_profile_start)
                deepSNV_output_cdc_handle.write(for_out)

                if each_snv_pos_v == '-':
                    current_seq = each_snv_seq
                    current_pos_start = each_snv_pos
                    current_profile_start = '\t'.join(each_snv.strip().split('\t')[1:])
                    current_pos = each_snv_pos
                    deleted_ncs = each_snv_pos_wt
                else:
                    deepSNV_output_cdc_handle.write(each_snv)
                    current_seq = ''
                    current_pos_start = 0
                    current_profile_start = ''
                    current_pos = 0
                    deleted_ncs = ''

            elif (current_seq == '') and (current_pos == 0) and (deleted_ncs == '') and (each_snv_pos_v != '-'):
                deepSNV_output_cdc_handle.write(each_snv)

deepSNV_output_cdc_handle.close()

