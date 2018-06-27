import os


wd = '/Users/songweizhi/Desktop/FC0'

deepSNV_output = 'deepSNV_output_summary_all_frequency.txt'
#deepSNV_output = 'deepSNV_output_summary_210_existence.txt'
#deepSNV_output = 'deepSNV_output_summary_210_frequency.txt'
#deepSNV_output = 'deepSNV_output_summary_D2_existence.txt'
#deepSNV_output = 'deepSNV_output_summary_D2_frequency.txt'


# prepare output file name
deepSNV_output_basename, ext = os.path.splitext(deepSNV_output)
deepSNV_output_cdc = '%s_cdc.txt' % deepSNV_output_basename

os.chdir(wd)

deepSNV_output_cdc_handle = open(deepSNV_output_cdc, 'w')
current_seq = ''
current_pos_start = 0
current_profile_start = ''
current_pos = 0
changed_ncs_WT = ''
changed_ncs_V = ''
current_mutation_id = ''
processed_mutations = []
for each_snv in open(deepSNV_output):
    if each_snv.startswith('\t'):
        deepSNV_output_cdc_handle.write(each_snv)

    if not each_snv.startswith('\t'):
        each_snv_seq = each_snv.strip().split('\t')[0].split('|')[0]
        each_snv_pos = int(each_snv.strip().split('\t')[0].split('|')[1])
        each_snv_pos_wt = each_snv.strip().split('\t')[0].split('|')[2]
        each_snv_pos_v = each_snv.strip().split('\t')[0].split('|')[3]
        mutation_id = '%s|%s' % (each_snv_seq, each_snv_pos)


        if (current_seq == '') and (current_pos == 0) and (changed_ncs_WT == '') and (changed_ncs_V == ''):
            current_seq = each_snv_seq
            current_pos = each_snv_pos
            current_pos_start = each_snv_pos
            changed_ncs_WT = each_snv_pos_wt
            changed_ncs_V = each_snv_pos_v
            current_mutation_id = mutation_id
            current_profile_start = '\t'.join(each_snv.strip().split('\t')[1:])

        elif (each_snv_seq == current_seq) and (each_snv_pos == current_pos + 1):

            # if current_mutation_id in processed_mutations:
            #     deepSNV_output_cdc_handle.write('%s|%s|%s|%s\n' % (current_seq, current_pos, changed_ncs_WT, changed_ncs_V))
            #     processed_mutations.append(current_mutation_id)
            #
            # else:
            #     current_pos = each_snv_pos
            #     changed_ncs_WT += each_snv_pos_wt
            #     changed_ncs_V += each_snv_pos_v

            current_pos = each_snv_pos
            changed_ncs_WT += each_snv_pos_wt
            changed_ncs_V += each_snv_pos_v


        elif (each_snv_seq != current_seq) or (each_snv_pos != current_pos + 1):

            # write out
            if len(changed_ncs_WT) == 1:
                #deepSNV_output_cdc_handle.write('%s|%s|%s|%s\t%s\n' % (current_seq, current_pos, changed_ncs_WT, changed_ncs_V, current_profile_start))
                deepSNV_output_cdc_handle.write('%s|%s|%s|%s\n' % (current_seq, current_pos, changed_ncs_WT, changed_ncs_V))
                processed_mutations.append(current_mutation_id)

            elif len(changed_ncs_WT) > 1:
                #deepSNV_output_cdc_handle.write('%s|%s-%s|%s|%s\t%s\n' % (current_seq, current_pos_start, current_pos, changed_ncs_WT, changed_ncs_V, current_profile_start))
                deepSNV_output_cdc_handle.write('%s|%s-%s|%s|%s\n' % (current_seq, current_pos_start, current_pos, changed_ncs_WT, changed_ncs_V))
                processed_mutations.append(current_mutation_id)

            # reset variables
            current_seq = each_snv_seq
            current_pos = each_snv_pos
            current_pos_start = each_snv_pos
            changed_ncs_WT = each_snv_pos_wt
            changed_ncs_V = each_snv_pos_v
            current_profile_start = '\t'.join(each_snv.strip().split('\t')[1:])

# write out
if len(changed_ncs_WT) == 1:
#    deepSNV_output_cdc_handle.write('%s|%s|%s|%s\t%s\n' % (current_seq, current_pos, changed_ncs_WT, changed_ncs_V, current_profile_start))
    deepSNV_output_cdc_handle.write('%s|%s|%s|%s\n' % (current_seq, current_pos, changed_ncs_WT, changed_ncs_V))
    processed_mutations.append(current_mutation_id)

elif len(changed_ncs_WT) > 1:
#    deepSNV_output_cdc_handle.write('%s|%s-%s|%s|%s\t%s\n' % (current_seq, current_pos_start, current_pos, changed_ncs_WT, changed_ncs_V, current_profile_start))
    deepSNV_output_cdc_handle.write('%s|%s-%s|%s|%s\n' % (current_seq, current_pos_start, current_pos, changed_ncs_WT, changed_ncs_V))
    processed_mutations.append(current_mutation_id)

#         if (each_snv_pos_v == '-') and (current_seq == '') and (current_pos == 0) and (deleted_ncs == ''):
#             current_seq = each_snv_seq
#             current_pos_start = each_snv_pos
#             current_profile_start = '\t'.join(each_snv.strip().split('\t')[1:])
#             current_pos = each_snv_pos
#             deleted_ncs = each_snv_pos_wt
#         elif (each_snv_seq == current_seq) and (each_snv_pos == (current_pos + 1)) and (each_snv_pos_v != '-'):
#             for_out = '%s|%s-%s|%s|-\t%s\n' % (current_seq, current_pos_start, current_pos, deleted_ncs, current_profile_start)
#             deepSNV_output_cdc_handle.write(for_out)
#             deepSNV_output_cdc_handle.write(each_snv)
#             current_seq = ''
#             current_pos_start = 0
#             current_profile_start = ''
#             current_pos = 0
#             deleted_ncs = ''
#
#         elif (each_snv_seq == current_seq) and (each_snv_pos == (current_pos + 1)) and (each_snv_pos_v == '-'):
#             current_pos = each_snv_pos
#             deleted_ncs += each_snv_pos_wt
#         elif each_snv_pos != (current_pos + 1):
#             if (current_seq != '') and (current_pos != 0) and (deleted_ncs != ''):
#                 for_out = ''
#                 if len(deleted_ncs) == 1:
#                     for_out = '%s|%s|%s|-\t%s\n' % (current_seq, current_pos, deleted_ncs, current_profile_start)
#                 elif len(deleted_ncs) > 1:
#                     for_out = '%s|%s-%s|%s|-\t%s\n' % (current_seq, current_pos_start, current_pos, deleted_ncs, current_profile_start)
#                 deepSNV_output_cdc_handle.write(for_out)
#
#                 if each_snv_pos_v == '-':
#                     current_seq = each_snv_seq
#                     current_pos_start = each_snv_pos
#                     current_profile_start = '\t'.join(each_snv.strip().split('\t')[1:])
#                     current_pos = each_snv_pos
#                     deleted_ncs = each_snv_pos_wt
#                 else:
#                     deepSNV_output_cdc_handle.write(each_snv)
#                     current_seq = ''
#                     current_pos_start = 0
#                     current_profile_start = ''
#                     current_pos = 0
#                     deleted_ncs = ''
#
#             elif (current_seq == '') and (current_pos == 0) and (deleted_ncs == '') and (each_snv_pos_v != '-'):
#                 deepSNV_output_cdc_handle.write(each_snv)
#
# deepSNV_output_cdc_handle.close()

