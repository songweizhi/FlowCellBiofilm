import os


wd = '/Users/songweizhi/Desktop/FC'
os.chdir(wd)

plot_sam_depth_script = '~/PycharmProjects/FlowCellBiofilm/plot_sam_depth.py'
snv_cdc_file = 'deepSNV_output_summary_all_frequency_cdc.txt'
flanking_length = 10000
kmer = 1000


summary_file_header_all = ['1D9', '1D18', '1D27', '1D42', '5D9', '5D18', '5D27', '5D42', '9D9', '9D18', '9D27', '9D42', '2D9', '2D18', '2D27', '2D42', '6D9', '6D18', '6D27', '6D42', '10D9', '10D18', '10D27', '10D42', '4D9', '4D18', '4D27', '4D42', '8D9', '8D18', '8D27', '8D42', '12D9', '12D18', '12D27', '12D42']
reference_dict = {'1': '2.10wt_illumina.fasta',
                  '5': '2.10wt_illumina.fasta',
                  '9': '2.10wt_illumina.fasta',
                  '2': 'D2_pacbio.fasta',
                  '6': 'D2_pacbio.fasta',
                  '10': 'D2_pacbio.fasta',
                  '4': 'combined_ref.fasta',
                  '8': 'combined_ref.fasta',
                  '12': 'combined_ref.fasta'}


for each_snv in open(snv_cdc_file):
    if not each_snv.startswith('	'):
        each_snv_id = each_snv.strip().split('\t')[0]
        each_snv_id_split = each_snv_id.split('|')
        each_snv_pos = each_snv_id_split[1]
        ref_seq_id = each_snv_id_split[0]
        time_points = each_snv.strip().split('\t')[1:]

        n = 0
        for each_time_point in time_points:

            if each_time_point != '0':
                current_time_point = summary_file_header_all[n]
                depth_file = '%s.depth' % current_time_point
                reference_file = reference_dict[current_time_point.split('D')[0]]
                reference_seq = ref_seq_id
                plot_name = '%s_%s_f%sbp_%smer' % (each_snv_id, current_time_point, flanking_length, kmer)
                cmd = ''

                if '-' in each_snv_pos:
                    each_snv_pos_start = int(each_snv_pos.split('-')[0])
                    each_snv_pos_end = int(each_snv_pos.split('-')[1])
                    plot_start = each_snv_pos_start - flanking_length
                    plot_end = each_snv_pos_end + flanking_length
                    cmd = 'python3 %s -r %s -d %s -i %s -s %s -e %s -k %s -o %s' % (plot_sam_depth_script, reference_file, depth_file, reference_seq, plot_start, plot_end, kmer, plot_name)
                else:
                    plot_start = int(each_snv_pos) - flanking_length
                    plot_end = int(each_snv_pos) + flanking_length
                    cmd = 'python3 %s -r %s -d %s -i %s -s %s -e %s -k %s -o %s' % (plot_sam_depth_script, reference_file, depth_file, reference_seq, plot_start, plot_end, kmer, plot_name)

                print(cmd)

            n += 1

