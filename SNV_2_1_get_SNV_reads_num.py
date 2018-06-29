import os
import glob
import pysam
import argparse


# python3 ~/PycharmProjects/FlowCellBiofilm/SNV_2_1_deepSNV_remove_false_positive.py -bf /Users/songweizhi/Desktop/pysam_wd -df /Users/songweizhi/Dropbox/Research/Flow_cell/deepSNV_output_summary/4_deepSNV_output -i 8D27


def get_depth_and_snv_reads_num(input_bam_file, detected_snv):

    detected_snv_split = detected_snv.split('|')
    refseq_id = detected_snv_split[0]
    snv_pos = int(detected_snv_split[1])
    snv_pos_wt = detected_snv_split[2]
    snv_pos_v = detected_snv_split[3]

    samfile = pysam.AlignmentFile(input_bam_file, "rb")

    current_pos_depth = 0
    snv_pos_v_num = 0

    for pileupcolumn in samfile.pileup(refseq_id):

        if pileupcolumn.pos + 1 == snv_pos:
            mapped_ncs = []
            for pileupread in pileupcolumn.pileups:
                current_bp = ''
                if pileupread.is_del:
                    current_bp = '-'
                elif pileupread.is_refskip:
                    current_bp = 'S'
                else:
                    current_bp = pileupread.alignment.query_sequence[pileupread.query_position]
                mapped_ncs.append(current_bp)

            current_pos_depth = pileupcolumn.n
            snv_pos_v_num = mapped_ncs.count(snv_pos_v)

    samfile.close()

    return current_pos_depth, snv_pos_v_num


################################################# input #################################################

parser = argparse.ArgumentParser(description='', add_help=False)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')

optional.add_argument('-h', action='help', help='Show this help message and exit')
required.add_argument('-bf', dest='Bam_folder', nargs='?', required=True,  type=str, help='reference sequence file')
required.add_argument('-df', dest='deepSNV_folder', nargs='?', required=True,  type=str, help='depth file')
required.add_argument('-i', dest='Sample_ID', nargs='?', required=False, type=str, default=None, help='id of sequence to plot')

args = vars(parser.parse_args())
bam_file_folder = args['Bam_folder']
deepSNV_output_folder = args['deepSNV_folder']
sample_id = args['Sample_ID']


#########################################################################################################

# output file'
output_file = '%s_snv_reads_num.txt' % sample_id


# get deepSNV output file list
deepSNV_output_file_re = '%s/%s_*.txt' % (deepSNV_output_folder, sample_id)
deepSNV_output_file_list = [os.path.basename(file_name) for file_name in glob.glob(deepSNV_output_file_re)]


output_file_handle = open(output_file, 'w')
for each_deepSNV_output in deepSNV_output_file_list:
    pwd_each_deepSNV_output = '%s/%s' % (deepSNV_output_folder, each_deepSNV_output)
    treatment_id = sample_id
    print('Processing %s' % each_deepSNV_output)
    for each_snv in open(pwd_each_deepSNV_output):

        # ignore the first line
        if not each_snv.startswith('chr,pos,ref,var'):
            each_snv_split = each_snv.strip().split(',')
            snv = '|'.join(each_snv_split[0:4])
            pwd_bam_file = '%s/%s.bam' % (bam_file_folder, treatment_id)

            # run get_depth_and_snv_reads_num
            snv_depth, snv_reads_count = get_depth_and_snv_reads_num(pwd_bam_file, snv)

            if snv_reads_count < 10:
                output_file_handle.write('%s\t%s\t%s\t%s\t%s\n' % (treatment_id, snv, snv_depth, snv_reads_count, 'FP'))

            elif snv_reads_count >= 10:
                output_file_handle.write('%s\t%s\t%s\t%s\t%s\n' % (treatment_id, snv, snv_depth, snv_reads_count, 'P'))

output_file_handle.close()
