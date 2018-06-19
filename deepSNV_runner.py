import os
import argparse
from Bio import SeqIO
from datetime import datetime


# example command
# python ~/PycharmProjects/Flow_cell_biofilm/deep_SNV_runner.py -r D2_pacbio.fasta -q 2D9.bam -c D2D0.bam


##################################################### CONFIGURATION ####################################################

parser = argparse.ArgumentParser()

parser.add_argument('-r',
                    required=True,
                    help='reference sequences')

parser.add_argument('-q',
                    required=True,
                    help='query bam file')

parser.add_argument('-c',
                    required=True,
                    help='control bam file')

parser.add_argument('-R',
                    required=False,
                    default='~/R_scripts/deepSNV/get_SNVs.R',
                    help='get_SNVs.R')

args = vars(parser.parse_args())
ref = args['r']
query_bam = args['q']
control_bam = args['c']
get_SNVs_R = args['R']


########################################################## MAIN ########################################################

# get file name
query_bam_file_path, query_bam_file_name =      os.path.split(query_bam)
control_bam_file_path, control_bam_file_name =  os.path.split(control_bam)
query_bam_basename, query_bam_extension =       os.path.splitext(query_bam_file_name)
control_bam_basename, control_bam_extension =   os.path.splitext(control_bam_file_name)


# get seq_id_to_length_dict
seq_id_to_length_dict = {}
for each in SeqIO.parse(ref, 'fasta'):
    seq_id_to_length_dict[each.id] = len(each.seq)


# run deepSNV for all reference sequences
start_pos = 1
for each_seq in seq_id_to_length_dict:
    each_seq_len = seq_id_to_length_dict[each_seq]
    R_cmd = 'Rscript %s -q %s -c %s -i %s -s %s -e %s' % (get_SNVs_R, query_bam, control_bam, each_seq, start_pos, each_seq_len)
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Running: ' + R_cmd)
    os.system(R_cmd)

    # report
    output_file = '%s_vs_%s_%s_%s_%s_SNVs.txt' % (query_bam_basename, control_bam_basename, each_seq, start_pos, each_seq_len)
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Output exported to: ' + output_file)

