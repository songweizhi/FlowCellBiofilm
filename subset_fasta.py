import argparse
from Bio import SeqIO


parser = argparse.ArgumentParser()

parser.add_argument('-i', required=True, help='fasta in ')
parser.add_argument('-o', required=True, help='fasta out')
parser.add_argument('-s', required=True, type=int, help='subset step')

args = vars(parser.parse_args())

fasta_in    = args['i']
fasta_out   = args['o']
subset_step = args['s']


fasta_out_handle = open(fasta_out, 'w')
n = 1
for seq_record in SeqIO.parse(fasta_in, 'fasta'):

    if n%subset_step == 0:
        fasta_out_handle.write('>%s\n' % seq_record.id)
        fasta_out_handle.write('%s\n' % str(seq_record.seq))

    n += 1

fasta_out_handle.close()

