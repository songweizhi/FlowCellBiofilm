import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pylab
import numpy as np
from datetime import datetime


def get_refseq_id_and_len(ref_in):

    refseq_id_list = []
    refseq_len_dict = {}
    for each_seq in SeqIO.parse(ref_in, 'fasta'):
        refseq_id_list.append(each_seq.id)
        refseq_len_dict[each_seq.id] = len(each_seq.seq)

    return refseq_id_list, refseq_len_dict


def cigar_splitter(cigar):

    # get the position of letters
    letter_pos_list = []
    n = 0
    for each_element in cigar:
        if each_element.isalpha() == True:
            letter_pos_list.append(n)
        n += 1

    # split cigar
    index = 0
    cigar_splitted = []
    while index <= len(letter_pos_list) - 1:
        if index == 0:
            cigar_splitted.append(cigar[:(letter_pos_list[index] + 1)])
        else:
            cigar_splitted.append(cigar[(letter_pos_list[index - 1] + 1):(letter_pos_list[index] + 1)])
        index += 1

    return cigar_splitted


def get_cigar_seg_seq(cigar, cigar_seg, sequence):

    # get length and type of provided cigar segment
    cigar_seg_len = int(cigar_seg[:-1])
    cigar_seg_type = cigar_seg[-1]

    # split cigar at provided segment
    cigar_split_M = cigar.split('%s%s' % (cigar_seg_len, cigar_seg_type))

    # get left side length of provided segment
    left_len = 0
    if cigar_split_M[0] == '':
        left_len = 0
    else:
        left_split = cigar_splitter(cigar_split_M[0])
        for each in left_split:
            each_len = int(each[:-1])
            left_len += each_len

    # get sequence
    longest_M_seq = sequence[left_len: left_len + cigar_seg_len]

    # return sequence
    return longest_M_seq


def get_percent(list_in):

    list_sum = sum(list_in)
    list_p = []
    for each in list_in:
        each_p = each/list_sum*100
        each_p = float("{0:.2f}".format(each_p))
        list_p.append(each_p)

    return list_p


def unique_list_elements(list_input):

    list_output = []
    for each_element in list_input:
        if each_element not in list_output:
            list_output.append(each_element)

    return list_output


##################################################### CONFIGURATION ####################################################

parser = argparse.ArgumentParser()

parser.add_argument('-ref',
                    required=True,
                    help='reference sequences')

parser.add_argument('-sam',
                    required=True,
                    help='sam file')

args = vars(parser.parse_args())
refseq_file = args['ref']
sam_in = args['sam']

# input files
pwd_blastn_exe = 'blastn'

sequencing_depth = 23

################################################ extract clipping reads ################################################

print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Extract clipping reads')
output_reads_file = 'clipping_reads.fasta'
output_reads_file_handle = open(output_reads_file, 'w')
for each_aln in open(sam_in):
    if not each_aln.startswith('@'):
        each_split = each_aln.strip().split('\t')
        reads_id = each_split[0]

        reads_seq = each_split[9]
        reads_len = len(reads_seq)
        full_length_M = '%sM' % str(reads_len)
        cigar = each_split[5]

        # if not perfect match and CIGAR string available
        if (cigar != full_length_M) and (cigar != '*'):

            # split CIGAR string
            cigar_splitted = cigar_splitter(cigar)
            #print(cigar_splitted)

            # get longest aligned sequence (M) and soft clipping sequence (S)
            longest_M = 0
            longest_S = 0
            for each in cigar_splitted:

                if each[-1] == 'M':
                    if int(each[:-1]) > longest_M:
                        longest_M = int(each[:-1])

                if each[-1] == 'S':
                    if int(each[:-1]) > longest_S:
                        longest_S = int(each[:-1])

            # get M + S
            total_MS = longest_M + longest_S

            # if the percentage of matched and soft clipping sequences higher than 20% and their total length higher than 90%,
            # export the matched and soft clipping sequences separately.
            if (longest_M/reads_len >= 0.1) and (longest_S/reads_len >= 0.1) and ((longest_M + longest_S)/reads_len >= 0.8):

                # get cigar segment ID
                cigar_seg_longest_M = '%sM' % longest_M
                cigar_seg_longest_S = '%sS' % longest_S

                longest_M_seq = get_cigar_seg_seq(cigar, cigar_seg_longest_M, reads_seq)
                longest_S_seq = get_cigar_seg_seq(cigar, cigar_seg_longest_S, reads_seq)

                longest_M_seq_id = '%s_M' % reads_id
                longest_S_seq_id = '%s_S' % reads_id
                output_reads_file_handle.write('>%s\n' % longest_M_seq_id)
                output_reads_file_handle.write('%s\n' % longest_M_seq)
                output_reads_file_handle.write('>%s\n' % longest_S_seq_id)
                output_reads_file_handle.write('%s\n' % longest_S_seq)

output_reads_file_handle.close()


################################### blast clipping reads against reference sequences ###################################

# report
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Blast clipping reads against reference sequences')

# blast clipping reads against reference sequences
blast_output = 'blast_output.tab'
blast_parameters = '-evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn'
cmd_blast = '%s -query %s -subject %s -out %s %s' % (pwd_blastn_exe, output_reads_file, refseq_file, blast_output, blast_parameters)
os.system(cmd_blast)


################################################# parse blast results ##################################################

# report
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Get break points for each reference sequence')

# get refseq_id_list and refseq_len_dict
refseq_id_list, refseq_len_dict = get_refseq_id_and_len(refseq_file)

# initialize a dict to hold break points for each reference sequences
ctg_break_point_dict = {}
for each_ctg in refseq_id_list:
    ctg_break_point_dict[each_ctg] = []

# get ctg_break_point_dict
for match in open(blast_output):
    match_split = match.strip().split('\t')
    query = match_split[0]
    refseq_file = match_split[1]
    identity = float(match_split[2])
    align_len = int(match_split[3])
    query_len = int(match_split[12])
    subject_len = int(match_split[13])
    q_start = int(match_split[6])
    q_end = int(match_split[7])
    s_start = int(match_split[8])
    s_end = int(match_split[9])
    match_patern = []

    if (identity >= 99) and (align_len/query_len >= 0.9):
        ctg_break_point_dict[refseq_file].append(s_start)
        ctg_break_point_dict[refseq_file].append(s_end)

# export match ends
tmp1 = 'tmp1.txt'
tmp1_handle = open(tmp1, 'w')
for each_ref in ctg_break_point_dict:

    current_ref_pos_list = ctg_break_point_dict[each_ref]
    current_ref_pos_list_uniq = unique_list_elements(current_ref_pos_list)

    for each_uniq_pos in sorted(current_ref_pos_list_uniq):
        if current_ref_pos_list.count(each_uniq_pos) > sequencing_depth/10:
            for_out = '%s\t%s\t%s\n' % (each_ref, each_uniq_pos, current_ref_pos_list.count(each_uniq_pos))
            tmp1_handle.write(for_out)

tmp1_handle.close()

# Combine continuous ends (potential insertion point)
tmp2 = 'tmp2.txt'
tmp2_handle = open(tmp2, 'w')
current_id = ''
current_seq = ''
current_loc = 0
currenr_num = 0
for each_pos in open(tmp1):

    each_pos_split = each_pos.strip().split('\t')
    each_pos_seq = each_pos_split[0]
    each_pos_loc = int(each_pos_split[1])
    each_pos_num = int(each_pos_split[2])

    if current_loc == 0:
        current_seq = each_pos_seq
        current_loc = each_pos_loc
        currenr_num = each_pos_num
        current_id = each_pos_loc

    elif (current_seq == each_pos_seq) and (each_pos_loc > current_loc + 1):

        if '-' in str(current_id):
            tmp2_handle.write('%s\t%s\t%s\n' % (current_seq, current_id, currenr_num))

        # reset
        current_seq = each_pos_seq
        current_loc = each_pos_loc
        currenr_num = each_pos_num
        current_id = each_pos_loc

    elif (current_seq == each_pos_seq) and (each_pos_loc < current_loc + 5):

        current_loc = each_pos_loc
        currenr_num += each_pos_num
        current_id = '%s-%s' % (current_id, each_pos_loc)

    elif current_seq != each_pos_seq:

        if '-' in str(current_id):
            tmp2_handle.write('%s\t%s\t%s\n' % (current_seq, current_id, currenr_num))

        # reset
        current_seq = each_pos_seq
        current_loc = each_pos_loc
        currenr_num = each_pos_num
        current_id = each_pos_loc

if '-' in str(current_id):
    tmp2_handle.write('%s\t%s\t%s\n' % (current_seq, current_id, currenr_num))

tmp2_handle.close()


# report
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Done!')






