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


def export_dna_record(gene_seq, gene_id, gene_description, output_handle):

    seq_object = Seq(gene_seq, IUPAC.unambiguous_dna)
    seq_record = SeqRecord(seq_object)
    seq_record.id = gene_id
    seq_record.description = gene_description
    SeqIO.write(seq_record, output_handle, 'fasta')


def get_percent(list_in):

    list_sum = sum(list_in)
    list_p = []
    for each in list_in:
        each_p = each/list_sum*100
        each_p = float("{0:.2f}".format(each_p))
        list_p.append(each_p)

    return list_p


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


################################################ extract clipping reads ################################################

print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Extract clipping reads')
output_reads_file = 'clipping_reads.fasta'
output_reads_file_handle = open(output_reads_file, 'w')
for each in open(sam_in):
    if not each.startswith('@'):
        each_split = each.strip().split('\t')
        reads_id = each_split[0]

        reads_seq = each_split[9]
        reads_len = len(reads_seq)
        full_length_M = '%sM' % str(reads_len)
        cigar = each_split[5]

        if cigar != full_length_M:
            if cigar == '*':
                # extract reads
                output_reads_file_handle.write('>%s\n' % reads_id)
                output_reads_file_handle.write('%s\n' % reads_seq)
            else:
                cigar_splitted = cigar_splitter(cigar)
                # get longest M
                current_M = 0
                for each in cigar_splitted:
                    if each[-1] == 'M':
                        if int(each[:-1]) > current_M:
                            current_M = int(each[:-1])
                percent_longest_M = current_M/reads_len

                if percent_longest_M >= 0.75:
                    # extract reads
                    output_reads_file_handle.write('>%s\n' % reads_id)
                    output_reads_file_handle.write('%s\n' % reads_seq)

output_reads_file_handle.close()


################################### blast clipping reads against reference sequences ###################################

print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Blast clipping reads against reference sequences')
# blast clipping reads against reference sequences
blast_output = 'blast_output.tab'
blast_parameters = '-evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn'
cmd_blast = '%s -query %s -subject %s -out %s %s' % (pwd_blastn_exe, output_reads_file, refseq_file, blast_output, blast_parameters)
os.system(cmd_blast)


########################### get break point list for each reference sequence and store in dict #########################

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
    if identity >= 99:
        if q_start == 1:
            ctg_break_point_dict[refseq_file].append(s_end)
        if q_end == query_len:
            ctg_break_point_dict[refseq_file].append(s_start)


# get plot
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Get plot')
for each_ref in refseq_id_list:
    break_points = ctg_break_point_dict[each_ref]
    break_points_sorted = sorted(break_points)
    break_points_uniq = []
    for each_pos in break_points_sorted:
        if each_pos not in break_points_uniq:
            break_points_uniq.append(each_pos)

    break_points_uniq_count = []
    for each_uniq_pos in break_points_uniq:
        each_uniq_pos_count = break_points_sorted.count(each_uniq_pos)
        break_points_uniq_count.append(each_uniq_pos_count)

    # turn count number into percent
    break_points_uniq_count_p = get_percent(break_points_uniq_count)

    # plot
    pwd_image_file = 'plot_%s.png' % each_ref
    x = np.array(break_points_uniq)
    y = np.array(break_points_uniq_count_p)
    pylab.scatter(x, y, s=6)
    plt.title('Reference ID: %s' % each_ref)
    plt.xlabel('Position (bp)')
    plt.ylabel('Percent (%)')
    plt.axis([0, refseq_len_dict[each_ref], 0, None])  # set the range of X/Y-axis
    plt.tight_layout(True)
    plt.savefig(pwd_image_file, dpi=300)
    plt.close()

    # print(break_points_uniq)
    # print(break_points_uniq_count)
