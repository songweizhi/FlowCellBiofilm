import pysam
import os
import argparse
from Bio import SeqIO
from datetime import datetime


def get_reads_mapped_to_a_region(sam_file, seq_id, start_pos, end_pos):

    samfile = pysam.AlignmentFile(sam_file, "rb")
    iter = samfile.fetch(seq_id, start_pos, end_pos)

    mapped_read_id_list = []
    for mapped_read in iter:
        read_id = mapped_read.query_name
        mapped_read_id_list.append(read_id)

    return mapped_read_id_list


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


os.chdir('/Users/songweizhi/Desktop/888')

sam_file = '9D42_sorted.bam'
pwd_blastn_exe = 'blastn'
refseq_file = '2.10wt_illumina.fasta'
refseq_id = '2.10_chromosome'
soft_reads_len_cutoff = 40
M_region_len = 200
region_span = '866233-941242'


# get dict
region_edge_dict = {}
region_span_split = region_span.split('-')
region_edge_dict[int(region_span_split[0])] = 'left'
region_edge_dict[int(region_span_split[1])] = 'right'

bridging_reads_full_len_db_file = 'bridging_reads_full_len_db.fasta'
bridging_reads_full_len_db_file_handle = open(bridging_reads_full_len_db_file, 'w')
bridging_reads_list = []
for each_edge in region_edge_dict:

    extracted_reads_file_soft_seg = 'extracted_reads_%s_soft_seg.fasta' % region_edge_dict[each_edge]

    # get M region
    start_pos = 0
    end_pos = 0
    if region_edge_dict[each_edge] == 'left':
        start_pos = each_edge
        end_pos = each_edge + M_region_len

    if region_edge_dict[each_edge] == 'right':
        start_pos = each_edge - M_region_len
        end_pos = each_edge

    samfile = pysam.AlignmentFile(sam_file, "rb")
    iter = samfile.fetch(refseq_id, start_pos, end_pos)

    extracted_reads_file_soft_seg_handle = open(extracted_reads_file_soft_seg, 'w')
    wrote_list = []
    mapped_read_id_list = []
    for each_aln in iter:
        each_aln = str(each_aln)
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

            # if the percentage of matched and soft clipping sequences higher than 10% and their total length higher than 80%,
            # export the matched and soft clipping sequences separately.
            if (longest_M / reads_len >= 0.2) and (longest_S / reads_len >= 0.2) and ((longest_M + longest_S) / reads_len >= 0.9):

                seq_to_keep = ''
                cigar_seg_to_keep = ''
                if region_edge_dict[each_edge] == 'left':
                    cigar_seg_to_keep = cigar_splitted[0]
                    seq_to_keep = get_cigar_seg_seq(cigar, cigar_splitted[0], reads_seq)

                if region_edge_dict[each_edge] == 'right':
                    cigar_seg_to_keep = cigar_splitted[-1]
                    seq_to_keep = get_cigar_seg_seq(cigar, cigar_splitted[-1], reads_seq)

                # ignore reads shorter than defined cutoff
                if len(seq_to_keep) >= soft_reads_len_cutoff:

                    # write out full length
                    reads_id_for_write = ''
                    if reads_id not in wrote_list:
                        reads_id_for_write = '%s_1' % reads_id
                        bridging_reads_full_len_db_file_handle.write('>%s\n' % reads_id_for_write)
                        bridging_reads_full_len_db_file_handle.write('%s\n' % reads_seq)
                        wrote_list.append(reads_id)
                    else:
                        reads_id_for_write = '%s_2' % reads_id
                        bridging_reads_full_len_db_file_handle.write('>%s\n' % reads_id_for_write)
                        bridging_reads_full_len_db_file_handle.write('%s\n' % reads_seq)

                    # write out segment
                    seq_to_keep_id = '%s_%s_%s' % (reads_id_for_write, region_edge_dict[each_edge], cigar_seg_to_keep)
                    extracted_reads_file_soft_seg_handle.write('>%s\n' % seq_to_keep_id)
                    extracted_reads_file_soft_seg_handle.write('%s\n' % seq_to_keep)

    extracted_reads_file_soft_seg_handle.close()


    ################################### blast extracted reads against reference sequences ##################################

    blast_output = 'blast_output_%s.tab' % region_edge_dict[each_edge]
    blast_parameters = '-evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn'
    cmd_blast = '%s -query %s -subject %s -out %s %s' % (pwd_blastn_exe, extracted_reads_file_soft_seg, refseq_file, blast_output, blast_parameters)
    os.system(cmd_blast)


    ################################################# parse blast results ##################################################

    for match in open(blast_output):
        match_split = match.strip().split('\t')
        query = match_split[0]
        query_no_suffixes = '_'.join(query.split('_')[:2])
        refseq = match_split[1]
        identity = float(match_split[2])
        align_len = int(match_split[3])
        query_len = int(match_split[12])
        subject_len = int(match_split[13])
        q_start = int(match_split[6])
        q_end = int(match_split[7])
        s_start = int(match_split[8])
        s_end = int(match_split[9])

        if (identity >= 99) and (align_len / query_len >= 0.9):

            if region_edge_dict[each_edge] == 'left':
                if s_end == 941242:
                    bridging_reads_list.append(query_no_suffixes)

            if region_edge_dict[each_edge] == 'right':
                if s_start == 866233:
                    bridging_reads_list.append(query_no_suffixes)

    os.system('rm %s' % blast_output)
    os.system('rm %s' % extracted_reads_file_soft_seg)

bridging_reads_full_len_db_file_handle.close()


# extract bridgingreads
bridging_reads_file = 'bridging_reads_for_assembly.fasta'
bridging_reads_file_handle = open(bridging_reads_file, 'w')

bridging_reads_list_uniq = unique_list_elements(bridging_reads_list)

# export
for each_seq in SeqIO.parse(bridging_reads_full_len_db_file, 'fasta'):

    if each_seq.id in bridging_reads_list_uniq:
        bridging_reads_file_handle.write('>%s\n' % each_seq.id)
        bridging_reads_file_handle.write('%s\n' % str(each_seq.seq))

bridging_reads_file_handle.close()


os.system('rm %s' % bridging_reads_full_len_db_file)

