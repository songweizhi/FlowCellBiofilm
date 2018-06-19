import os
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna


def export_dna_record(gene_seq, gene_id, gene_description, output_handle):
    seq_object = Seq(gene_seq, IUPAC.unambiguous_dna)
    seq_record = SeqRecord(seq_object)
    seq_record.id = gene_id
    seq_record.description = gene_description
    SeqIO.write(seq_record, output_handle, 'fasta')


def simulate_reads(pwd_ref_file, read_number, read_length, insert_size, pwd_output_folder):

    # get the id and length list of reference sequences
    seq_id_list = []
    seq_length_list = []
    for each_seq in SeqIO.parse(pwd_ref_file, 'fasta'):
        seq_id_list.append(each_seq.id)
        seq_length_list.append(len(each_seq.seq))

    # get the number of reads need to be simulated from each sequence
    read_number_list = []
    n = 1
    for each_length in seq_length_list:
        if n < len(seq_length_list):
            current_reads_number = round(read_number*each_length/sum(seq_length_list))
            n += 1
            read_number_list.append(current_reads_number)
        elif n == len(seq_length_list):
            current_reads_number = read_number - sum(read_number_list)
            read_number_list.append(current_reads_number)

    # prepare output file name
    path, file_name = os.path.split(pwd_ref_file)
    genome_name, ext = os.path.splitext(file_name)
    output_r1 = '%s/%s_R1.fasta' % (pwd_output_folder, genome_name)
    output_r2 = '%s/%s_R2.fasta' % (pwd_output_folder, genome_name)
    output_r12 = '%s/%s_R12.fasta' % (pwd_output_folder, genome_name)

    # create output reads file
    output_r1_handle = open(output_r1, 'w')
    output_r2_handle = open(output_r2, 'w')

    # simulate reads
    fragment_length = 2 * read_length + insert_size
    m = 0
    total_simulated = 1
    for each_seq2 in SeqIO.parse(pwd_ref_file, 'fasta'):
        reads_num_to_simulate = read_number_list[m]
        current_seq = str(each_seq2.seq)
        current_seq_len = seq_length_list[m]
        simulated = 1
        while simulated <= reads_num_to_simulate:
            random_position = random.randint(1, current_seq_len) # random_position
            current_fragment = ''

            # if simulated reads located at the middle
            if (random_position + fragment_length) <= current_seq_len:
                current_fragment = current_seq[random_position - 1: random_position + fragment_length - 1]

            # get forward and reverse reads
            current_fragment_r1 = current_fragment[:read_length]
            current_fragment_r2 = current_fragment[-read_length:]
            current_fragment_r2_reverse_complement = str(Seq(current_fragment_r2, generic_dna).reverse_complement())
            current_read_r1_id = '%s_||_%s_r1' % (genome_name, total_simulated)
            current_read_r2_id = '%s_||_%s_r2' % (genome_name, total_simulated)

            # write out sequence
            if current_fragment != '':
                export_dna_record(current_fragment_r1, current_read_r1_id, '', output_r1_handle)
                export_dna_record(current_fragment_r2_reverse_complement, current_read_r2_id, '', output_r2_handle)
                simulated += 1
                total_simulated += 1
        m += 1

    # close file
    output_r1_handle.close()
    output_r2_handle.close()


wd = '/Users/songweizhi/Desktop/MapHGT/ref_with_1HGT'
os.chdir(wd)
seqin_file = '/Users/songweizhi/Desktop/MapHGT/ref_with_1HGT/combined_ref_with_HGT.fasta'


simulate_reads(seqin_file, 500000, 250, 500, '.')
