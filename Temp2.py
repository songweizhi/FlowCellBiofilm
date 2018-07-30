import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord


def export_dna_record(gene_seq, gene_id, gene_description, output_handle):
    seq_object = Seq(gene_seq, IUPAC.unambiguous_dna)
    seq_record = SeqRecord(seq_object)
    seq_record.id = gene_id
    seq_record.description = gene_description
    SeqIO.write(seq_record, output_handle, 'fasta')

wd = '/Users/songweizhi/Desktop'
seq_file = 'combined_references.fasta'
seq_id = '2.10_chromosome'
start_pos = 866221
end_pos = 941272

os.chdir(wd)


output_file = '%s_%s_%s.fasta' % (seq_id, start_pos, end_pos)
output_file_handle = open(output_file, 'w')
for each_seq in SeqIO.parse(seq_file, 'fasta'):

    if each_seq.id == seq_id:
        gene_id = '%s_%s_%s' % (seq_id, start_pos, end_pos)
        need_seq = str(each_seq.seq)[start_pos-1: end_pos]
        export_dna_record(need_seq, gene_id, '', output_file_handle)

output_file_handle.close()








