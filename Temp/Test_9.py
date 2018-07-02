
from Bio import SeqIO

ref_length_dict = {}
for each in SeqIO.parse('/Users/songweizhi/Desktop/combined_references.fasta', 'fasta'):
    ref_length_dict[each.id] = len(each.seq)

print(ref_length_dict)

ref_length_dict = {'2.10_chromosome': 3758219, '2.10_plasmid1': 237747, '2.10_plasmid2': 94490, '2.10_plasmid3': 70353, 'D2_c': 4010148, 'D2_p': 1010395}