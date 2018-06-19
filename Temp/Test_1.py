import os
import glob
from Bio import SeqIO


os.chdir('/Users/songweizhi/Desktop/refseq_genomes/archaea_253')
files = '/Users/songweizhi/Desktop/refseq_genomes/archaea_253/*.fna'

file_list = [os.path.basename(file_name) for file_name in glob.glob(files)]
print(file_list)

for each_genome in file_list:
    name_list = []

    length_list = []
    length_list_sum = 0
    for each_seq in SeqIO.parse(each_genome, 'fasta'):
        name_list.append(each_seq.description)
        length_list.append(str(len(str(each_seq.seq))))
        length_list_sum += len(str(each_seq.seq))
    print(name_list)
    #print('%s\t%s' % (each_genome, '\t'.join(length_list)))
    length_list_percent = []
    for each_len in length_list:
        each_len_percent = float("{0:.2f}".format(int(each_len)/length_list_sum))
        length_list_percent.append(str(each_len_percent))
    print('%s\t%s' % ('', '\t'.join(length_list_percent)))








