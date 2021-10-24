import argparse
from Bio import SeqIO

'''
module load python/3.7.3
source ~/mypython3env/bin/activate

cd /srv/scratch/z5039045/Flow_cell_biofilm/3_mapping_mapHGT
python3 get_mean_ref_depth.py -r ../0_References_for_mapHGT/combined_refs.fasta -s 4D9.sam &
python3 get_mean_ref_depth.py -r ../0_References_for_mapHGT/combined_refs.fasta -s 4D18.sam &
python3 get_mean_ref_depth.py -r ../0_References_for_mapHGT/combined_refs.fasta -s 4D27.sam &
python3 get_mean_ref_depth.py -r ../0_References_for_mapHGT/combined_refs.fasta -s 4D42.sam &

python3 get_mean_ref_depth.py -r ../0_References_for_mapHGT/combined_refs.fasta -s 8D9.sam &
python3 get_mean_ref_depth.py -r ../0_References_for_mapHGT/combined_refs.fasta -s 8D18.sam &
python3 get_mean_ref_depth.py -r ../0_References_for_mapHGT/combined_refs.fasta -s 8D27.sam &
python3 get_mean_ref_depth.py -r ../0_References_for_mapHGT/combined_refs.fasta -s 8D42.sam &

python3 get_mean_ref_depth.py -r ../0_References_for_mapHGT/combined_refs.fasta -s 12D9.sam &
python3 get_mean_ref_depth.py -r ../0_References_for_mapHGT/combined_refs.fasta -s 12D18.sam &
python3 get_mean_ref_depth.py -r ../0_References_for_mapHGT/combined_refs.fasta -s 12D27.sam &
python3 get_mean_ref_depth.py -r ../0_References_for_mapHGT/combined_refs.fasta -s 12D42.sam &


Sample	2.10	    D2
4D9	    732.28	    3016.83
4D18	715.28	    1395.56
4D27	1229.44	    863.58
4D42	1622.12	    971.38

8D9	    777.86	    1129.35
8D18	1136.92	    1282.73
8D27	2091.82	    686.11
8D42	1435.96	    662.88

12D9	1285.9	    667.66
12D18	447.98	    1539.01
12D27	1043.35	    891.18
12D42	2442.73	    13.66

'''


parser = argparse.ArgumentParser()

parser.add_argument('-r', required=True, help='reference')
parser.add_argument('-s', required=True, help='sam file')

args = vars(parser.parse_args())
ref_file     = args['r']
sam_file     = args['s']


ref_length_dict = {}
for seq_record in SeqIO.parse(ref_file, 'fasta'):
    ref_length_dict[seq_record.id] = len(seq_record.seq)

ref_to_mapped_reads_total_length = {}
for each_read in open(sam_file):

    if not each_read.startswith('@'):
        each_read_split = each_read.strip().split('\t')
        ref_id          = each_read_split[2]
        read_len        = len(each_read_split[9])

        if ref_id not in ref_to_mapped_reads_total_length:
            ref_to_mapped_reads_total_length[ref_id] = read_len
        else:
            ref_to_mapped_reads_total_length[ref_id] += read_len

print('Timepoint\tRef\tDepth')
for ref in ref_length_dict:
    ref_len = ref_length_dict[ref]
    ref_mapped_reads_len = ref_to_mapped_reads_total_length[ref]
    ref_seq_depth = float("{0:.2f}".format(ref_mapped_reads_len/ref_len))
    print('%s\t%s\t%s' % (sam_file, ref, ref_seq_depth))
