
sam_file        = '4D42.sam'
sam_file_cigar  = '4D42_cigar.txt'


sam_file_cigar_handle = open(sam_file_cigar, 'w')
for each_read in open(sam_file):

    if not each_read.startswith('@'):
        each_read_split = each_read.strip().split('\t')
        read_id         = each_read_split[0]
        read_id_base    = '_'.join(read_id.split('_')[:-1])
        read_strand     = read_id.split('_')[-1]
        ref_id          = each_read_split[2]
        ref_pos         = int(each_read_split[3])
        cigar           = each_read_split[5]
        read_seq        = each_read_split[9]
        read_len        = len(read_seq)

        sam_file_cigar_handle.write('%s\t%s\n' % (cigar, read_len))

sam_file_cigar_handle.close()
