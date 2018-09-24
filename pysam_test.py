import pysam


def get_reads_mapped_to_a_region(sam_file, seq_id, start_pos, end_pos):

    samfile = pysam.AlignmentFile(sam_file, "rb")
    iter = samfile.fetch(seq_id, start_pos, end_pos)

    mapped_read_id_list = []
    for mapped_read in iter:
        read_id = mapped_read.query_name
        mapped_read_id_list.append(read_id)

    return mapped_read_id_list


def get_mates_of_reads_mapped_to_a_region(sam_file, seq_id, start_pos, end_pos):

    samfile = pysam.AlignmentFile(sam_file, "rb")
    iter = samfile.fetch(seq_id, start_pos, end_pos)

    for mapped_read in iter:

        if (mapped_read.is_paired) and (not mapped_read.mate_is_unmapped) and (not mapped_read.is_duplicate):
            print(mapped_read)
            print(samfile.mate(mapped_read))


        # if mapped_read.is_paired and not mapped_read.mate_is_unmapped and not mapped_read.is_duplicate:
        #
        #     mapped_read_mate = samfile.mate(mapped_read)
        #
        #     print(mapped_read_mate)



sam_file = '/Users/songweizhi/Desktop/pysam_test/ref_seq.bam'

mapped_read_list = get_reads_mapped_to_a_region(sam_file, 'ref_seq', 278, 281)
print(mapped_read_list)


get_mates_of_reads_mapped_to_a_region(sam_file, 'ref_seq', 278, 281)






