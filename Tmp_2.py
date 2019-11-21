



# check get_read_common_part function
def get_read_common_part(read_id):

    read_id_common_part = '_'.join(read_id.split('_')[:-1])

    return read_id_common_part


sam_file = '/Users/songweizhi/Desktop/simulate_reads/combined_ref.sam'


read_common_part_to_ref_dict = {}
for each in open(sam_file):
    if not each.startswith('@'):
        each_split = each.strip().split('\t')
        read_id = each_split[0]

        # check read common part
        read_id_common_part = get_read_common_part(read_id)
        ref_name = each_split[2]
        query_seq = each_split[9]

        if read_id_common_part not in read_common_part_to_ref_dict:
            read_common_part_to_ref_dict[read_id_common_part] = {ref_name}
        else:
            read_common_part_to_ref_dict[read_id_common_part].add(ref_name)


multi_ref_matched_reads = set()
for read_to_ref in read_common_part_to_ref_dict:
    ref_set = read_common_part_to_ref_dict[read_to_ref]
    if len(ref_set) > 1:
        multi_ref_matched_reads.add(read_to_ref)


print(multi_ref_matched_reads)



for each in open(sam_file):
    if not each.startswith('@'):
        each_split = each.strip().split('\t')
        read_id = each_split[0]
        read_id_common_part = get_read_common_part(read_id)


        if read_id_common_part in multi_ref_matched_reads:

            print(each_split)







