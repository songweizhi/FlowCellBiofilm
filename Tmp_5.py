
def cigar_splitter(cigar):

    # get the position of letters
    letter_pos_list = []
    n = 0
    for each_element in cigar:
        if each_element.isalpha() is True:
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


sam_file = '4D9.sam'

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
        cigar_splitted  = cigar_splitter(cigar)

        if ('S' in cigar) and (len(cigar_splitted) == 2):
            print(cigar)
