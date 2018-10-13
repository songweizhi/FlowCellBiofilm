
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


sequence = 'AAAAATTTTTGGGGGCCCCC'

cigar = '2M1I2D7S8M'
cigar_seg = '8M'



print(get_cigar_seg_seq(cigar, cigar_seg, sequence))