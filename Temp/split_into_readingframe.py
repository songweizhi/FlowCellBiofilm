gene_len = 2280


gene_seq = 'ATGAAGATCAGCATCGAACGCGGCACCCTGTTGAAAGCTGTGGCTCAGGCCCAGTCAGTGAAACCCGGG'


n = 0
rf_list = []
rf_pos_list = []
while n < len(gene_seq)/3:
    start_pos = 1 + 3 * n
    end_pos = 3 * (n + 1)
    current_rf = gene_seq[(start_pos - 1): end_pos]
    rf_list.append(current_rf)
    rf_pos_list.append([start_pos, start_pos + 1, start_pos + 2])
    n += 1

print(rf_list)
print(rf_pos_list)




def get_rf_pos_list(gene_len):

    rf_pos_list = []
    n = 0
    while n < (gene_len / 3):
        start_pos = 1 + 3 * n
        rf_pos_list.append([start_pos, start_pos + 1, start_pos + 2])
        n += 1

    return rf_pos_list


print(get_rf_pos_list(180))




