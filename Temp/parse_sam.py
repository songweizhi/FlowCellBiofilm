# get the number of reads with mutation for identified SNV
import os
import pysam


os.chdir('/Users/songweizhi/Desktop/pysam_wd')
input_bam_file = "D2_ref_sorted.bam"
input_bam_file = "/Users/songweizhi/Dropbox/Research/Flow_cell/deepSNV_output_summary/8D27.bam"
#detected_snv = 'D2_c|24|C|G'
detected_snv = 'D2_c|107077|T|-'


def get_depth_and_snv_reads_num(input_bam_file, detected_snv):

    detected_snv_split = detected_snv.split('|')
    refseq_id = detected_snv_split[0]
    snv_pos = int(detected_snv_split[1])
    snv_pos_wt = detected_snv_split[2]
    snv_pos_v = detected_snv_split[3]

    samfile = pysam.AlignmentFile(input_bam_file, "rb")

    current_pos_depth = 0
    snv_pos_v_num = 0

    for pileupcolumn in samfile.pileup(refseq_id):

        if pileupcolumn.pos + 1 == snv_pos:
            mapped_ncs = []
            for pileupread in pileupcolumn.pileups:
                current_bp = ''
                if pileupread.is_del:
                    current_bp = '-'
                elif pileupread.is_refskip:
                    current_bp = 'S'
                else:
                    current_bp = pileupread.alignment.query_sequence[pileupread.query_position]
                mapped_ncs.append(current_bp)

            current_pos_depth = pileupcolumn.n
            snv_pos_v_num = mapped_ncs.count(snv_pos_v)

    samfile.close()

    return current_pos_depth, snv_pos_v_num


snv_depth, snv_reads_count = get_depth_and_snv_reads_num(input_bam_file, detected_snv)


print(snv_depth)
print(snv_reads_count)

