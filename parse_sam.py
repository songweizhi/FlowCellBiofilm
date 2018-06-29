# get the number of reads with mutation for identified SNV
import os
import pysam


os.chdir('/Users/songweizhi/Desktop/pysam_wd')
input_bam_file = "D2_ref_sorted.bam"
snv_pos = 24


samfile = pysam.AlignmentFile(input_bam_file, "rb")
for pileupcolumn in samfile.pileup("D2_c"):
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

        print('%s\t%s\t%s' % (pileupcolumn.pos + 1, pileupcolumn.n, '\t'.join(mapped_ncs)))

samfile.close()

