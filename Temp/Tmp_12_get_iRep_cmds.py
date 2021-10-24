
sample_id_list = ['1D18', '1D27', '1D42', '1D9', '5D18', '5D27', '5D42', '5D9', '9D18', '9D27', '9D42', '9D9', '10D18', '10D27', '10D42', '10D9', '2D18', '2D27', '2D42', '2D9', '6D18', '6D27', '6D42', '6D9', '12D18', '12D27', '12D42', '12D9', '4D18', '4D27', '4D42', '4D9', '8D18', '8D27', '8D42', '8D9']

for sample_id in sample_id_list:

    file_out = '/Users/songweizhi/Desktop/js/js_%s_iRep.sh' % sample_id
    file_out_handle = open(file_out, 'w')
    file_out_handle.write('#!/bin/bash\n#PBS -l nodes=1:ppn=12\n#PBS -l mem=60gb\n#PBS -l walltime=11:59:00\n#PBS -M weizhi.song@unsw.edu.au\n#PBS -j oe\n#PBS -m ae\n\n')
    file_out_handle.write('module load python/3.7.3\nsource ~/mypython3env/bin/activate\nmodule load bowtie/2.3.5.1\n')
    file_out_handle.write('\n')

    if sample_id in ['1D18', '1D27', '1D42', '1D9', '5D18', '5D27', '5D42', '5D9', '9D18', '9D27', '9D42', '9D9']:
        file_out_handle.write('cd /srv/scratch/z5039045/Flow_cell_biofilm/7_iRep/quality_filtered\n')
        file_out_handle.write('gunzip %s_R1_Q30_P.fastq.gz\n' % sample_id)
        file_out_handle.write('gunzip %s_R2_Q30_P.fastq.gz\n' % sample_id)
        file_out_handle.write('\n')
        file_out_handle.write('cd /srv/scratch/z5039045/Flow_cell_biofilm/7_iRep\n')
        file_out_handle.write('bowtie2 -x refs/2.10wt_illumina_c -1 quality_filtered/%s_R1_Q30_P.fastq -2 quality_filtered/%s_R2_Q30_P.fastq -S %s_c.sam --reorder -p 12\n' % (sample_id, sample_id, sample_id))
        file_out_handle.write('cd /srv/scratch/z5039045/Flow_cell_biofilm/7_iRep/quality_filtered\n')
        file_out_handle.write('rm %s_R1_Q30_P.fastq\n' % sample_id)
        file_out_handle.write('rm %s_R2_Q30_P.fastq\n' % sample_id)
        file_out_handle.write('cd /srv/scratch/z5039045/Flow_cell_biofilm/7_iRep\n')
        file_out_handle.write('bPTR -f refs/2.10wt_illumina_c.fasta -s %s_c.sam -o %s_c.bPTR.tsv -plot %s_c.bPTR.pdf -m coverage\n' % (sample_id, sample_id, sample_id))
        file_out_handle.write('rm %s_c.sam\n' % sample_id)

    if sample_id in ['10D18', '10D27', '10D42', '10D9', '2D18', '2D27', '2D42', '2D9', '6D18', '6D27', '6D42', '6D9']:
        file_out_handle.write('cd /srv/scratch/z5039045/Flow_cell_biofilm/7_iRep/quality_filtered\n')
        file_out_handle.write('gunzip %s_R1_Q30_P.fastq.gz\n' % sample_id)
        file_out_handle.write('gunzip %s_R2_Q30_P.fastq.gz\n' % sample_id)
        file_out_handle.write('\n')
        file_out_handle.write('cd /srv/scratch/z5039045/Flow_cell_biofilm/7_iRep\n')
        file_out_handle.write('bowtie2 -x refs/D2_pacbio_c -1 quality_filtered/%s_R1_Q30_P.fastq -2 quality_filtered/%s_R2_Q30_P.fastq -S %s_c.sam --reorder -p 12\n' % (sample_id, sample_id, sample_id))
        file_out_handle.write('cd /srv/scratch/z5039045/Flow_cell_biofilm/7_iRep/quality_filtered\n')
        file_out_handle.write('rm %s_R1_Q30_P.fastq\n' % sample_id)
        file_out_handle.write('rm %s_R2_Q30_P.fastq\n' % sample_id)
        file_out_handle.write('cd /srv/scratch/z5039045/Flow_cell_biofilm/7_iRep\n')
        file_out_handle.write('bPTR -f refs/D2_pacbio_c.fasta -s %s_c.sam -o %s_c.bPTR.tsv -plot %s_c.bPTR.pdf -m coverage\n' % (sample_id, sample_id, sample_id))
        file_out_handle.write('rm %s_c.sam\n' % sample_id)

    if sample_id in ['12D18', '12D27', '12D42', '12D9', '4D18', '4D27', '4D42', '4D9', '8D18', '8D27', '8D42', '8D9']:
        file_out_handle.write('cd /srv/scratch/z5039045/Flow_cell_biofilm/7_iRep/quality_filtered\n')
        file_out_handle.write('gunzip %s_R1_Q30_P.fastq.gz\n' % sample_id)
        file_out_handle.write('gunzip %s_R2_Q30_P.fastq.gz\n' % sample_id)
        file_out_handle.write('\n')
        file_out_handle.write('cd /srv/scratch/z5039045/Flow_cell_biofilm/7_iRep\n')
        file_out_handle.write('bowtie2 -x refs/2.10wt_illumina_c -1 quality_filtered/%s_R1_Q30_P.fastq -2 quality_filtered/%s_R2_Q30_P.fastq -S %s_210_c.sam --reorder -p 12\n' % (sample_id, sample_id, sample_id))
        file_out_handle.write('bowtie2 -x refs/D2_pacbio_c -1 quality_filtered/%s_R1_Q30_P.fastq -2 quality_filtered/%s_R2_Q30_P.fastq -S %s_D2_c.sam --reorder -p 12\n' % (sample_id, sample_id, sample_id))
        file_out_handle.write('\n')
        file_out_handle.write('cd /srv/scratch/z5039045/Flow_cell_biofilm/7_iRep/quality_filtered\n')
        file_out_handle.write('rm %s_R1_Q30_P.fastq\n' % sample_id)
        file_out_handle.write('rm %s_R2_Q30_P.fastq\n' % sample_id)
        file_out_handle.write('\n')
        file_out_handle.write('cd /srv/scratch/z5039045/Flow_cell_biofilm/7_iRep\n')
        file_out_handle.write('bPTR -f refs/2.10wt_illumina_c.fasta -s %s_210_c.sam -o %s_210_c.bPTR.tsv -plot %s_210_c.bPTR.pdf -m coverage\n' % (sample_id, sample_id, sample_id))
        file_out_handle.write('bPTR -f refs/D2_pacbio_c.fasta -s %s_D2_c.sam -o %s_D2_c.bPTR.tsv -plot %s_D2_c.bPTR.pdf -m coverage\n' % (sample_id, sample_id, sample_id))
        file_out_handle.write('\n')
        file_out_handle.write('rm %s_210_c.sam\n' % sample_id)
        file_out_handle.write('rm %s_D2_c.sam\n' % sample_id)

    file_out_handle.close()

