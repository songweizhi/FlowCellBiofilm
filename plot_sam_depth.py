
import os
import argparse
from Bio import SeqIO
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def take_kmer_mean(num_list, k_mer):

    k_mer_average_list = []
    n = 0

    while n < len(num_list):

        if n <= len(num_list) - k_mer:

            to_add_value = list(range(1, k_mer + 1))

            pos_list = []
            for each in to_add_value:
                pos_list.append(num_list[n + each - 1])

            # get average
            current_average = sum(pos_list) / k_mer
            k_mer_average_list.append(float("{0:.2f}".format(current_average)))

            n += 1

        elif (len(num_list) - k_mer) < n < len(num_list):
            k_mer_average_list.append(num_list[n])

            n += 1

    return k_mer_average_list


def main(depth_file, seq_to_plot, start_pos, end_pos, plot_filename):
    #print('Extracting absolute depth from input file')
    x = []
    y = []
    bp_num = 1
    current_pos = 0
    for each_base in open(depth_file):
        each_base_split = each_base.strip().split('\t')
        seq_id = each_base_split[0]
        pos = int(each_base_split[1])
        depth = int(each_base_split[2])

        if seq_id == seq_to_plot:

            if pos < start_pos:
                pass

            elif pos == start_pos:
                x.append(pos)
                y.append(depth)
                current_pos = pos

            else:
                start_0 = None
                end_0 = None
                to_add = []

                if (pos == current_pos + 1) and (pos <= end_pos):
                    x.append(pos)
                    y.append(depth)
                    current_pos = pos

                elif (pos > current_pos + 1) and (pos <= end_pos):

                    # add zero
                    start_0 = current_pos + 1
                    end_0 = pos - 1
                    to_add = list(range(start_0, end_0 + 1))
                    for each_0 in to_add:
                        x.append(each_0)
                        y.append(0)

                    x.append(pos)
                    y.append(depth)
                    current_pos = pos

                elif (pos > current_pos + 1) and (pos > end_pos):

                    # add zero
                    start_0 = current_pos + 1
                    end_0 = end_pos
                    to_add = list(range(start_0, end_0 + 1))
                    for each_0 in to_add:
                        x.append(each_0)
                        y.append(0)

                    current_pos = pos

    #print('Calculating k-mer means')
    y = take_kmer_mean(y, k_mer)

    #print('Plotting')
    # Change the color and its transparency
    plt.plot(x, y, color="skyblue", alpha=0.7, linewidth=0.7)

    # titles
    plt.title(plot_filename)
    plt.xlabel('Position (bp)', fontsize=10)
    plt.ylabel('Depth (X)', fontsize=10)

    # Get plot
    plt.savefig('%s.png' % plot_filename, dpi=300)
    plt.close()


################################################# input #################################################

parser = argparse.ArgumentParser()

parser.add_argument('-r', required=True, help='reference sequence file')
parser.add_argument('-d', required=True, help='depth file')
parser.add_argument('-i', required=False, default=None, help='id of sequence to plot')
parser.add_argument('-s', required=False, type=int, default=None, help='start position to plot')
parser.add_argument('-e', required=False, type=int, default=None, help='end position to plot')
parser.add_argument('-k', required=False, type=int, default=100, help='k-mer mean depth')

args = vars(parser.parse_args())
sequence_file = args['r']
depth_file = args['d']
seq_to_plot = args['i']
start_pos = args['s']
end_pos = args['e']
k_mer = args['k']


#########################################################################################################

# get sequence length dict
seq_id_length_dict = {}
for each_seq in SeqIO.parse(sequence_file, 'fasta'):
    seq_id_length_dict[each_seq.id] = len(each_seq.seq)


if seq_to_plot not in seq_id_length_dict:
    print('Reference sequence %s not found in bam file, program exited!' % seq_to_plot)
    exit()


# get depth_file base name
depth_file_basename, depth_file_extension = os.path.splitext(depth_file)


if seq_to_plot != None:

    # get the start and end of the region to plot
    if start_pos == None:
        start_pos = 1
    if end_pos == None:
        end_pos = seq_id_length_dict[seq_to_plot]

    print('Processing %s' % seq_to_plot)
    plot_filename = '%s__%s__%s-%sbp__%smer' % (depth_file_basename, seq_to_plot, start_pos, end_pos, k_mer)
    main(depth_file, seq_to_plot, start_pos, end_pos, plot_filename)

if seq_to_plot == None:

    for each_ctg in SeqIO.parse(sequence_file, 'fasta'):
        print('Processing %s' % each_ctg.id)
        plot_filename = '%s__%s__%s-%sbp__%smer' % (depth_file_basename, each_ctg.id, 1, seq_id_length_dict[each_ctg.id], k_mer)
        main(depth_file, each_ctg.id, 1, seq_id_length_dict[each_ctg.id], plot_filename)
