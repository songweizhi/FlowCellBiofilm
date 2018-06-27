# library
import os
import matplotlib.pyplot as plt


os.chdir('/Users/songweizhi/Desktop/')

depth_file = '9D27.cov'
seq_to_plot = '2.10_chromosome'
start_pos = 1
end_pos = 2000
title = '10_chromosome'

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

# Change the color and its transparency
plt.fill_between(x, y, color="skyblue", alpha=0.4)

# Get plot
plt.savefig('%s.png' % (title), dpi=300)
plt.close()
