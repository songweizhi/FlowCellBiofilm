import numpy as np
from sklearn.cluster import KMeans
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from sklearn.datasets import make_blobs
import math
import random
import seaborn as sns


def get_color_list(color_num):

    if color_num <= 8:
        color_list_combined = ['#3787c0', '#39399f', '#ffb939', '#399f39', '#9f399f', '#fb694a', '#9f9f39', '#959595']

    elif 8 < color_num <= 16:
        color_list_combined = ['#2b7bba', '#89bedc', '#2e2e99', '#8a8acc', '#ffa500', '#ffc55c', '#2e992e', '#8acc8a', '#992e99', '#cc8acc', '#d52221', '#fc8161', '#99992e', '#cccc8a', '#5c5c5c', '#adadad']

    else:
        color_num_each = math.ceil(color_num/8) + 2

        color_list_1 = sns.color_palette('Blues',  n_colors=color_num_each).as_hex()
        color_list_2 = sns.light_palette('navy',   n_colors=color_num_each).as_hex()
        color_list_3 = sns.light_palette('orange', n_colors=color_num_each).as_hex()
        color_list_4 = sns.light_palette('green',  n_colors=color_num_each).as_hex()
        color_list_5 = sns.light_palette('purple', n_colors=color_num_each).as_hex()
        color_list_6 = sns.color_palette('Reds',   n_colors=color_num_each).as_hex()
        color_list_7 = sns.light_palette('olive',  n_colors=color_num_each).as_hex()
        color_list_8 = sns.color_palette('Greys', n_colors=color_num_each).as_hex()

        color_list_combined = []
        for color_list in [color_list_1, color_list_2, color_list_3, color_list_4, color_list_5, color_list_6, color_list_7, color_list_8]:
            for color in color_list[2:][::-1]:
                color_list_combined.append(color)

    color_list_to_return = random.sample(color_list_combined, color_num)

    color_list_to_return_sorted = []
    for color_to_return in color_list_combined:
        if color_to_return in color_list_to_return:
            color_list_to_return_sorted.append(color_to_return)

    return color_list_to_return_sorted


# create dataset
dots_coordinate_list = [[23, 14], [25, 15], [26, 16], [26, 17], [23, 17], [25, 16], [26, 14], [26, 15], [39, 51], [40, 52], [41, 51], [40, 49], [39, 50], [40, 52], [41, 49], [42, 51], [54, 34], [55, 35], [56, 36], [57, 37], [54, 35], [55, 37], [56, 34], [57, 36], [74, 54], [75, 55], [76, 56], [74, 55], [75, 54], [76, 56], [39, 81], [40, 82], [41, 81], [40, 79], [39, 80], [40, 82], [41, 79], [42, 81]]
coordinate_array = np.array(dots_coordinate_list)

# plot
# plt.scatter(coordinate_array[:, 0], coordinate_array[:, 1], c='white', marker='o', edgecolor='black', s=12)
# plt.show()

km = KMeans(n_clusters=5, init='random', n_init=10, max_iter=300, tol=1e-04, random_state=0)
y_km = km.fit_predict(coordinate_array)


print('km')
print(km)

print('y_km')
print(y_km)

print('cluster_centers_')
print(km.cluster_centers_)

# get uniq_cluster_id_list
uniq_cluster_id_list = []
for i in y_km:
    if i not in uniq_cluster_id_list:
        uniq_cluster_id_list.append(i)

# get color list
color_list = get_color_list(len(uniq_cluster_id_list))

# plot clusters
for cluster_id in sorted(uniq_cluster_id_list):
    plt.scatter(coordinate_array[y_km == cluster_id, 0], coordinate_array[y_km == cluster_id, 1], s=12, c=color_list[cluster_id], marker='o', edgecolor='', label='cluster %s' % (cluster_id + 1))

# plot centroids
plt.scatter(km.cluster_centers_[:, 0], km.cluster_centers_[:, 1], s=20, marker='*', c='red', edgecolor='red', label='centroids')

plt.xlim(0, 100)
plt.ylim(0, 100)

plt.legend(scatterpoints=1)
plt.tight_layout()
plt.savefig('/Users/songweizhi/Desktop/demo_cluster.pdf')
plt.close()
