import matplotlib.pyplot as plt

import numpy as np
from scipy import stats

#from sklearn import datasets, linear_model
# from sklearn.metrics import mean_squared_error, r2_score

from scipy.stats import linregress


def get_affect(frequency_list, plot_filename):
    time_point = [9, 18, 27, 42]
    time_point_rescaled = []
    for each_tp in time_point:
        each_tp_rescaled = each_tp / 42
        time_point_rescaled.append(each_tp_rescaled)
    time_point_rescaled_arrary = np.array(time_point_rescaled)

    frequency_arrary = np.array(frequency_list)
    min_frequency = frequency_arrary.min()
    max_frequency = frequency_arrary.max()
    slope, intercept, r_value, p_value, std_err = stats.linregress(time_point_rescaled_arrary, frequency_arrary)
    print(slope)

    if (slope > 0) and (max_frequency >= 10):
        affect = 'Beneficial'
    elif (slope < 0) and (max_frequency <= 10):
        affect = 'Harmful'
    else:
        affect = 'Neutral'

    # plot
    plt.plot(time_point_rescaled_arrary, frequency_arrary, 'o')
    plt.plot(time_point_rescaled_arrary, intercept + slope*time_point_rescaled_arrary, 'r')

    # add text
    x_min = plt.xlim()[0]  # get the x-axes minimum value
    x_max = plt.xlim()[1]  # get the x-axes maximum value
    y_min = plt.ylim()[0]  # get the y-axes minimum value
    y_max = plt.ylim()[1]  # get the y-axes maximum value

    # set text position
    text_x = 0.2
    text_y_slope = y_min + (y_max - y_min) / 5 * 4.4
    text_y_p_value = y_min + (y_max - y_min) / 5 * 4.1
    text_y_affect = y_min + (y_max - y_min) / 5 * 3.8
    plt.text(text_x, text_y_slope, 'Slope: %s' % float("{0:.2f}".format(slope)))
    plt.text(text_x, text_y_p_value, 'P_value: %s' % float("{0:.2f}".format(p_value)))
    plt.text(text_x, text_y_affect, 'Affect: %s' % affect)

    plt.xlabel('Frequency')
    plt.ylabel('Time point')
    plt.savefig('%s.png' % plot_filename, dpi=300)
    plt.close()

    return affect





frequency_list = [1, 3, 2, 12]

affect = get_affect(frequency_list)

print(affect)



# time_point_rescaled = []
# for each_tp in time_point:
#     each_tp_rescaled = each_tp/42
#     time_point_rescaled.append(each_tp_rescaled)
#
#
# frequency_arrary = np.array(frequency_list)
# time_point_rescaled_arrary = np.array(time_point_rescaled)
# slope, intercept, r_value, p_value, std_err = stats.linregress(time_point_rescaled_arrary, frequency_arrary)
#
#
# print('slope: %s' % slope)
# print('intercept: %s' % intercept)
# print('r_value: %s' % r_value)
# print('p_value: %s' % p_value)
# print('std_err: %s' % std_err)
#
#
# # plot
# plt.plot(time_point_rescaled_arrary, frequency_arrary, 'o')
# plt.plot(time_point_rescaled_arrary, intercept + slope*time_point_rescaled_arrary, 'r')
#
# add text
# x_min = plt.xlim()[0]  # get the x-axes minimum value
# x_max = plt.xlim()[1]  # get the x-axes maximum value
# y_min = plt.ylim()[0]  # get the y-axes minimum value
# y_max = plt.ylim()[1]  # get the y-axes maximum value
#
# # set text position
# text_x = 0.2
# text_y_slope = y_min + (y_max - y_min) / 5 * 4.4
# text_y_p_value = y_min + (y_max - y_min) / 5 * 4.1
#
# plt.text(text_x, text_y_slope, 'Slope: %s' % slope)
# plt.text(text_x, text_y_p_value, 'P_value: %s' % p_value)
#
#
# plt.legend()
# plt.xlabel('Frequency')
# plt.ylabel('Time point')
# plt.show()




