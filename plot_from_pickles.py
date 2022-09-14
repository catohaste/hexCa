import pickle
import matplotlib.pyplot as plt
import random

from lib import *
from functions import *
from plot import animate_var_by_color, animate_var_over_x_avg_y, plot_var_over_time_fixed_x_avg_y


##################################################################################################
# results folder
results_dir = "results/"

save_dir = results_dir + '2022-04-27_1704 random ICs const V_PLC D0/'

##################################################################################################
# LOAD PICKLES

var_strings = ['Ca_cyt', 'IP3']
color_strings = ['Blues', 'Oranges']

pickle_dir = save_dir + "pickles/"

with open(pickle_dir + 'Ca_cyt' + '.pickle', 'rb') as handle:
    Ca_cyt_val_loc = pickle.load(handle)
with open(pickle_dir + 'ip3' + '.pickle', 'rb') as handle:
    ip3_val_loc = pickle.load(handle)
with open(pickle_dir + 'layout_dict.pickle', 'rb') as handle:
    layout_dict = pickle.load(handle)
pointy = create_layout_from_dict(layout_dict)

# loading hex_tuples not working, weirdly don't need it
# with open(pickle_dir + 'hex_tuples.pickle', 'wb') as handle:
#     hex_tuples = pickle.load(handle)
# hex_array = [tuple_to_hex(hexa) for hexa in hex_tuples]

hex_array_a, Ca_cyt = unpack_val_loc_tuple_std_layout(Ca_cyt_val_loc, pointy)
hex_array, ip3 = unpack_val_loc_tuple_std_layout(ip3_val_loc, pointy)

variables = [Ca_cyt, ip3]

store_timepoint_N = 2400
store_timepoint_N = len(Ca_cyt[hex_array[0]])

hex_x_N,hex_y_N = 40, 6

##################################################################################################
# PLOT

# plt.figure(figsize=(6, 6))
# plt.subplot(211)
# plt.plot(ip3[hex_array[0]])
# plt.subplot(212)
# plt.plot(ip3[hex_array[39]])
#
# plt.show()

for var, var_str, color_str in zip(variables, var_strings, color_strings):
    animate_var_by_color(var, store_timepoint_N, hex_array, (hex_x_N,hex_y_N), pointy, 12, color_str, save_dir + var_str)
    # animate_var_over_x_avg_y(var, store_timepoint_N, hex_array, (hex_x_N,hex_y_N), pointy, 12, color_str, var_str, save_dir + var_str)
    # plot_var_over_time_fixed_x_avg_y(var, hex_array, pointy, 12, color_str, var_str, save_dir + var_str)
    