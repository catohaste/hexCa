import time
import pickle
import matplotlib.pyplot as plt
import random

from lib import *
from functions import *
from plot import thesis_vary_constant_connections_animation

##################################################################################################
# results folder
results_dir = "results/vary_connections/"
results_dir = "results/vary_connections_short/"

folder_list = ['connect_frac_0_0', 'connect_frac_0_1', 'connect_frac_0_2', 'connect_frac_0_4', 'connect_frac_0_8', 'connect_frac_1_0']
# folder_list = ['connect_frac_0_8', 'connect_frac_1_0']

save_dir_list = [results_dir + x + '/' for x in folder_list]

variations = [0, 0.1, 0.2, 0.4, 0.8, 1.0]

hex_x_N,hex_y_N = 60, 6
hex_grid_dim = hex_x_N, hex_y_N
dimension_prefix_str = str(hex_x_N) + '_' + str(hex_y_N) + '_'

##################################################################################################
# LOAD PICKLES

pickle_dir_list = [save_dir + "pickles/" for save_dir in save_dir_list]

# sort out layout and hexes with dummy load
dummy_pickle_dir = pickle_dir_list[0]
with open(dummy_pickle_dir + 'layout_dict.pickle', 'rb') as handle:
    layout_dict = pickle.load(handle)
pointy = create_layout_from_dict(layout_dict)
with open(dummy_pickle_dir + 'Ca_cyt' + '.pickle', 'rb') as handle:
    dummy_Ca_cyt_val_loc = pickle.load(handle)
hex_array, dummy_Ca_cyt = unpack_val_loc_tuple_std_layout(dummy_Ca_cyt_val_loc, pointy)

varied_Ca_cyt_val_loc = {}
varied_ip3_val_loc = {}
varied_connect_params = {}
varied_connections = {}
for variation, pickle_dir in zip(variations, pickle_dir_list):
    with open(pickle_dir + 'Ca_cyt' + '.pickle', 'rb') as handle:
        varied_Ca_cyt_val_loc[variation] = pickle.load(handle)
    with open(pickle_dir + 'ip3' + '.pickle', 'rb') as handle:
        varied_ip3_val_loc[variation] = pickle.load(handle)
    with open(pickle_dir + 'connect_params.pickle', 'rb') as handle:
        varied_connect_params[variation] = pickle.load(handle)
        connect_str = str(varied_connect_params[variation]['init_avg_degree_fraction']) + '_'
        varied_connections[variation] = load_connections_from_file(hex_array, path.join(pickle_dir, dimension_prefix_str + connect_str))

# hex_array_a, Ca_cyt = unpack_val_loc_tuple_std_layout(Ca_cyt_val_loc, pointy)
# hex_array, ip3 = unpack_val_loc_tuple_std_layout(ip3_val_loc, pointy)

# variables = [Ca_cyt, ip3]

store_timepoint_N = 1200
store_timepoint_N = len(dummy_Ca_cyt[hex_array[0]])
store_dt = 0.5

##################################################################################################

var_strings = ['Ca_cyt', 'IP3']
color_strings = ['Oranges', 'Blues']
varied_var_loc_list = [varied_Ca_cyt_val_loc, varied_ip3_val_loc]

var_strings = ['IP3']
color_strings = ['Blues']
varied_var_loc_list = [varied_ip3_val_loc]

plot_time_indicies = [0]

start_time = time.time()

for varied_var_loc, var_str, color_str in zip(varied_var_loc_list, var_strings, color_strings):
    thesis_vary_constant_connections_animation(varied_var_loc, varied_connections, store_timepoint_N, store_dt, hex_array, hex_grid_dim, pointy, 8, var_str, color_str, results_dir + var_str + '_vary_connect_anim')
    
end_time = time.time()
print('Time animating', end_time - start_time)

