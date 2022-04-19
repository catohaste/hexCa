from os import mkdir, path
import random
import pickle
import datetime
from shutil import copy2
import time

import matplotlib.pyplot as plt

from lib import *
from functions import *
from models import *
from plot import plot_hexes, plot_var_by_color, animate_var_by_color

from params import params

##################################################################################################
# set up results folder
results_dir = "results/"
if not path.isdir(results_dir):
    mkdir(results_dir)

now = datetime.datetime.now()
now_str = now.strftime("%Y-%m-%d_%H%M/")
save_dir = results_dir + now_str
save_dir = results_dir + 'dev/'
if not path.isdir(save_dir):
    mkdir(save_dir)
    
# save code
code_dir = save_dir + "code/"
if not path.isdir(code_dir):
    mkdir(code_dir)
filenames = ['main.py', 'functions.py', 'plot.py']
for filename in filenames:
    copy2(filename, code_dir + filename)

################################################################################################## 
# set up layout
layout_dict = {
    "layout_str": 'pointy',
    "center_x": 0,
    "center_y": 0,
    "radius": 10.0
}
pointy = create_layout_from_dict(layout_dict)

##################################################################################################
# set up hexagonal grid with (q,r,s) coordinates
hex_x_N = 10
hex_y_N = 3
hex_array = []
# PARALLELOGRAM MAP
# for x in range(hex_x_N):
#     for y in range(hex_y_N):
#         hex_array.append(Hex(-x-y,y,x))

# RECTANGLE MAP
for y in range(hex_y_N):
    y_offset = int(np.floor(y/2))
    for x in range(-y_offset, hex_x_N - y_offset):
        hex_array.append(Hex(-x-y,y,x))
        
##################################################################################################

t_endpoint = 100

dt = 0.001
store_dt = 0.5
time_scaling = store_dt / dt

run_t = np.arange(0,t_endpoint,dt)
store_t = np.arange(0,t_endpoint,store_dt)

run_timepoint_N = len(run_t)
store_timepoint_N = len(store_t)

# initialize V_PLC, different value in each hex
params["V_PLC"] = allocate_var_dict(hex_array, 1, 0.787)
params["V_PLC"] = initialize_column_of_hexes_to_value_2(params["V_PLC"], hex_array, 0.9, 0, 1, pointy)
params["V_PLC"] = initialize_column_of_hexes_to_value_2(params["V_PLC"], hex_array, 0.9, 0.5, 1, pointy)
params["V_PLC"] = initialize_column_of_hexes_to_value_2(params["V_PLC"], hex_array, 0.9, 1, 1, pointy)

# temporarily turn off cell-cell communication
params["D_IP3"] = 0

# allocation initial conditions for variables
Ca_cyt_0 = 2
Ca_cyt = allocate_var_dict(hex_array, store_timepoint_N, Ca_cyt_0)
ip3 = allocate_var_dict(hex_array, store_timepoint_N, 0.2)
Ca_stored = allocate_var_dict(hex_array, store_timepoint_N, params["c_tot"] - Ca_cyt_0)
ip3R_act = allocate_var_dict(hex_array, store_timepoint_N, 0.6)

start = time.time()

variables = Ca_cyt, ip3, Ca_stored, ip3R_act,
Ca_cyt, ip3, Ca_stored, ip3R_act, = politi(variables, run_t, store_t, hex_array, params)

end = time.time()

print('Time solving', end - start)
    
################################################################################################## 
# PLOT

plot_vars = [Ca_cyt, ip3]
var_strings = ['Ca_cyt', 'ip3']
color_strings = ['Blues', 'Reds']

# save figs
# plot_hexes(hex_array, (hex_x_N,hex_y_N), pointy, 12, save_dir)
# selected_t_idx = 0
# plot_var_by_color(var_dict, selected_t_idx, hex_array, (hex_x_N,hex_y_N), pointy, 12, save_dir)
for var, var_str, color_str in zip(variables, var_strings, color_strings):
    animate_var_by_color(var, store_timepoint_N, hex_array, (hex_x_N,hex_y_N), pointy, 3, color_str, save_dir + var_str)

anim_time = time.time()
print('Time animating', anim_time - end)

################################################################################################## 
# PICKLE
pickle_dir = save_dir + "pickles/"
if not path.isdir(pickle_dir):
    mkdir(pickle_dir)
for var, var_str in zip(variables, var_strings):
    value_loc_tuple = create_val_loc_tuple_std_layout(var, hex_array, pointy)
    with open(pickle_dir + var_str + '.pickle', 'wb') as handle:
        pickle.dump(value_loc_tuple, handle)
with open(pickle_dir + 'layout_dict.pickle', 'wb') as handle:
    pickle.dump(layout_dict, handle)
hex_tuples = [hex_to_tuple(hexa) for hexa in hex_array]   
with open(pickle_dir + 'hex_tuples.pickle', 'wb') as handle:
    pickle.dump(hex_tuples, handle)
    
total = time.time()
print('Total time', total - start)
