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
# save_dir = results_dir + 'dev/'
if not path.isdir(save_dir):
    mkdir(save_dir)
    
# save code
code_dir = save_dir + "code/"
if not path.isdir(code_dir):
    mkdir(code_dir)
filenames = ['main.py', 'functions.py', 'plot.py', 'models.py', 'params.py']
for filename in filenames:
    copy2(filename, code_dir + filename)

################################################################################################## 
# set up layout
radius = 10.0
layout_dict = {
    "layout_str": 'pointy',
    "center_x": 0,
    "center_y": 0,
    "radius": radius
}
flat_layout_dict = {
    "layout_str": 'flat',
    "center_x": 0,
    "center_y": 0,
    "radius": radius / np.sqrt(3/4)
}

pointy = create_layout_from_dict(layout_dict)
flat = create_layout_from_dict(flat_layout_dict)

##################################################################################################
# set up hexagonal grid with (q,r,s) coordinates
hex_x_N = 40
hex_y_N = 6
hex_array = []

# set up POINTY hexagonal grid with (q,r,s) coordinates
# PARALLELOGRAM MAP
# for x in range(hex_x_N):
#     for y in range(hex_y_N):
#         hex_array.append(Hex(-x-y,y,x))

# RECTANGLE MAP
for y in range(hex_y_N):
    y_offset = int(np.floor(y/2))
    for x in range(-y_offset, hex_x_N - y_offset):
        hex_array.append(Hex(-x-y,y,x))
        
# set up FLAT hexagonal grid with (q,r,s) coordinates
# RECTANGLE MAP
# for x in range(hex_x_N):
#     x_offset = int(np.floor(x/2))
#     for y in range(-x_offset, hex_y_N - x_offset):
#         hex_array.append(Hex(x,y,-x-y))
        
##################################################################################################

t_endpoint = 300

dt = 0.05
store_dt = 0.5
time_scaling = store_dt / dt

run_t = np.arange(0,t_endpoint,dt)
store_t = np.arange(0,t_endpoint,store_dt)

run_timepoint_N = len(run_t)
store_timepoint_N = len(store_t)

# initialize V_PLC, different value in each hex
params["V_PLC"] = allocate_var_dict(hex_array, 1, 0.787)
# params["V_PLC"] = initialize_column_of_hexes_to_value_2(params["V_PLC"], hex_array, 0.85, 0, 1, pointy)

# temporarily turn off cell-cell communication
params["D_IP3"] = 0

# allocation initial conditions for variables
Ca_cyt_0 = 2
Ca_cyt = allocate_var_dict(hex_array, store_timepoint_N, Ca_cyt_0)
Ca_stored = allocate_var_dict(hex_array, store_timepoint_N, params["c_tot"] - Ca_cyt_0)

# # gradient in calcium initial conditions
# Ca_cyt = initialize_var_dict_to_x_gradient(Ca_cyt, hex_array, (0,Ca_cyt_0), pointy)
# for hexa in hex_array:
#     Ca_stored[hexa][0] = params["c_tot"] - Ca_cyt[hexa][0]
    
ip3 = allocate_var_dict(hex_array, store_timepoint_N, 0.2)
# gradient in calcium initial conditions
# ip3 = initialize_var_dict_to_x_gradient(Ca_cyt, hex_array, (0,1.2), pointy)

ip3R_act = allocate_var_dict(hex_array, store_timepoint_N, 0.6)
ip3R_act = initialize_var_dict_to_x_gradient(ip3R_act, hex_array, (0,1.2), pointy)

start = time.time()

variables = Ca_cyt, ip3, Ca_stored, ip3R_act,

Ca_cyt_new, ip3_new, Ca_stored_new, ip3R_act_new, = Ca_cyt, ip3, Ca_stored, ip3R_act,
# Ca_cyt_new, ip3_new, Ca_stored_new, ip3R_act_new, = politi(variables, run_t, store_t, hex_array, params)

end = time.time()

print('Time solving', end - start)

##################################################################################################
# PLOT

plot_vars = [Ca_cyt_new, ip3_new]
plot_var_strings = ['Ca_cyt', 'ip3']
color_strings = ['Blues', 'Reds']

# plt.plot(store_t, ip3_new[hex_array[0]])
# plt.show()

# save figs
plot_hexes(hex_array, (hex_x_N,hex_y_N), flat, 12, save_dir)
# selected_t_idx = 0

# animate_var_by_color(Ca_cyt_new, store_timepoint_N, hex_array, (hex_x_N,hex_y_N), pointy, 12, 'Blues', save_dir + 'Ca_cyt')
# animate_var_by_color(ip3_new, store_timepoint_N, hex_array, (hex_x_N,hex_y_N), pointy, 12, 'Reds', save_dir + 'ip3')

anim_time = time.time()
print('Time animating', anim_time - end)

##################################################################################################
# PICKLE

pickle_vars = [Ca_cyt_new, ip3_new, Ca_stored_new, ip3R_act_new,]
pickle_var_strings = ['Ca_cyt', 'ip3', 'Ca_stored', 'ip3R_act']

pickle_dir = save_dir + "pickles/"
if not path.isdir(pickle_dir):
    mkdir(pickle_dir)
for var, var_str in zip(pickle_vars, pickle_var_strings):
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
