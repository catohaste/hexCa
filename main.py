from os import mkdir, path
import random
import pickle
import datetime
from shutil import copy2

from lib import *
from functions import *
from models import *
from plot import plot_hexes, plot_var_by_color, animate_var_by_color

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
hex_x_N = 40
hex_y_N = 6
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

timepoint_N = 960
value_range = (0,100)

var_dict = allocate_var_dict(hex_array, timepoint_N, 0)
var_dict = initialize_leftmost_hexes_to_value(var_dict, hex_array, 100, pointy)

var_dict = diffusion(var_dict, hex_array, timepoint_N)
    
################################################################################################## 
# PLOT

# save figs
# plot_hexes(hex_array, (hex_x_N,hex_y_N), pointy, 12, save_dir)
# selected_t_idx = 0
# plot_var_by_color(var_dict, selected_t_idx, hex_array, (hex_x_N,hex_y_N), pointy, 12, save_dir)
animate_var_by_color(var_dict, timepoint_N, hex_array, (hex_x_N,hex_y_N), pointy, 12, save_dir)

################################################################################################## 
# PICKLE
pickle_dir = save_dir + "pickles/"
if not path.isdir(pickle_dir):
    mkdir(pickle_dir)
value_loc_tuples = create_val_loc_tuple_std_layout(var_dict, hex_array, pointy)
with open(pickle_dir + 'value_loc_tuples.pickle', 'wb') as handle:
    pickle.dump(value_loc_tuples, handle)
with open(pickle_dir + 'layout_dict.pickle', 'wb') as handle:
    pickle.dump(layout_dict, handle)
