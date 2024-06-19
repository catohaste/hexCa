##################################################################################################
"""
Potential issues
What if I want to make more than one connection at each time point?

Tasks
Make politi model use SDEs

"""
##################################################################################################
import sys
from os import mkdir, path
import random
import pickle
import datetime
from shutil import copy2
import time
import pandas as pd
import networkx as nx
from copy import deepcopy

import matplotlib.pyplot as plt

from lib import *
from functions import *
from models import *
from plot import *

from params import params

start = time.time()

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
filenames = ['main.py', 'functions.py', 'plot.py', 'models.py', 'params.py', 'lib.py']
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
# hex_x_N = 60
# hex_y_N = 1
# hex_x_N = 120
# hex_y_N = 12
# hex_x_N = 60
# hex_y_N = 6
hex_x_N = 50
hex_y_N = 30
hex_array = []

dimension_suffix_str = '_' + str(hex_x_N) + '_' + str(hex_y_N)
dimension_prefix_str = str(hex_x_N) + '_' + str(hex_y_N) + '_'

# set up POINTY hexagonal grid with (q,r,s) coordinates
# PARALLELOGRAM MAP
# for x in range(hex_x_N):
#     for y in range(hex_y_N):
#         hex_array.append(Hex(-x-y,y,x))

# RECTANGLE MAP
# ODD-R
for x in range(hex_x_N):
    for y in range(hex_y_N):
        h = OffsetCoord(x,y)
        hexa = roffset_to_cube(-1, h)
        hex_array.append(hexa)
# for y in range(hex_y_N):
#     y_offset = int(np.floor(y/2))
#     for x in range(-y_offset, hex_x_N - y_offset):
#         hex_array.append(Hex(-x-y,y,x))
        
# set up FLAT hexagonal grid with (q,r,s) coordinates
# RECTANGLE MAP
# for x in range(hex_x_N):
#     x_offset = int(np.floor(x/2))
#     for y in range(-x_offset, hex_y_N - x_offset):
#         hex_array.append(Hex(x,y,-x-y))

##################################################################################################
t_endpoint = 1200

dt = 0.05
store_dt = 0.5
time_scaling = store_dt / dt

run_t = np.arange(0,t_endpoint,dt)
store_t = np.arange(0,t_endpoint,store_dt)

run_timepoint_N = len(run_t)
store_timepoint_N = len(store_t)

# initialize V_PLC, different value in each hex
# params["V_PLC"] = allocate_var_dict(hex_array, 1, 0.787)
# params["V_PLC"] = initialize_column_of_hexes_to_value_2(params["V_PLC"], hex_array, 0.9, 0, 1, pointy)
params["V_PLC"] = initialize_var_dict_to_x_gradient(params["V_PLC"], hex_array, (0.9,0.787), pointy)
# print(params["V_PLC"])

plot_V_PLC(params["V_PLC"], hex_array, (hex_x_N,hex_y_N), pointy, 12, 'Reds', save_dir + 'V_PLC' + '_initial', show_time=False)

# set cell-cell communication, 0 => OFF, standard 0.02
params["D_IP3"] = 0.08

# allocation initial conditions for variables
Ca_cyt_0 = 2
Ca_cyt = allocate_var_dict(hex_array, store_timepoint_N, Ca_cyt_0)
Ca_stored = allocate_var_dict(hex_array, store_timepoint_N, params["c_tot"] - Ca_cyt_0)
ip3 = allocate_var_dict(hex_array, store_timepoint_N, 0.2)
ip3R_act = allocate_var_dict(hex_array, store_timepoint_N, 0.6)

variables = Ca_cyt, ip3, Ca_stored, ip3R_act,

##################################################################################################
""" SET INITIAL CONDITIONS """

# set ICs randomly from V_PLC 0.787 df 
ICs_df = pd.read_csv("ICs.csv", delimiter=',', header=0, index_col=0)
# variables = set_initial_conditions_from_df(ICs_df, variables)
# create_less_random_initial_conditions_df(hex_array, ICs_df, 'not_random_ICs' + dimension_suffix_str + '.csv')

# ip3R_act = initialize_var_dict_to_x_gradient(ip3R_act, hex_array, (0.435,0.845), pointy)
# ip3R_act = initialize_var_dict_to_random_val_in_range(ip3R_act, (0,1.2))
# ip3R_act (min, max) = (0.43397934441757985 0.8460908021227952) for V_PLC 0.787

# set ICs less randomly from V_PLC 0.787 df
# these conditions were selected randomly but now are the same for every run
less_random_ICs_df = pd.read_csv('not_random_ICs' + dimension_suffix_str + '.csv', delimiter=',', header=0, index_col=0)
variables = set_initial_conditions_from_df_less_random(less_random_ICs_df, variables)

initialize_time = time.time()

print('Time initializing', initialize_time - start)

##################################################################################################
""" SET CONNECTIONS """
""" Connections are an INPUT """
""" A connection between cells allows IP3 to travel between these cells. See connection demo. """

connection_params = {
    'neighbour_dist_limit': 1, # how far away can I connect
    'init_avg_degree_fraction' : 1, # average fraction of potential connections at start 
    'boundary_conditions' : 'no-flux', # flux or no-flux
    
    # both values should be a integers, and relate to the number of store_dt s
    # they should be less than store_timepoint_N
    'constant_connections' : True,
    'birth_connect_dt' : t_endpoint, # connection birth rate: new connection every x timesteps
    'death_connect_dt' : t_endpoint # connection death rate: lose connection every x timesteps
    
    # both values should be a integers, and relate to the number of store_dt s
    # they should be less than store_timepoint_N
    # 'constant_connections' : False,
    # 'birth_connect_dt' : 1, # connection birth rate: new connection every x timesteps
    # 'death_connect_dt' : 4 # connection death rate: lose connection every x timesteps
}

# create connections randomly and write to file
connections_dict = create_connections(hex_array, store_timepoint_N, connection_params)
dir_name = 'connections'
if not path.isdir(dir_name):
    mkdir(dir_name)
# write_connections_to_file(connections_dict, path.join(dir_name, dimension_prefix_str + 'full_'))

connections_dict = load_connections_from_file(hex_array, path.join('connections', dimension_prefix_str + 'full_'))

initial_connections = connections_dict['initial_connections']
birth_connections = connections_dict['birth_connections']
death_connections = connections_dict['death_connections']

# beginning connections
plot_graph_fixed_time(initial_connections, hex_array, (hex_x_N,hex_y_N), pointy, 12, "Blues", save_dir + 'connections_start')

# end connections
if connection_params['constant_connections']:
    end_connections = initial_connections
else:
    birth_t, death_t, birth_and_death_t = get_birth_and_death_t(store_timepoint_N, connection_params)
    end_connect_timepoint = max(birth_and_death_t)
    end_connections = get_current_connections(end_connect_timepoint, initial_connections, birth_connections, death_connections)
    
plot_graph_fixed_time(end_connections, hex_array, (hex_x_N,hex_y_N), pointy, 12, "Blues", save_dir + 'connections_end')

print('initial', initial_connections)
print('end', end_connections)

connect_time = time.time()

print('Time connecting', connect_time - initialize_time)

# print("size initial", sys.getsizeof(initial_connections))
# print("size potential", sys.getsizeof(potential_connections))
# print("size full", sys.getsizeof(store_cell_connections))

##################################################################################################
""" RUN """

start = time.time()

Ca_cyt_new, ip3_new, Ca_stored_new, ip3R_act_new, = Ca_cyt, ip3, Ca_stored, ip3R_act,
# Ca_cyt_new, ip3_new, Ca_stored_new, ip3R_act_new, = politi(variables, run_t, store_t, hex_array, params)
Ca_cyt_new, ip3_new, Ca_stored_new, ip3R_act_new, = politi_reduced_connectivity(variables, run_t, store_t, hex_array, initial_connections, birth_connections, death_connections, params)

solv_time = time.time()

print('Time solving', solv_time - connect_time)

##################################################################################################
""" PICKLE """

pickle_dir = save_dir + "pickles/"
if not path.isdir(pickle_dir):
    mkdir(pickle_dir)

time_params = {
    "endpoint": t_endpoint,
    'dt': dt,
    'store_dt': store_dt,
    'birth_connect_dt': connection_params['birth_connect_dt'],
    'death_connect_dt': connection_params['death_connect_dt']
}

with open(pickle_dir + 'time_params.pickle', 'wb') as handle:
    pickle.dump(time_params, handle)
pickle_vars = [Ca_cyt_new, ip3_new, Ca_stored_new, ip3R_act_new,]
pickle_var_strings = ['Ca_cyt', 'ip3', 'Ca_stored', 'ip3R_act']

for var, var_str in zip(pickle_vars, pickle_var_strings):
    value_loc_tuple = create_val_loc_tuple_std_layout(var, hex_array, pointy)
    with open(pickle_dir + var_str + '.pickle', 'wb') as handle:
        pickle.dump(value_loc_tuple, handle)
        
with open(pickle_dir + 'layout_dict.pickle', 'wb') as handle:
    pickle.dump(layout_dict, handle)
hex_tuples = [hex_to_tuple(hexa) for hexa in hex_array]
with open(pickle_dir + 'hex_tuples.pickle', 'wb') as handle:
    pickle.dump(hex_tuples, handle)
    
pickle_params = deepcopy(params)
pickle_params.pop('V_PLC')
with open(pickle_dir + 'params.pickle', 'wb') as handle:
    pickle.dump(pickle_params, handle)
    
V_PLC_value_loc_tuple = create_val_loc_tuple_std_layout(params['V_PLC'], hex_array, pointy)
with open(pickle_dir + 'V_PLC' + '.pickle', 'wb') as handle:
    pickle.dump(V_PLC_value_loc_tuple, handle)
    
with open(pickle_dir + 'connect_params.pickle', 'wb') as handle:
    pickle.dump(connection_params, handle)
    
connect_str = str(connection_params['init_avg_degree_fraction']) + '_'
write_connections_to_file(connections_dict, path.join(pickle_dir, dimension_prefix_str + connect_str))

pickling_time = time.time()
print('Time pickling', pickling_time - solv_time)

##################################################################################################
""" MAKE LINKS """
""" Links are an OUTPUT """
""" A link is created when neighbouring cells 'fire' (amplitude above threshold) within a given time interval. """

# Ca_cyt_links = make_links(Ca_cyt_new, store_t)
# ip3_links = make_links(ip3_new, store_t)

# print(len(ip3_links))
# print(ip3_links)
# for t_idx, t_links in enumerate(links):
#     print(t_idx, len(t_links))

link_time = time.time()
print('Time linking', link_time - pickling_time)

##################################################################################################
""" PLOT """

# animate connections
animate_connections(initial_connections, birth_connections, death_connections, hex_array, (hex_x_N,hex_y_N), pointy, 12, "Blues", save_dir + 'connections')

plot_vars = [Ca_cyt_new, ip3_new]
# link_vars = [Ca_cyt_links, ip3_links]
plot_var_strings = ['Ca_cyt', 'IP3']
color_strings = ['Oranges', 'Blues']

plot_var_strings = ['IP3']
color_strings = ['Blues']
plot_time_indicies = [int(x/store_dt)for x in [0, 488, 494, 500, 506, 512, 518, 524, 530, 536, 542, 548]]

for var, var_str, color_str in zip(plot_vars, plot_var_strings, color_strings):
    animate_var_by_color(var, store_timepoint_N, store_dt, hex_array, (hex_x_N,hex_y_N), pointy, 5, var_str, color_str, save_dir + var_str, show_time=True, show_colorbar=True, D_str=str(params['D_IP3']))
    plot_colorbar(var, 5, var_str, color_str, save_dir + var_str + '_colorbar')
    # for plot_time in plot_time_indicies:
    #     plot_var_by_color(var, plot_time, store_dt, hex_array, (hex_x_N,hex_y_N), pointy, 12, color_str, save_dir+ var_str + '_' + str(int(plot_time*store_dt)), show_time=True)

# plot_vars = [Ca_cyt_new, ip3_new, Ca_stored_new, ip3R_act_new]
# plot_var_strings = ['Ca_cyt', 'IP3', 'Ca_ER', 'IP3R_active']
# color_strings = ['Blues', 'Oranges', 'Greens', 'Reds']
 
# save figs
# plot_hexes(hex_array, (hex_x_N,hex_y_N), flat, 12, save_dir)
# plot_var_by_color(Ca_cyt_new, 0, hex_array, (hex_x_N,hex_y_N), pointy, 12, save_dir + 'Ca_cyt' + '_intial')

new_variables = Ca_cyt_new, ip3_new, Ca_stored_new, ip3R_act_new,
all_var_strings = ['Ca_cyt', 'IP3', 'Ca_ER', 'IP3R_active']
all_color_strings = ['Blues', 'Oranges', 'Greens', 'Reds']

# chosen_cells = [(4,4), (12,4), (20,4), (28,4), (36,4)]
# chosen_cells = [(0,0)]
# plot_hexes_highlight_cells(chosen_cells, hex_array, (hex_x_N,hex_y_N), pointy, 12, save_dir)
# plot_all_vars_over_time_single_cells(chosen_cells, new_variables, all_var_strings, all_color_strings, save_dir)
#
# for var, link_var, var_str, color_str in zip(plot_vars, link_vars, plot_var_strings, color_strings):
#     animate_var_by_color(var, store_timepoint_N, hex_array, (hex_x_N,hex_y_N), pointy, 12, color_str, save_dir + var_str)
#     animate_var_over_x_avg_y(var, store_timepoint_N, hex_array, (hex_x_N,hex_y_N), pointy, 12, color_str, var_str, save_dir + var_str)
#     plot_var_over_time_fixed_x_avg_y(var, hex_array, pointy, 12, color_str, var_str, save_dir + var_str)
#     plot_links(link_var, hex_array, (hex_x_N,hex_y_N), pointy, 12, color_str, save_dir + var_str + "_links" )
#     plot_var_running_time_avg_single_cells(chosen_cells, 400, var, var_str, color_str, save_dir + var_str)

# plot_initial_graph(store_cell_connections[0], hex_array, (hex_x_N,hex_y_N), pointy, 12, "Oranges", save_dir + 'connections_' + str(0))
# demo_connections(connection_params, pointy)

anim_time = time.time()
print('Time animating', anim_time - link_time)

##################################################################################################

total_time = time.time()
print('Total time', total_time - start)


