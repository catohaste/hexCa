import sys
from os import mkdir, path
import random
import pickle
import datetime
from shutil import copy2
import time
import pandas as pd
import networkx as nx

import matplotlib.pyplot as plt

from lib import *
from functions import *
from models import *
from plot_a_model_in_development import *

from params import params

##################################################################################################
# save directory
save_dir = 'a_model_in_development/'
if not path.isdir(save_dir):
    mkdir(save_dir)
    
# save code
code_dir = save_dir + "code/"
if not path.isdir(code_dir):
    mkdir(code_dir)
filenames = ['a_model_in_development.py', 'functions.py', 'plot.py','plot_a_model_in_development.py',  'models.py', 'params.py', 'lib.py']
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

pointy = create_layout_from_dict(layout_dict)

##################################################################################################
# set up hexagonal grid with (q,r,s) coordinates
hex_x_N = 65
hex_y_N = 5

stripes = {
    'black': [],
    'brown': [],
    'red': [],
    'orange': [],
    'yellow': [],
    'green': [],
    'blue': [],
    'purple': [],
}

stripes = {
    'black': [],
    'brown': [],
    'red': [],
    'orange': [],
    'yellow': [],
    'green': []
}

# stripes = {
#     'green': []
# }

# symmetric
stripe_offset = 0
for stripe in stripes:
    
    for x in range(hex_x_N - 1):
        for y in range(stripe_offset, stripe_offset + hex_y_N):
            h = OffsetCoord(x,y)
            hexa = roffset_to_cube(-1, h)
            stripes[stripe].append(hexa)
            
    if stripe_offset % 2 == 0:

        for y in range(stripe_offset, stripe_offset + hex_y_N, 2):
            h = OffsetCoord(hex_x_N - 1,y)
            hexa = roffset_to_cube(-1, h)
            stripes[stripe].append(hexa)
            
    else:
        
        for y in range(stripe_offset + 1, stripe_offset + hex_y_N, 2):
            h = OffsetCoord(hex_x_N - 1,y)
            hexa = roffset_to_cube(-1, h)
            stripes[stripe].append(hexa)
        
    stripe_offset -= hex_y_N

##################################################################################################

t_endpoint = 600

dt = 0.05
store_dt = 0.5
time_scaling = store_dt / dt

run_t = np.arange(0,t_endpoint,dt)
store_t = np.arange(0,t_endpoint,store_dt)

run_timepoint_N = len(run_t)
store_timepoint_N = len(store_t)

##################################################################################################
# allocate all stripes

stripes_var = {}

##################################################################################################

# constant (draw hexagons)
if 'black' in stripes:
    black_var = allocate_var_dict(stripes['black'], store_timepoint_N, 0.8)
    
    stripes_var['black'] = black_var
    print("Done BLACK")

# constant random (colour hexagons)
if 'brown' in stripes:
    brown_var = allocate_var_dict(stripes['brown'], store_timepoint_N, 0.6)
    brown_var = set_var_dict_to_random_val_in_range_all_t(brown_var, [0,1])
    
    stripes_var['brown'] = brown_var
    print("Done BROWN")

# random flashing (animate hexagons)
if 'red' in stripes:
    red_var = allocate_var_dict(stripes['red'], store_timepoint_N, 0.6)
    red_var = initialize_var_dict_to_random_val_in_range(red_var, [0,1])
    red_var = random_pulsing(red_var, stripes['red'], store_timepoint_N, (0,1))
    
    stripes_var['red'] = red_var
    print("Done RED")

# diffusion (something meaningful)
if 'orange' in stripes:
    orange_var = allocate_var_dict(stripes['orange'], store_timepoint_N, 0)
    orange_var = initialize_column_of_hexes_to_value(orange_var, stripes['orange'], 1, 0, 3, pointy)
    orange_var = initialize_column_of_hexes_to_value(orange_var, stripes['orange'], 1, 0.2, 3, pointy)
    orange_var = initialize_column_of_hexes_to_value(orange_var, stripes['orange'], 1, 0.4, 3, pointy)
    orange_var = initialize_column_of_hexes_to_value(orange_var, stripes['orange'], 1, 0.6, 3, pointy)
    orange_var = initialize_column_of_hexes_to_value(orange_var, stripes['orange'], 1, 0.8, 3, pointy)
    orange_var = initialize_column_of_hexes_to_value(orange_var, stripes['orange'], 1, 1, 3, pointy)
    orange_var = diffusion(orange_var, stripes['orange'], store_timepoint_N)
    
    stripes_var['orange'] = orange_var
    print("Done ORANGE")

##################################################################################################
""" YELLOW """
# calcium no diffusion (more meaningful) with graded starting conditions forming sweep
if 'yellow' in stripes:
    
    # initialize V_PLC, different value in each hex
    params["V_PLC"] = allocate_var_dict(stripes['yellow'], 1, 0.787)

    # set cell-cell communication, 0 => OFF, standard 0.02
    params["D_IP3"] = 0

    # allocation initial conditions for variables
    Ca_cyt_0 = 2
    Ca_cyt = allocate_var_dict(stripes['yellow'], store_timepoint_N, Ca_cyt_0)
    Ca_stored = allocate_var_dict(stripes['yellow'], store_timepoint_N, params["c_tot"] - Ca_cyt_0)
    ip3 = allocate_var_dict(stripes['yellow'], store_timepoint_N, 0.2)
    ip3R_act = allocate_var_dict(stripes['yellow'], store_timepoint_N, 0.6)

    variables = Ca_cyt, ip3, Ca_stored, ip3R_act,

    # set ICs randomly from V_PLC 0.787 df
    ICs_df = pd.read_csv("ICs.csv", delimiter=',', header=0, index_col=0)
    variables = initialize_var_dict_to_x_gradient_from_ICs_df(ICs_df, variables, stripes['yellow'], pointy)

    Ca_cyt_new, ip3_new, Ca_stored_new, ip3R_act_new, = Ca_cyt, ip3, Ca_stored, ip3R_act,
    Ca_cyt_new, ip3_new, Ca_stored_new, ip3R_act_new, = politi(variables, run_t, store_t, stripes['yellow'], params)

    yellow_var = Ca_cyt_new
    
    stripes_var['yellow'] = yellow_var
    print("Done YELLOW")
    
##################################################################################################
""" GREEN """
# calcium no diffusion with constant starting conditions but varied frewquency, forming sweep
if 'green' in stripes:

    # initialize V_PLC, different value in each hex
    params["V_PLC"] = initialize_var_dict_to_x_gradient(params["V_PLC"], stripes['green'], (0.787,1.1), pointy)

    # set cell-cell communication, 0 => OFF, standard 0.02
    params["D_IP3"] = 0

    # allocation initial conditions for variables
    Ca_cyt_0 = 2
    Ca_cyt = allocate_var_dict(stripes['green'], store_timepoint_N, Ca_cyt_0)
    Ca_stored = allocate_var_dict(stripes['green'], store_timepoint_N, params["c_tot"] - Ca_cyt_0)
    ip3 = allocate_var_dict(stripes['green'], store_timepoint_N, 0.2)
    ip3R_act = allocate_var_dict(stripes['green'], store_timepoint_N, 0.6)

    variables = Ca_cyt, ip3, Ca_stored, ip3R_act,

    # set ICs randomly from V_PLC 0.787 df
    ICs_df = pd.read_csv("ICs.csv", delimiter=',', header=0, index_col=0)
    variables = set_constant_initial_conditions_from_df_first_row(ICs_df, variables)

    Ca_cyt_new, ip3_new, Ca_stored_new, ip3R_act_new, = Ca_cyt, ip3, Ca_stored, ip3R_act,
    Ca_cyt_new, ip3_new, Ca_stored_new, ip3R_act_new, = politi(variables, run_t, store_t, stripes['green'], params)

    green_var = Ca_cyt_new

    stripes_var['green'] = green_var
    print("Done GREEN")

##################################################################################################
# """ GREEN """
# # calcium with diffusion (turing pattern)
# if 'green' in stripes:
#
#     # initialize V_PLC, same value in each hex
#     params["V_PLC"] = allocate_var_dict(stripes['green'], 1, 0.787)
#
#     # set cell-cell communication, 0 => OFF, standard 0.02
#     params["D_IP3"] = 0.02
#
#     # allocation initial conditions for variables
#     Ca_cyt_0 = 2
#     Ca_cyt = allocate_var_dict(stripes['green'], store_timepoint_N, Ca_cyt_0)
#     Ca_stored = allocate_var_dict(stripes['green'], store_timepoint_N, params["c_tot"] - Ca_cyt_0)
#     ip3 = allocate_var_dict(stripes['green'], store_timepoint_N, 0.2)
#     ip3R_act = allocate_var_dict(stripes['green'], store_timepoint_N, 0.6)
#
#     variables = Ca_cyt, ip3, Ca_stored, ip3R_act,
#
#     # set ICs randomly from V_PLC 0.787 df
#     ICs_df = pd.read_csv("ICs.csv", delimiter=',', header=0, index_col=0)
#     variables = set_initial_conditions_from_df(ICs_df, variables)
#
#     Ca_cyt_new, ip3_new, Ca_stored_new, ip3R_act_new, = Ca_cyt, ip3, Ca_stored, ip3R_act,
#     Ca_cyt_new, ip3_new, Ca_stored_new, ip3R_act_new, = politi(variables, run_t, store_t, stripes['green'], params)
#
#     green_var = Ca_cyt_new
#
#     stripes_var['green'] = green_var
#     print("Done GREEN")

##################################################################################################
""" BLUE """
# animated reduced connectivity
if 'blue' in stripes:
    
    blue_var = allocate_var_dict(stripes['blue'], store_timepoint_N, 0.6)
    
    # connection params
    neighbour_dist_limit = 1 # how far away can I connect
    init_avg_degree_fraction = 1 # average fraction of potential connections
    boundary_conditions = 'no-flux' # flux or no-flux

    # both values should be a multiple of store_dt and less than t_endpoint
    constant_connections = False
    birth_connect_dt = 2 # connection birth rate: new connection every x timesteps
    death_connect_dt = 4 # connection death rate: lose connection every x timesteps

    connection_params = {
        'dist_limit': neighbour_dist_limit,
        'init_avg_degree_fraction': init_avg_degree_fraction,
        'boundary_conditions': boundary_conditions
    }

    """potential"""
    potential_connections = nx.Graph()
    potential_connections.add_nodes_from(hex_array)
    for hexa in hex_array:
        neighbors = hex_neighbors_cumulative_distance(hexa, neighbour_dist_limit)
        for neighbor in neighbors:
            if neighbor in hex_array:
                potential_connections.add_edge(hexa, neighbor)

    potential_deg = list(dict(potential_connections.degree()).values())
    # print(min(potential_deg), max(potential_deg), np.mean(potential_deg))
    print("potential", potential_connections)

    """initial"""
    if init_avg_degree_fraction == 1:
        initial_connections = potential_connections
        birth_connections = {}
        death_connections = {}
    else:
        initial_connections = nx.Graph()
        initial_connections.add_nodes_from(hex_array)

        print('before', get_mean_degree_fraction(initial_connections, potential_connections))
        while get_mean_degree_fraction(initial_connections, potential_connections) < init_avg_degree_fraction:
            remaining_possible_connections = nx.difference(potential_connections, initial_connections)
            random_edge = random.sample(list(remaining_possible_connections.edges), 1)[0]
            initial_connections.add_edge(random_edge[0], random_edge[1])
        print('after',get_mean_degree_fraction(initial_connections, potential_connections))


        """birth and death"""
        birth_connections = {}
        death_connections = {}

        # allocate
        birth_t = range(birth_connect_dt, t_endpoint, birth_connect_dt)
        for t in birth_t:
            birth_connections[t] = []
        death_t = range(death_connect_dt, t_endpoint, death_connect_dt)
        for t in death_t:
            death_connections[t] = []
    
        # set up connections
        birth_and_death_t = list(set(list(birth_t) + list(death_t)))
        birth_and_death_t.sort()
        for current_t in birth_and_death_t:
            current_connections = get_current_connections(current_t, initial_connections, birth_connections, death_connections)
            if current_t in birth_t:
                possible_birth_connections = nx.difference(potential_connections, current_connections)
                random_edge = random.sample(list(possible_birth_connections.edges), 1)[0]
                birth_connections[current_t].append(random_edge)
            if current_t in death_t:
                random_edge = random.sample(list(current_connections.edges), 1)[0]
                death_connections[current_t].append(random_edge)
    
    stripes_var['blue'] = blue_var
    print("Done BLUE")


##################################################################################################
""" PURPLE """
# ????????

if 'purple' in stripes:

    purple_var = allocate_var_dict(stripes['purple'], store_timepoint_N, 0.6)
    purple_var = set_var_dict_to_random_val_in_range_all_t(purple_var, [0,1])
    
    stripes_var['purple'] = purple_var
    print("Done PURPLE")

##################################################################################################
""" PLOT """

flag_colours = ''
for stripe in stripes:
    flag_colours = flag_colours + '_' + stripe
flag_name = 'flag' + flag_colours
    
plot_var_flag(stripes_var, 150, stripes, (hex_x_N, hex_y_N), pointy, 12, save_dir + flag_name + '_initial')

animate_var_flag(stripes_var, store_timepoint_N, stripes, (hex_x_N,hex_y_N), pointy, 12, save_dir + flag_name)

##################################################################################################
# """ SET INITIAL CONDITIONS """
#
# # set ICs randomly from V_PLC 0.787 df
# ICs_df = pd.read_csv("ICs.csv", delimiter=',', header=0, index_col=0)
# # variables = set_initial_conditions_from_df(ICs_df, variables)
#
# # ip3R_act = initialize_var_dict_to_x_gradient(ip3R_act, hex_array, (0.435,0.845), pointy)
# # ip3R_act = initialize_var_dict_to_random_val_in_range(ip3R_act, (0,1.2))
# #ip3R_act (min, max) = (0.43397934441757985 0.8460908021227952) for V_PLC 0.787
#
# # set ICs less randomly from V_PLC 0.787 df
# # these conditions were selected randomly but now are the same for every run
# less_random_ICs_df = pd.read_csv("not_random_ICs.csv", delimiter=',', header=0, index_col=0)
# variables = set_initial_conditions_from_df_less_random(less_random_ICs_df, variables)
#
# ##################################################################################################
# """ SET CONNECTIONS """
#
# # connection params
# neighbour_dist_limit = 2 # how far away can I connect
# init_avg_degree = 4 # average number of connections I'm aiming for
# connect_birth_rate = 4 # connection birth rate
# connect_death_rate = 2 # connection death rate
#
# # allocate
# store_cell_connections = []
# for t in range(store_timepoint_N):
#     store_cell_connections.append(nx.Graph())
#
# for t in range(store_timepoint_N):
#     store_cell_connections[t].add_nodes_from(hex_array)
#
# # initialize
# init_cell_connections = []
#
#
#
# print(cell_connections[0])
#
# ##################################################################################################
# """ RUN """
#
# start = time.time()
#

#
# solv_time = time.time()
#
# print('Time solving', solv_time - start)
#
# ##################################################################################################
# """ PICKLE """
#
# pickle_vars = [Ca_cyt_new, ip3_new, Ca_stored_new, ip3R_act_new,]
# pickle_var_strings = ['Ca_cyt', 'ip3', 'Ca_stored', 'ip3R_act']
#
# pickle_dir = save_dir + "pickles/"
# if not path.isdir(pickle_dir):
#     mkdir(pickle_dir)
# for var, var_str in zip(pickle_vars, pickle_var_strings):
#     value_loc_tuple = create_val_loc_tuple_std_layout(var, hex_array, pointy)
#     with open(pickle_dir + var_str + '.pickle', 'wb') as handle:
#         pickle.dump(value_loc_tuple, handle)
# with open(pickle_dir + 'layout_dict.pickle', 'wb') as handle:
#     pickle.dump(layout_dict, handle)
# hex_tuples = [hex_to_tuple(hexa) for hexa in hex_array]
# with open(pickle_dir + 'hex_tuples.pickle', 'wb') as handle:
#     pickle.dump(hex_tuples, handle)
#
# pickling_time = time.time()
# print('Time pickling', pickling_time - solv_time)
#
# ##################################################################################################
# """ MAKE LINKS """
#
# # Ca_cyt_links = make_links(Ca_cyt_new, store_t)
# # ip3_links = make_links(ip3_new, store_t)
#
# # print(len(links))
# # for t_idx, t_links in enumerate(links):
# #     print(t_idx, len(t_links))
#
# link_time = time.time()
# print('Time linking', link_time - pickling_time)
#
#
#
# ##################################################################################################
# """ PLOT """
#
# plot_vars = [Ca_cyt_new, ip3_new]
# link_vars = [Ca_cyt_links, ip3_links]
# plot_var_strings = ['Ca_cyt', 'IP3']
# color_strings = ['Blues', 'Oranges']
#
# # plot_vars = [Ca_cyt_new, ip3_new, Ca_stored_new, ip3R_act_new]
# # plot_var_strings = ['Ca_cyt', 'IP3', 'Ca_ER', 'IP3R_active']
# # color_strings = ['Blues', 'Oranges', 'Greens', 'Reds']
#
# # save figs
# # plot_hexes(hex_array, (hex_x_N,hex_y_N), flat, 12, save_dir)
# # plot_var_by_color(Ca_cyt_new, 0, hex_array, (hex_x_N,hex_y_N), pointy, 12, save_dir+ 'Ca_cyt' + '_intial')
#
# new_variables = Ca_cyt_new, ip3_new, Ca_stored_new, ip3R_act_new,
# all_var_strings = ['Ca_cyt', 'IP3', 'Ca_ER', 'IP3R_active']
# all_color_strings = ['Blues', 'Oranges', 'Greens', 'Reds']
#
# chosen_cells = [(4,4), (12,4), (20,4), (28,4), (36,4)]
# # plot_hexes_highlight_cells(chosen_cells, hex_array, (hex_x_N,hex_y_N), pointy, 12, save_dir)
# # plot_all_vars_over_time_single_cells(chosen_cells, new_variables, all_var_strings, all_color_strings, save_dir)
# #
# # for var, link_var, var_str, color_str in zip(plot_vars, link_vars, plot_var_strings, color_strings):
# #     animate_var_by_color(var, store_timepoint_N, hex_array, (hex_x_N,hex_y_N), pointy, 12, color_str, save_dir + var_str)
# #     animate_var_over_x_avg_y(var, store_timepoint_N, hex_array, (hex_x_N,hex_y_N), pointy, 12, color_str, var_str, save_dir + var_str)
# #     plot_var_over_time_fixed_x_avg_y(var, hex_array, pointy, 12, color_str, var_str, save_dir + var_str)
# #     plot_links(link_var, hex_array, (hex_x_N,hex_y_N), pointy, 12, color_str, save_dir + var_str + "_links" )
# #     plot_var_running_time_avg_single_cells(chosen_cells, 400, var, var_str, color_str, save_dir + var_str)
#
# anim_time = time.time()
# print('Time animating', anim_time - link_time)
#
# ##################################################################################################
# """ STORE """
# # # create dataframe for initial conditions
# # csv_out = {
# #     "time": store_t,
# #     "Ca_cyt": Ca_cyt_new[Hex(0,0,0)] ,
# #     "IP3": ip3_new[Hex(0,0,0)],
# #     "Ca_stored": Ca_stored_new[Hex(0,0,0)] ,
# #     "IP3R_act": ip3R_act_new[Hex(0,0,0)]
# # }
# #
# # pd.DataFrame(csv_out).to_csv(save_dir + "ICs_V_PLC_0787.csv")
# #
# # # create dataframe for less random initial conditions
# # new_col_names = ['q','r','s','time','Ca_cyt','IP3','Ca_stored','IP3R_act']
# # not_random_ICs_df = pd.DataFrame(columns=new_col_names)
# #
# # counter = 0
# # for hexa in hex_array:
# #     row = ICs_df.sample(1)
# #     new_row_values = [int(hexa.q), int(hexa.r), int(hexa.s)] +  row.values.flatten().tolist()
# #     new_dict = dict(zip(new_col_names, new_row_values))
# #     not_random_ICs_df = not_random_ICs_df.append(new_dict, ignore_index=True)
# #
# # not_random_ICs_df.to_csv('not_random_ICs.csv')
#
# total_time = time.time()
# print('Total time', total_time - start)
#
#
