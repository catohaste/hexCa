import sys
from os import mkdir, path
import random
import pickle
import datetime
from shutil import copy2
import time
import pandas as pd

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

# run_selection = ['black', 'brown', 'red', 'orange', 'yellow', 'green', 'blue', 'purple']
run_selection = ['yellow']

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
animation_type = {} # set to 'colour' or 'connect' for each stripe
connect_var = {}

##################################################################################################

# constant (draw hexagons)
if 'black' in run_selection:
    black_var = allocate_var_dict(stripes['black'], store_timepoint_N, 0.8)
    
    stripes_var['black'] = black_var
    animation_type['black'] = 'colour'
    print("Done BLACK")

# constant random (colour hexagons)
if 'brown' in run_selection:
    brown_var = allocate_var_dict(stripes['brown'], store_timepoint_N, 0.6)
    
    ##### INITIAL CONDITIONS ############
    # # set ICs randomly
    # brown_var = set_var_dict_to_random_val_in_range_all_t(brown_var, [0,1])
    #
    # # create less random ICs
    # brown_less_random_ICs_df = pd.DataFrame(columns=['q','r','s','value'])
    # for idx, hexa in enumerate(stripes['brown']):
    #      brown_less_random_ICs_df.loc[idx] = [hexa.q, hexa.r, hexa.s] + [brown_var[hexa][0]]
    # brown_less_random_ICs_df.to_csv(path_or_buf='not_random_ICs_brown.csv')
    
    # load less random ICs
    brown_ICs_df = pd.read_csv("not_random_ICs_brown.csv", delimiter=',', header=0, index_col=0)
    for idx, row in brown_ICs_df.iterrows():
        hexa = Hex(row['q'],row['r'],row['s'])
        for t in range(store_timepoint_N):
            brown_var[hexa][t] = row['value']
    
    stripes_var['brown'] = brown_var
    animation_type['brown'] = 'colour'
    print("Done BROWN")

# random flashing (animate hexagons)
if 'red' in run_selection:
    red_var = allocate_var_dict(stripes['red'], store_timepoint_N, 0.6)
    
    #### INITIAL CONDITIONS ############
    # # set ICs randomly
    # red_var = initialize_var_dict_to_random_val_in_range(red_var, [0,1])
    #
    # # create less random ICs
    # red_less_random_ICs_df = pd.DataFrame(columns=['q','r','s','value'])
    # for idx, hexa in enumerate(stripes['red']):
    #      red_less_random_ICs_df.loc[idx] = [hexa.q, hexa.r, hexa.s] + [red_var[hexa][0]]
    # red_less_random_ICs_df.to_csv(path_or_buf='not_random_ICs_red.csv')
    
    # load less random ICs
    red_ICs_df = pd.read_csv("not_random_ICs_red.csv", delimiter=',', header=0, index_col=0)
    for idx, row in red_ICs_df.iterrows():
        hexa = Hex(row['q'],row['r'],row['s'])
        red_var[hexa][0] = row['value']
    
    red_var = random_pulsing(red_var, stripes['red'], store_timepoint_N, (0,1))
    
    stripes_var['red'] = red_var
    animation_type['red'] = 'colour'
    print("Done RED")

# diffusion (something meaningful)
if 'orange' in run_selection:
    orange_var = allocate_var_dict(stripes['orange'], store_timepoint_N, 0)
    orange_var = initialize_column_of_hexes_to_value(orange_var, stripes['orange'], 1, 0, 3, pointy)
    orange_var = initialize_column_of_hexes_to_value(orange_var, stripes['orange'], 1, 0.2, 3, pointy)
    orange_var = initialize_column_of_hexes_to_value(orange_var, stripes['orange'], 1, 0.4, 3, pointy)
    orange_var = initialize_column_of_hexes_to_value(orange_var, stripes['orange'], 1, 0.6, 3, pointy)
    orange_var = initialize_column_of_hexes_to_value(orange_var, stripes['orange'], 1, 0.8, 3, pointy)
    orange_var = initialize_column_of_hexes_to_value(orange_var, stripes['orange'], 1, 1, 3, pointy)
    orange_var = diffusion(orange_var, stripes['orange'], store_timepoint_N)
    
    stripes_var['orange'] = orange_var
    animation_type['orange'] = 'colour'
    print("Done ORANGE")

##################################################################################################
""" YELLOW """
# calcium no diffusion (more meaningful) with graded starting conditions but constant V_PLC forming sweep
if 'yellow' in run_selection:
    
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
    animation_type['yellow'] = 'colour'
    print("Done YELLOW")
    
##################################################################################################
""" GREEN """
# calcium no diffusion with constant starting conditions but varied frequency (V_PLC), forming sweep
if 'green' in run_selection:

    # initialize V_PLC, different value in each hex
    params["V_PLC"] = initialize_var_dict_to_x_gradient(params["V_PLC"], stripes['green'], (1.1,0.787), pointy)

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
    animation_type['green'] = 'colour'
    print("Done GREEN")

##################################################################################################
# """ GREEN """
# # calcium with diffusion (turing pattern)
# if 'green' in run_selection:
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
# animated connections
if 'blue' in run_selection:
    
    blue_var = allocate_var_dict(stripes['blue'], store_timepoint_N, 0.6)
    
    connection_params = {
        'neighbour_dist_limit': 1, # how far away can I connect
        'init_avg_degree_fraction' : 0.1, # average fraction of potential connections at start 
        'boundary_conditions' : 'no-flux', # flux or no-flux
        
        'constant_connections' : False,
        
        # both values should be a integers, and relate to the number of store_dt s
        # they should be less than store_timepoint_N
        'birth_connect_dt' : 1, # connection birth rate: new connection every x timesteps
        'death_connect_dt' : 4 # connection death rate: lose connection every x timesteps
    }
    
    # create connections and write to file
    # connections_blue = create_connections(stripes['blue'], store_timepoint_N, connection_params)
    # dir_name = 'connections'
    # if not path.isdir(dir_name):
    #     mkdir(dir_name)
    # write_connections_to_file(connections_blue, path.join(dir_name, 'blue_'))
    
    # load connections from file
    connections_blue = load_connections_from_file(stripes['blue'], path.join('connections', 'blue_'))
    
    blue_var = allocate_var_dict(stripes['blue'], store_timepoint_N, 0.4)
    
    stripes_var['blue'] = blue_var
    animation_type['blue'] = 'connect'
    connect_var['blue'] = connections_blue
    
    print("Done BLUE")

##################################################################################################
""" PURPLE """
# random ICs, increasing connections and diffusing calcium, heads towards Turing pattern
if 'purple' in run_selection:
    
    # create connections and write to file
    # connections_purple = create_connections(stripes['purple'], store_timepoint_N, connection_params)
    # dir_name = 'connections'
    # if not path.isdir(dir_name):
    #     mkdir(dir_name)
    # write_connections_to_file(connections_purple, path.join(dir_name, 'purple_'))
    
    # load connections from file
    connections_purple = load_connections_from_file(stripes['purple'], path.join('connections', 'purple_'))
    
    initial_connections_purple = connections_purple['initial_connections']
    birth_connections_purple = connections_purple['birth_connections']
    death_connections_purple = connections_purple['death_connections']
    
    purple_var = allocate_var_dict(stripes['purple'], store_timepoint_N, 0.6)
    stripes_var['purple'] = purple_var
    # purple_var = set_var_dict_to_random_val_in_range_all_t(purple_var, [0,1])
    
    # initialize V_PLC
    params["V_PLC"] = allocate_var_dict(stripes['purple'], 1, 0.787)

    # set cell-cell communication, 0 => OFF, standard 0.02
    params["D_IP3"] = 0.02

    # allocation initial conditions for variables
    Ca_cyt_0 = 2
    Ca_cyt = allocate_var_dict(stripes['purple'], store_timepoint_N, Ca_cyt_0)
    Ca_stored = allocate_var_dict(stripes['purple'], store_timepoint_N, params["c_tot"] - Ca_cyt_0)
    ip3 = allocate_var_dict(stripes['purple'], store_timepoint_N, 0.2)
    ip3R_act = allocate_var_dict(stripes['purple'], store_timepoint_N, 0.6)

    variables = Ca_cyt, ip3, Ca_stored, ip3R_act,

    # set ICs randomly from V_PLC 0.787 df
    ICs_df = pd.read_csv("ICs.csv", delimiter=',', header=0, index_col=0)
    # create_less_random_initial_conditions_df(stripes['purple'], ICs_df, 'not_random_ICs_purple.csv')
    less_random_ICs_df = pd.read_csv("not_random_ICs_purple.csv", delimiter=',', header=0, index_col=0)
    variables = set_initial_conditions_from_df_less_random(less_random_ICs_df, variables)
    
    # variables = set_initial_conditions_from_df(ICs_df, variables) # random initial conditions
    # variables = set_constant_initial_conditions_from_df_first_row(ICs_df, variables) # constant initial conditions (in space)
    
    Ca_cyt_new, ip3_new, Ca_stored_new, ip3R_act_new, = Ca_cyt, ip3, Ca_stored, ip3R_act,
    Ca_cyt_new, ip3_new, Ca_stored_new, ip3R_act_new, = politi_reduced_connectivity(variables, run_t, store_t, stripes['purple'], initial_connections_purple, birth_connections_purple, death_connections_purple, params)

    stripes_var['purple'] = Ca_cyt_new
    
    animation_type['purple'] = 'colour'
    
    # animation_type['purple'] = 'connect'
    # connect_var['purple'] = connections_purple
    
    print("Done PURPLE")

##################################################################################################
""" PLOT """

flag_colours = ''
for stripe in run_selection:
    flag_colours = flag_colours + '_' + stripe
flag_name = 'flag' + flag_colours
    
plot_var_flag(run_selection, stripes_var, connect_var, 150, stripes, animation_type, (hex_x_N, hex_y_N), pointy, 12, save_dir + flag_name + '_initial')

# animate_var_flag(run_selection, stripes_var, connect_var, store_timepoint_N, stripes, animation_type, (hex_x_N,hex_y_N), pointy, 12, save_dir + flag_name)

##################################################################################################
