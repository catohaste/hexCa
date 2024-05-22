import numpy as np
from collections import Counter
from copy import deepcopy
import random
import networkx as nx
import pickle
import json
from os import path
import pandas as pd

from lib import *

def Hill(a,K,var,pow):
    """Hill equation"""
    return (a * ( (var ** pow) / (K ** pow + var ** pow) ))
    
def tau_p(k_3K, k_5P):
    return 1 / (k_3K + k_5P)

def get_lim_of_var_dict(var_dict):
    
    min_val = 0
    max_val = 0
    
    for hexa in var_dict:
        hex_min = min(var_dict[hexa])
        hex_max = max(var_dict[hexa])
        min_val = min(min_val, hex_min)
        max_val = max(max_val, hex_max)
        
        return min_val, max_val

def create_layout_from_dict(layout_dict):
    
    if layout_dict["layout_str"] == "pointy":
        orientation = layout_pointy
    elif layout_dict["layout_str"] == "flat":
        orientation = layout_flat
    else:
        print("layout_dict['layout_str'] should be 'pointy' or 'flat'.\nRevert to 'pointy' as default." )
        orientation = layout_pointy
        
    layout = Layout(orientation, Point(layout_dict["radius"], layout_dict["radius"]), Point(layout_dict["center_x"], layout_dict["center_y"]))
    
    return layout
    
def allocate_var_dict(hex_array, timepoint_N, value):
    
    var_dict = {}
    if timepoint_N == 1:
        for hexa in hex_array:
            var_dict[hexa] = value
    else:
        for hexa in hex_array:
            var_dict[hexa] = value * np.ones((timepoint_N,), dtype=float)
            
    return var_dict
    
def initialize_column_of_hexes_to_value(var_dict, hex_array, value, x_coord_frac, half_n_cols, pointy_layout):
    """
    Initialize (t=0 only) n_cols of hexs to given value.
    Give the x_coord of hexs to initialize as fraction of whole (i.e. between 0,1)
    Strictly speaking, layout shouldn't be require. But makes implementation here somewhat easier.
    """
    
    radius = pointy_layout.size[0]
    
    centers = {}
    for hexa in hex_array:
        centers[hexa] = hex_to_pixel(pointy_layout, hexa)
        
    x_coords = [centers[hexa][0] for hexa in centers]

    x_middle = min(x_coords) + x_coord_frac * (max(x_coords) - min(x_coords)) 
    x_middle_up = x_middle + radius*half_n_cols
    x_middle_down = x_middle - radius*half_n_cols
    
    middle_hexes = [hexa for hexa in centers if centers[hexa][0] > x_middle_down and centers[hexa][0] < x_middle_up]
    for hexa in middle_hexes:
        var_dict[hexa][0] = value
    
    return var_dict
    
## FIX ME, there don't need to be 2 of these functions
def initialize_column_of_hexes_to_value_2(var_dict, hex_array, value, x_coord_frac, half_n_cols, pointy_layout):
    """
    Initialize (t=0 only) n_cols of hexs to given value.
    Give the x_coord of hexs to initialize as fraction of whole (i.e. between 0,1)
    Strictly speaking, layout shouldn't be require. But makes implementation here somewhat easier.
    """
    
    radius = pointy_layout.size[0]
    
    centers = {}
    for hexa in hex_array:
        centers[hexa] = hex_to_pixel(pointy_layout, hexa)
        
    x_coords = [centers[hexa][0] for hexa in centers]

    x_middle = min(x_coords) + x_coord_frac * (max(x_coords) - min(x_coords)) 
    x_middle_up = x_middle + radius*half_n_cols
    x_middle_down = x_middle - radius*half_n_cols
    
    middle_hexes = [hexa for hexa in centers if centers[hexa][0] > x_middle_down and centers[hexa][0] < x_middle_up]
    for hexa in middle_hexes:
        var_dict[hexa] = value
    
    return var_dict

def initialize_var_dict_to_x_gradient(var_dict, hex_array, value_range, pointy_layout):
    
    radius = pointy_layout.size[0]
    
    centers = {}
    for hexa in hex_array:
        centers[hexa] = hex_to_pixel(pointy_layout, hexa)
        
    x_coords = [centers[hexa][0] for hexa in centers]
    
    x_cols = list(set(x_coords))
    x_cols.sort()
    x_cols_N = len(x_cols)
    
    x_vals = np.linspace(value_range[0], value_range[1], x_cols_N)
    x_col_val_zipped = zip(x_cols, x_vals)
    
    for hexa in hex_array:
        current_x = centers[hexa][0]
        x_idx = x_cols.index(current_x)
        var_dict[hexa] = x_vals[x_idx]
    
    return var_dict
    
def initialize_var_dict_to_x_gradient_from_ICs_df(df, variables, hex_array, pointy_layout):
    
    Ca_cyt, ip3, Ca_stored, ip3R_act, = variables
    
    radius = pointy_layout.size[0]
    
    centers = {}
    for hexa in hex_array:
        centers[hexa] = hex_to_pixel(pointy_layout, hexa)
        
    x_coords = [centers[hexa][0] for hexa in centers]
    
    x_cols = list(set(x_coords))
    x_cols.sort()
    x_cols_N = len(x_cols)
    
    row_indices = list(df.index)
    rows_N = len(df)
    
    if x_cols_N > rows_N:
        row_copies_N = int(np.ceil(x_cols_N / rows_N))
        out_row_indices = []
        for i in range(row_copies_N):
            for index in row_indices:
                out_row_indices.append(index)
        out_row_indices.sort(reverse=True) # reverse sort means stripe moves left to right
        out_row_indices = out_row_indices[:x_cols_N]
    else:
        out_row_indices = row_indices[:x_cols_N]
        out_row_indices.sort(reverse=True) # reverse sort means stripe moves left to right
    
    for hexa in hex_array:
        current_x = centers[hexa][0]
        x_idx = x_cols.index(current_x)
        row_idx = out_row_indices[x_idx]
        row = df.loc[row_idx]
        Ca_cyt[hexa][0] = row['Ca_cyt']
        ip3[hexa][0] = row['IP3']
        Ca_stored[hexa][0] = row['Ca_stored']
        ip3R_act[hexa][0] = row['IP3R_act']
    
    return Ca_cyt, ip3, Ca_stored, ip3R_act,

def initialize_var_dict_to_random_val_in_range(var_dict, value_range):
    
    for hexa in var_dict:
        var_dict[hexa][0] = np.random.uniform(low=value_range[0], high=value_range[1])
    
    return var_dict
    
def set_var_dict_to_random_val_in_range_all_t(var_dict, value_range):
    
    keys = []
    for key in var_dict:
        keys.append(key)
        
    timepoint_N = len(var_dict[keys[0]])

    for hexa in var_dict:
        var_dict[hexa][0] = np.random.uniform(low=value_range[0], high=value_range[1])
        for t in range(1, timepoint_N):
            var_dict[hexa][t] = var_dict[hexa][0]

    return var_dict
    
def set_initial_conditions_from_df(df, variables):
    
    Ca_cyt, ip3, Ca_stored, ip3R_act, = variables
    
    for hexa in Ca_cyt:
        row = df.sample(1)
        Ca_cyt[hexa][0] = row['Ca_cyt'].iloc[0]
        ip3[hexa][0] = row['IP3'].iloc[0]
        Ca_stored[hexa][0] = row['Ca_stored'].iloc[0]
        ip3R_act[hexa][0] = row['IP3R_act'].iloc[0]
    
    return Ca_cyt, ip3, Ca_stored, ip3R_act,
    
def create_less_random_initial_conditions_df(hexes_list, ICs_df, filename):
    
    existing_column_names = ICs_df.columns.values.tolist()
    qrs = ['q','r','s']
    new_column_names = qrs + existing_column_names
    less_random_ICs_df = pd.DataFrame(columns=new_column_names)
    
    for idx, hexa in enumerate(hexes_list):
        row = ICs_df.sample(1, ignore_index=True)
        less_random_ICs_df.loc[idx] = [hexa.q, hexa.r, hexa.s] + list(row.loc[0])
    
    less_random_ICs_df.to_csv(path_or_buf=filename)
    
def set_initial_conditions_from_df_less_random(df, variables):
    """ instead of choosing a random row (or set of initial conditions) collect a number for each hex."""
    
    Ca_cyt, ip3, Ca_stored, ip3R_act, = variables
    
    for index, row in df.iterrows():
        hexa = Hex(row['q'],row['r'],row['s'])
        Ca_cyt[hexa][0] = row['Ca_cyt']
        ip3[hexa][0] = row['IP3']
        Ca_stored[hexa][0] = row['Ca_stored']
        ip3R_act[hexa][0] = row['IP3R_act']
    
    return Ca_cyt, ip3, Ca_stored, ip3R_act,
  
def set_constant_initial_conditions_from_df_first_row(df, variables):
    
    Ca_cyt, ip3, Ca_stored, ip3R_act, = variables
    
    row = df.iloc[[0]]
    
    for hexa in Ca_cyt:
        Ca_cyt[hexa][0] = row['Ca_cyt'].iloc[0]
        ip3[hexa][0] = row['IP3'].iloc[0]
        Ca_stored[hexa][0] = row['Ca_stored'].iloc[0]
        ip3R_act[hexa][0] = row['IP3R_act'].iloc[0]
    
    return Ca_cyt, ip3, Ca_stored, ip3R_act,
      
def make_links(var_dict, store_t):
    
    links = []
    
    var_dict_lim = get_lim_of_var_dict(var_dict)
    var_dict_diff = var_dict_lim[1] - var_dict_lim[0]
    on_threshold = var_dict_lim[0] + 0.6 * var_dict_diff
    close_threhold = 0.05 * var_dict_diff
    
    t_diff_range = (1, 1)
    t_diff_points = list(range(t_diff_range[0], t_diff_range[1] + 1))
    t_diff_min = min(t_diff_points)
    
    for t_idx, time in enumerate(store_t):
        links_current_time = []
        if t_idx < t_diff_min:
            # links.append(links_current_time)
            continue
        for hexa in var_dict:
            for direction in range(6):
                hexa_neighbor = hex_neighbor(hexa, direction)
                for t_diff in t_diff_points:
                    try:
                        if var_dict[hexa][t_idx] > on_threshold and abs(var_dict[hexa][t_idx] - var_dict[hexa_neighbor][t_idx - t_diff]) < close_threhold:
                            link = {
                                "time_lim_idx": (t_idx - t_diff, t_idx),
                                "hexes": (hexa_neighbor, hexa),
                                "var_lims": (var_dict[hexa_neighbor][t_idx - t_diff], var_dict[hexa][t_idx])
                            }
                            links_current_time.append(link)
                    except KeyError:
                        continue
        links.append(links_current_time)

    return links

def create_val_loc_dict_average_over_y(value_loc_tuples):
    
    value_loc_dict_averaged_over_y = {}
    
    x_locs = [np.round(loc[0], decimals=2) for val,loc in value_loc_tuples]
    x_loc_counter_dict = Counter(x_locs)
    
    timepoint_N = len(value_loc_tuples[0][0])
    
    for x_loc in Counter(x_locs).keys():
        value_loc_dict_averaged_over_y[x_loc] = np.zeros(timepoint_N, )

    for val, loc in value_loc_tuples:
        value_loc_dict_averaged_over_y[np.round(loc[0], decimals=2)] = np.add(value_loc_dict_averaged_over_y[np.round(loc[0], decimals=2)] , val)

    for loc in value_loc_dict_averaged_over_y:
        value_loc_dict_averaged_over_y[loc] = value_loc_dict_averaged_over_y[loc] / x_loc_counter_dict[loc]
        
    return value_loc_dict_averaged_over_y
    
def create_val_loc_tuple_std_layout(var_dict, hexes, pointy_layout):
    
    value_loc_tuples = []
    for hexa in hexes:
        value_loc_tuples.append((var_dict[hexa], hex_to_pixel(pointy_layout, hexa)))
        
    return value_loc_tuples
    
def unpack_val_loc_tuple_std_layout(value_loc_tuples, pointy):
    hex_array = []
    var_dict = {}
    for val,loc in value_loc_tuples:
        hexa = pixel_to_hex(pointy, loc)
        hex_array.append(hexa)
        var_dict[hexa] = val
        
    return hex_array, var_dict,
    
def hex_to_tuple(hex):
    return (hex.q,hex.r,hex.s)
    
def tuple_to_hex(hex_tuple):
    return Hex(hex_tuple[0], hex_tuple[1], hex_tuple[2])
    
def get_mean_degree_fraction(connections, potential_connections):
    
    degree_ratio = {}
    
    for hexa in connections:
         degree_ratio[hexa] = connections.degree(hexa) / potential_connections.degree(hexa)

    mean_degree_ratio = np.mean(list(degree_ratio.values()))
    
    return mean_degree_ratio
        
def get_current_connections(current_t, initial_connections, birth_connections, death_connections):
    
    current_connections = deepcopy(initial_connections)
    
    current_t_list = []
    for t in birth_connections:
        if t < current_t:
            current_t_list.append(t)
    for t in death_connections:
        if t < current_t and t not in current_t_list:
            current_t_list.append(t)
            
    current_t_list.sort()
    for t in current_t_list:
        try:
            new_edges = birth_connections[t]
        except KeyError:
            new_edges = []
        try:
            remove_edges = death_connections[t]
        except KeyError:
            remove_edges = []
        for new_edge in new_edges:
            current_connections.add_edge(new_edge[0], new_edge[1])
        for remove_edge in remove_edges:
            current_connections.remove_edge(remove_edge[0], remove_edge[1])
        
    return current_connections
    
def get_birth_and_death_t(timepoint_N, connection_params):
    
    birth_t = range(connection_params['birth_connect_dt'], timepoint_N, connection_params['birth_connect_dt'])
    death_t = range(connection_params['death_connect_dt'], timepoint_N, connection_params['death_connect_dt'])
    
    birth_and_death_t = list(set(list(birth_t) + list(death_t)))
    birth_and_death_t.sort()
    
    return birth_t, death_t, birth_and_death_t

def create_connections(hexes_list, timepoint_N, connection_params):
    
    """potential"""
    potential_connections = nx.Graph()
    potential_connections.add_nodes_from(hexes_list)
    for hexa in hexes_list:
        neighbors = hex_neighbors_cumulative_distance(hexa, connection_params['neighbour_dist_limit'])
        for neighbor in neighbors:
            if neighbor in hexes_list:
                potential_connections.add_edge(hexa, neighbor)
                
    # print(potential_connections)

    potential_deg = list(dict(potential_connections.degree()).values())
    # print(min(potential_deg), max(potential_deg), np.mean(potential_deg))
    # print("potential", potential_connections)

    """initial"""
    if connection_params['init_avg_degree_fraction'] == 1:
        initial_connections = potential_connections
        birth_connections = {}
        death_connections = {}
    else:
        initial_connections = nx.Graph()
        initial_connections.add_nodes_from(hexes_list)

        # print('before', get_mean_degree_fraction(initial_connections, potential_connections))
        while get_mean_degree_fraction(initial_connections, potential_connections) < connection_params['init_avg_degree_fraction']:
            remaining_possible_connections = nx.difference(potential_connections, initial_connections)
            random_edge = random.sample(list(remaining_possible_connections.edges), 1)[0]
            initial_connections.add_edge(random_edge[0], random_edge[1])
        # print('after',get_mean_degree_fraction(initial_connections, potential_connections))

        """birth and death"""
        birth_connections = {}
        death_connections = {}
        
        birth_t, death_t, birth_and_death_t = get_birth_and_death_t(timepoint_N, connection_params)

        # allocate
        # birth_t = range(connection_params['birth_connect_dt'], timepoint_N, connection_params['birth_connect_dt'])
        for t in birth_t:
            birth_connections[t] = []
        # death_t = range(connection_params['death_connect_dt'], timepoint_N, connection_params['death_connect_dt'])
        for t in death_t:
            death_connections[t] = []
    
        # set up connections
        # birth_and_death_t = list(set(list(birth_t) + list(death_t)))
        # birth_and_death_t.sort()
        for current_t in birth_and_death_t:
            current_connections = get_current_connections(current_t, initial_connections, birth_connections, death_connections)
            possible_birth_connections = nx.difference(potential_connections, current_connections)
            if len(list(possible_birth_connections.edges)) > 0:
                if current_t in birth_t:
                    random_edge = random.sample(list(possible_birth_connections.edges), 1)[0]
                    birth_connections[current_t].append(random_edge)
                if current_t in death_t:
                    random_edge = random.sample(list(current_connections.edges), 1)[0]
                    death_connections[current_t].append(random_edge)
            else:
                print('Fully connected at timepoint ' + str(current_t))
                break
                
        connections_dict = {
            'initial_connections'   : initial_connections,
            'birth_connections'     : birth_connections,
            'death_connections'     : death_connections
        }
        
        return connections_dict
        
def write_connections_to_file(connections_dict, filename_prefix):
    
    # you cannot pickle Hex class from lib.py
    # I've done a ridiculous work around so that you can store Hex
    # Hex is just named tuple from 'collections'. I have removed the name
    
    ### INITIAL #####
    initial_connections = connections_dict['initial_connections']
    initial_connections_edgelist = list(initial_connections.edges)
    
    non_hex_initial_edgelist = []
    for hex_edge in initial_connections_edgelist:
        non_hex_edge = hex_edge_to_non_hex_edge(hex_edge)
        non_hex_initial_edgelist.append(non_hex_edge)
    
    with open(filename_prefix + 'initial.p', 'wb') as initial_file:
        pickle.dump(non_hex_initial_edgelist, initial_file)
        
    ### BIRTH #####
    birth_connections = connections_dict['birth_connections']
    birth_connections_non_hex = {}
    for key in birth_connections:
        edge_list = birth_connections[key]
        birth_connections_non_hex[key] = [hex_edge_to_non_hex_edge(hex_edge) for hex_edge in edge_list]
        
    with open(filename_prefix + 'birth.p', 'wb') as birth_file:
        pickle.dump(birth_connections_non_hex, birth_file)
        
    ### DEATH #####
    death_connections = connections_dict['death_connections']
    death_connections_non_hex = {}
    for key in death_connections:
        edge_list = death_connections[key]
        death_connections_non_hex[key] = [hex_edge_to_non_hex_edge(hex_edge) for hex_edge in edge_list]

    with open(filename_prefix + 'death.p', 'wb') as death_file:
        pickle.dump(death_connections_non_hex, death_file)
    
    return
    
def load_connections_from_file(hexes_list, filename_prefix):
    
    ### INITIAL #####
    with open(filename_prefix + 'initial.p', 'rb') as initial_file:   # Unpickling
        non_hex_initial_edgelist = pickle.load(initial_file)
    
    initial_connections = nx.Graph()
    initial_connections.add_nodes_from(hexes_list)

    for non_hex_edge in non_hex_initial_edgelist:
        hex_edge = non_hex_edge_to_hex_edge(non_hex_edge)
        initial_connections.add_edge(hex_edge[0], hex_edge[1])
        
    ### BIRTH #####
    with open(filename_prefix + 'birth.p', 'rb') as birth_file:
        birth_connections_non_hex = pickle.load(birth_file)
    birth_connections = {}
    for key in birth_connections_non_hex:
        edge_list = birth_connections_non_hex[key]
        birth_connections[key] = [non_hex_edge_to_hex_edge(non_hex_edge) for non_hex_edge in edge_list]
        
    ### DEATH #####
    with open(filename_prefix + 'death.p', 'rb') as death_file:
        death_connections_non_hex = pickle.load(death_file)
    death_connections = {}
    for key in death_connections_non_hex:
        edge_list = death_connections_non_hex[key]
        death_connections[key] = [non_hex_edge_to_hex_edge(non_hex_edge) for non_hex_edge in edge_list]

    connections_dict = {}
    connections_dict['initial_connections'] = initial_connections
    connections_dict['birth_connections'] = birth_connections
    connections_dict['death_connections'] = death_connections
    
    return connections_dict
    
def hex_edge_to_non_hex_edge(hex_edge):
    
    non_hex_edge = ( (hex_edge[0].q, hex_edge[0].r, hex_edge[0].s) ,  (hex_edge[1].q, hex_edge[1].r, hex_edge[1].s) )
    
    return non_hex_edge
    
def non_hex_edge_to_hex_edge(non_hex_edge):
    
    hex_edge = ( Hex(non_hex_edge[0][0], non_hex_edge[0][1], non_hex_edge[0][2]), Hex(non_hex_edge[1][0], non_hex_edge[1][1], non_hex_edge[1][2]) )
    
    return hex_edge

def find_hexes_min_max(hexes_list):
    
    # [ min[q,r,s] , max[q,r,s] ]
    qrs_min_max = [ [float('inf'),float('inf'),float('inf')] , [float('-inf'),float('-inf'),float('-inf')] ]
    for hexa in hexes_list:
        if hexa.q < qrs_min_max[0][0]:
            qrs_min_max[0][0] = hexa.q
        if hexa.q > qrs_min_max[1][0]:
            qrs_min_max[1][0] = hexa.q
            
        if hexa.r < qrs_min_max[0][1]:
            qrs_min_max[0][1] = hexa.q
        if hexa.r > qrs_min_max[1][1]:
            qrs_min_max[1][1] = hexa.q
            
        if hexa.s < qrs_min_max[0][2]:
            qrs_min_max[0][2] = hexa.q
        if hexa.s > qrs_min_max[1][2]:
            qrs_min_max[1][2] = hexa.q
            
    return qrs_min_max

    