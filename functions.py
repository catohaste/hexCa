import numpy as np
from collections import Counter
from copy import deepcopy

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
        var_dict[hexa][0] = x_vals[x_idx]
    
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
        out_row_indices.sort()
        out_row_indices = out_row_indices[:x_cols_N]
    else:
        out_row_indices = row_indices[:x_cols_N]
    
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
        Ca_cyt[hexa][0] = row['Ca_cyt']
        ip3[hexa][0] = row['IP3']
        Ca_stored[hexa][0] = row['Ca_stored']
        ip3R_act[hexa][0] = row['IP3R_act']
    
    return Ca_cyt, ip3, Ca_stored, ip3R_act,
    
def set_initial_conditions_from_df_less_random(df, variables):
    """ instead of choosing a random row (or set of initial conditions) collect a number for each hex.
    """
    
    Ca_cyt, ip3, Ca_stored, ip3R_act, = variables
    
    for index, row in df.iterrows():
        # print(row['q'], row['r'],row['s'])
        hexa = Hex(row['q'],row['r'],row['s'] )
        Ca_cyt[hexa][0] = row['Ca_cyt']
        ip3[hexa][0] = row['IP3']
        Ca_stored[hexa][0] = row['Ca_stored']
        ip3R_act[hexa][0] = row['IP3R_act']
    
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
    