import numpy as np

from lib import *

def Hill(a,K,var,pow):
    """Hill equation"""
    return (a * ( (var ** pow) / (K ** pow + var ** pow) ))
    
def tau_p(k_3K, k_5P):
    return 1 / (k_3K + k_5P)

def create_layout_from_dict(layout_dict):
    
    if layout_dict["layout_str"] is "pointy":
        orientation = layout_pointy
    elif layout_dict["layout_str"] is "flat":
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

    