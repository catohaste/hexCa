import numpy as np

from lib import *

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
    for hexa in hex_array:
        var_dict[hexa] = value * np.ones((timepoint_N,), dtype=float)
            
    return var_dict
    
def initialize_leftmost_hexes_to_value(var_dict, hex_array, value, pointy_layout):
    
    radius = pointy_layout.size[0]
    
    centers = {}
    for hexa in hex_array:
        centers[hexa] = hex_to_pixel(pointy_layout, hexa)
        
    x_coords = [centers[hexa][0] for hexa in centers]
    rows = 3
    rows_edge = 4
    x_min_lim = min(x_coords) + radius*rows_edge
    x_max_lim = max(x_coords) - radius*rows_edge
    frac = 0.5
    x_middle = min(x_coords) + frac * (max(x_coords) - min(x_coords)) 
    x_middle_up = x_middle + radius*rows
    x_middle_down = x_middle - radius*rows
    
    middle_hexes = [hexa for hexa in centers if centers[hexa][0] > x_middle_down and centers[hexa][0] < x_middle_up]    
    right_hexes = [hexa for hexa in centers if centers[hexa][0] > x_max_lim]
    left_hexes = [hexa for hexa in centers if centers[hexa][0] < x_min_lim]
    for hexa in middle_hexes:
        var_dict[hexa][0] = value
    for hexa in left_hexes:
        var_dict[hexa][0] = value
    for hexa in right_hexes:
        var_dict[hexa][0] = value
    
    return var_dict

def create_val_loc_tuple_std_layout(var_dict, hexes, pointy_layout):
    
    value_loc_tuples = []
    for hexa in hexes:
        value_loc_tuples.append((var_dict[hexa], hex_to_pixel(pointy_layout, hexa)))
        
    return value_loc_tuples
    

    