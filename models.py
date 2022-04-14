import numpy as np
import random

from lib import *
from functions import *

def diffusion(var_dict, hex_array, timepoint_N):
    """
    zero-flux boundary conditions
    the try function is a weird fix, but it works
    """
    
    diff_const = 0.01
    for t_idx in range(1, timepoint_N):
        for hexa in hex_array:
            neighborhood_sum = 0
            neighbor_counter = 0
            for direction in range(6):
                neighbor = hex_neighbor(hexa, direction)
                try:
                    neighborhood_sum += var_dict[neighbor][t_idx - 1]
                    neighbor_counter += 1
                except KeyError:
                    continue
            neighborhood_avg = neighborhood_sum/neighbor_counter
            var_dict[hexa][t_idx] = var_dict[hexa][t_idx-1] + diff_const * (neighborhood_avg - var_dict[hexa][t_idx - 1])
            
    return var_dict


def random_pulsing(var_dict, hex_array, timepoint_N, value_range):
    """
    First timepoint, generate a random value within a given range, and direction (+/-).
    Subsequent timepoints, move value in appropriate direction.
    If max or min is reached, change direction.
    """
    
    # initialize data
    direction_dict = {}
    directions = [1,-1]
    for hexa in hex_array:
        var_dict[hexa] = np.ndarray((timepoint_N,), dtype=float)
        var_dict[hexa][0] = random.randint(value_range[0], value_range[1])
        direction_dict[hexa] = random.choice(directions)
        
    
    value_increment = (2 * (value_range[1] - value_range[0])) / (timepoint_N - 1) # for looping purposes
    for t in range(1, timepoint_N):
        for hexa in hex_array:
            if var_dict[hexa][t-1] <= value_range[0] or var_dict[hexa][t-1] >= value_range[1]:
                direction_dict[hexa] *= -1
            var_dict[hexa][t] = var_dict[hexa][t - 1] + (value_increment * direction_dict[hexa])
    
    return var_dict
    
    
    