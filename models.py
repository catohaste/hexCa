import numpy as np
import random

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
        
    value_increment = 0.2
    for t in range(1, timepoint_N):
        for hexa in hex_array:
            if var_dict[hexa][t-1] <= value_range[0] or var_dict[hexa][t-1] >= value_range[1]:
                direction_dict[hexa] *= -1     
            var_dict[hexa][t] = var_dict[hexa][t - 1] + (value_increment * direction_dict[hexa])
    
    return var_dict
    
    
    