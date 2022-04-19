import numpy as np
import random

from lib import *
from functions import *

def politi(variables, dt, run_timepoint_N, hex_array, params):
    """
    variables = Ca_cyt, ip3, Ca_stored, ip3R_act,
    """
    # upack params
    K_3K = params["K_3K"]
    k_3K = params["k_3K"]
    k_5P = params["k_5P"]
    K_PLC = params["K_PLC"]
    beta = params["beta"]
    V_SERCA = params["V_SERCA"]
    K_SERCA = params["K_SERCA"]
    V_pm = params["V_pm"]
    K_pm = params["K_pm"]
    v_0 = params["v_0"]
    phi = params["phi"]
    epsilon = params["epsilon"]
    k_1 = params["k_1"]
    k_2 = params["k_2"]
    K_a = params["K_a"]
    K_i = params["K_i"]
    K_p = params["K_p"]
    tau_r = params["tau_r"]
    
    V_PLC = params["V_PLC"]
    
    D_IP3 = params["D_IP3"]
    
    ##################################################################
    
    Ca_cyt, ip3, Ca_stored, ip3R_act, = variables
    
    ##################################################################
    
    tau_p = 1 / (k_3K + k_5P)
    eta = k_3K * tau_p
    
    # new_Ca_cyt = allocate_var_dict(hex_array, 1, 0)
    # new_ip3 = allocate_var_dict(hex_array, 1, 0)
    # new_Ca_stored = allocate_var_dict(hex_array, 1, 0)
    # new_ip3R_act = allocate_var_dict(hex_array, 1, 0)
    
    for t_idx in range(1, run_timepoint_N):
        
        Ca_temp = allocate_var_dict(hex_array, 1, 0)
        V_PLC_scaled = allocate_var_dict(hex_array, 1, 0)
    
        for hexa in hex_array:
        
            V_PLC_scaled[hexa] = V_PLC[hexa] * tau_p
        
            Ca_temp[hexa] = (k_1 * ((Hill(ip3R_act[hexa][t_idx-1], K_a, Ca_cyt[hexa][t_idx-1], 1 ) * Hill(1, K_p , ip3[hexa][t_idx-1], 1))**3) + k_2) * (Ca_stored[hexa][t_idx-1] - Ca_cyt[hexa][t_idx-1]) - Hill(V_SERCA, K_SERCA, Ca_cyt[hexa][t_idx-1], 2)
        
            ip3R_act[hexa][t_idx] = ip3R_act[hexa][t_idx-1] + dt * ( (1 / tau_r) * (1 - ip3R_act[hexa][t_idx-1] * ((K_i + Ca_cyt[hexa][t_idx-1]) / K_i) ) )
        
            Ca_stored[hexa][t_idx] = Ca_stored[hexa][t_idx-1] + dt * (((-1) * Ca_temp[hexa]) * (1 / beta) )
        
            Ca_cyt[hexa][t_idx] = Ca_cyt[hexa][t_idx-2] + dt * ( Ca_temp[hexa] + epsilon * (v_0 + phi * V_PLC_scaled[hexa] - Hill(V_pm, K_pm, Ca_cyt[hexa][t_idx-1], 2)) )
        
            # only IP3 travels between cells
            # this calculates an ip3 neighborhood average with no-flux boundary conditions
            ip3_neighborhood_sum = 0
            ip3_neighbor_counter = 0
            for direction in range(6):
                neighbor = hex_neighbor(hexa, direction)
                try:
                    ip3_neighborhood_sum += ip3[neighbor][t_idx-1]
                    ip3_neighbor_counter += 1
                except KeyError:
                    continue
            ip3_neighborhood_avg = ip3_neighborhood_sum/ip3_neighbor_counter
        
            ip3[hexa][t_idx] = ip3[hexa][t_idx-1] + dt * ( (1 / tau_p) * (Hill(V_PLC_scaled[hexa], K_PLC, Ca_cyt[hexa][t_idx-1], 2) - ( ip3[hexa][t_idx-1] * ( Hill(eta, K_3K, Ca_cyt[hexa][t_idx-1], 2) + 1 - eta) ) ) + D_IP3 * (ip3_neighborhood_avg - ip3[hexa][t_idx-1])
    
    return Ca_cyt, ip3, Ca_stored, ip3R_act,

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
    
    
    