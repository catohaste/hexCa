import numpy as np
import random

from lib import *
from functions import *

def politi_reduced_connectivity(variables, connections_over_t, run_t, store_t, hex_array, params):
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
    
    old_Ca_cyt = allocate_var_dict(hex_array, 1, 0)
    old_ip3 = allocate_var_dict(hex_array, 1, 0)
    old_Ca_stored = allocate_var_dict(hex_array, 1, 0)
    old_ip3R_act = allocate_var_dict(hex_array, 1, 0)
    
    for hexa in hex_array: 
        old_Ca_cyt[hexa] = Ca_cyt[hexa][0]
        old_ip3[hexa] = ip3[hexa][0]
        old_Ca_stored[hexa] = Ca_stored[hexa][0]
        old_ip3R_act[hexa] = ip3R_act[hexa][0]
        
    new_Ca_cyt = allocate_var_dict(hex_array, 1, 0)
    new_ip3 = allocate_var_dict(hex_array, 1, 0)
    new_Ca_stored = allocate_var_dict(hex_array, 1, 0)
    new_ip3R_act = allocate_var_dict(hex_array, 1, 0)
    
    ##################################################################
    
    tau_p = 1 / (k_3K + k_5P)
    eta = k_3K * tau_p
    
    Ca_temp = allocate_var_dict(hex_array, 1, 0)
    V_PLC_scaled = allocate_var_dict(hex_array, 1, 0)
    
    dt = run_t[1] - run_t[0]
    store_dt = store_t[1] - store_t[0]
    time_scaling = int(store_dt / dt)
    
    current_connection_graph = connections_over_t[0]
    
    for t_idx, t in enumerate(run_t[1:]):
    
        for hexa in hex_array:
        
            V_PLC_scaled[hexa] = V_PLC[hexa] * tau_p
        
            Ca_temp[hexa] = (k_1 * ((Hill(old_ip3R_act[hexa], K_a, old_Ca_cyt[hexa], 1 ) * Hill(1, K_p , old_ip3[hexa], 1))**3) + k_2) * (old_Ca_stored[hexa] - old_Ca_cyt[hexa]) - Hill(V_SERCA, K_SERCA, old_Ca_cyt[hexa], 2)
        
            new_ip3R_act[hexa] = old_ip3R_act[hexa] + dt * ( (1 / tau_r) * (1 - old_ip3R_act[hexa] * ((K_i + old_Ca_cyt[hexa]) / K_i) ) )
        
            new_Ca_stored[hexa] = old_Ca_stored[hexa] + dt * (((-1) * Ca_temp[hexa]) * (1 / beta) )
        
            new_Ca_cyt[hexa] = old_Ca_cyt[hexa] + dt * ( Ca_temp[hexa] + epsilon * (v_0 + phi * V_PLC_scaled[hexa] - Hill(V_pm, K_pm, old_Ca_cyt[hexa], 2)) )
        
            # only IP3 travels between cells
            # this calculates an ip3 neighborhood average with no-flux boundary conditions
            if len(current_connection_graph.edges(hexa)) == 0:
                
                new_ip3[hexa] = old_ip3[hexa] + dt * ( (1 / tau_p) * (Hill(V_PLC_scaled[hexa], K_PLC, old_Ca_cyt[hexa], 2) - ( old_ip3[hexa] * ( Hill(eta, K_3K, old_Ca_cyt[hexa], 2) + 1 - eta) ) ) )
                
            else:
            
                ip3_neighborhood_sum = 0
                ip3_neighbor_counter = 0
            
                for edge in current_connection_graph.edges(hexa):
                    if edge[0] == hexa:
                        neighbor = edge[1]
                    else:
                        neighbor = edge[0]
                    
                    ip3_neighborhood_sum += old_ip3[neighbor]
                    ip3_neighbor_counter += 1
                    
                if ip3_neighbor_counter != current_connection_graph.degree(hexa):
                    print("Something has gone awry")
                    
                ip3_neighborhood_avg = ip3_neighborhood_sum / ip3_neighbor_counter
        
                new_ip3[hexa] = old_ip3[hexa] + dt * ( (1 / tau_p) * (Hill(V_PLC_scaled[hexa], K_PLC, old_Ca_cyt[hexa], 2) - ( old_ip3[hexa] * ( Hill(eta, K_3K, old_Ca_cyt[hexa], 2) + 1 - eta) ) ) + D_IP3 * (ip3_neighborhood_avg - old_ip3[hexa]) )
            
            old_Ca_cyt[hexa] = new_Ca_cyt[hexa]
            old_ip3[hexa] = new_ip3[hexa]
            old_Ca_stored[hexa] = new_Ca_stored[hexa]
            old_ip3R_act[hexa] = new_ip3R_act[hexa]
            
            if t in store_t:
                
                t_idx = store_t.tolist().index(t)
                
                current_connection_graph = connections_over_t[t_idx]
                
                store_t_idx = int((t_idx+1) / time_scaling)
                
                Ca_cyt[hexa][store_t_idx] = new_Ca_cyt[hexa]
                ip3[hexa][store_t_idx] = new_ip3[hexa]
                Ca_stored[hexa][store_t_idx] = new_Ca_stored[hexa]
                ip3R_act[hexa][store_t_idx] = new_ip3R_act[hexa]
                
    
    return Ca_cyt, ip3, Ca_stored, ip3R_act,

def politi(variables, run_t, store_t, hex_array, params):
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
    
    old_Ca_cyt = allocate_var_dict(hex_array, 1, 0)
    old_ip3 = allocate_var_dict(hex_array, 1, 0)
    old_Ca_stored = allocate_var_dict(hex_array, 1, 0)
    old_ip3R_act = allocate_var_dict(hex_array, 1, 0)
    
    for hexa in hex_array: 
        old_Ca_cyt[hexa] = Ca_cyt[hexa][0]
        old_ip3[hexa] = ip3[hexa][0]
        old_Ca_stored[hexa] = Ca_stored[hexa][0]
        old_ip3R_act[hexa] = ip3R_act[hexa][0]
        
    new_Ca_cyt = allocate_var_dict(hex_array, 1, 0)
    new_ip3 = allocate_var_dict(hex_array, 1, 0)
    new_Ca_stored = allocate_var_dict(hex_array, 1, 0)
    new_ip3R_act = allocate_var_dict(hex_array, 1, 0)
    
    ##################################################################
    
    tau_p = 1 / (k_3K + k_5P)
    eta = k_3K * tau_p
    
    Ca_temp = allocate_var_dict(hex_array, 1, 0)
    V_PLC_scaled = allocate_var_dict(hex_array, 1, 0)
    
    dt = run_t[1] - run_t[0]
    store_dt = store_t[1] - store_t[0]
    time_scaling = int(store_dt / dt)
    
    for t_idx, t in enumerate(run_t[1:]):
    
        for hexa in hex_array:
        
            V_PLC_scaled[hexa] = V_PLC[hexa] * tau_p
        
            Ca_temp[hexa] = (k_1 * ((Hill(old_ip3R_act[hexa], K_a, old_Ca_cyt[hexa], 1 ) * Hill(1, K_p , old_ip3[hexa], 1))**3) + k_2) * (old_Ca_stored[hexa] - old_Ca_cyt[hexa]) - Hill(V_SERCA, K_SERCA, old_Ca_cyt[hexa], 2)
        
            new_ip3R_act[hexa] = old_ip3R_act[hexa] + dt * ( (1 / tau_r) * (1 - old_ip3R_act[hexa] * ((K_i + old_Ca_cyt[hexa]) / K_i) ) )
        
            new_Ca_stored[hexa] = old_Ca_stored[hexa] + dt * (((-1) * Ca_temp[hexa]) * (1 / beta) )
        
            new_Ca_cyt[hexa] = old_Ca_cyt[hexa] + dt * ( Ca_temp[hexa] + epsilon * (v_0 + phi * V_PLC_scaled[hexa] - Hill(V_pm, K_pm, old_Ca_cyt[hexa], 2)) )
        
            # only IP3 travels between cells
            # this calculates an ip3 neighborhood average with no-flux boundary conditions
            ip3_neighborhood_sum = 0
            ip3_neighbor_counter = 0
            for direction in range(6):
                neighbor = hex_neighbor(hexa, direction)
                try:
                    ip3_neighborhood_sum += old_ip3[neighbor]
                    ip3_neighbor_counter += 1
                except KeyError:
                    continue
            ip3_neighborhood_avg = ip3_neighborhood_sum/ip3_neighbor_counter
        
            new_ip3[hexa] = old_ip3[hexa] + dt * ( (1 / tau_p) * (Hill(V_PLC_scaled[hexa], K_PLC, old_Ca_cyt[hexa], 2) - ( old_ip3[hexa] * ( Hill(eta, K_3K, old_Ca_cyt[hexa], 2) + 1 - eta) ) ) + D_IP3 * (ip3_neighborhood_avg - old_ip3[hexa]) )
            
            old_Ca_cyt[hexa] = new_Ca_cyt[hexa]
            old_ip3[hexa] = new_ip3[hexa]
            old_Ca_stored[hexa] = new_Ca_stored[hexa]
            old_ip3R_act[hexa] = new_ip3R_act[hexa]
            
            if t in store_t:
                store_t_idx = int((t_idx+1) / time_scaling)
                
                Ca_cyt[hexa][store_t_idx] = new_Ca_cyt[hexa]
                ip3[hexa][store_t_idx] = new_ip3[hexa]
                Ca_stored[hexa][store_t_idx] = new_Ca_stored[hexa]
                ip3R_act[hexa][store_t_idx] = new_ip3R_act[hexa]
                
    
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
        # var_dict[hexa][0] = random.randint(value_range[0], value_range[1])
        var_dict[hexa][0] = np.random.uniform(low=value_range[0], high=value_range[1])
        direction_dict[hexa] = random.choice(directions)
        
    
    value_increment = (4 * (value_range[1] - value_range[0])) / (timepoint_N - 1) # for looping purposes
    for t in range(1, timepoint_N):
        for hexa in hex_array:
            if var_dict[hexa][t-1] <= value_range[0] or var_dict[hexa][t-1] >= value_range[1]:
                direction_dict[hexa] *= -1
            var_dict[hexa][t] = var_dict[hexa][t - 1] + (value_increment * direction_dict[hexa])
    
    return var_dict
    
    
    