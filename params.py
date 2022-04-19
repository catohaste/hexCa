''' negative feedback parameter values'''
# params = {
#
#     ## IP3 dynamics parameters
#     "K_3K" : 0.4, # microMolar
#     "k_3K" : 0.1, # seconds^(-1)
#     "k_5P" : 0, # seconds^(-1)
#     "K_PLC" : 0, # microMolar
#
#     ## Ca2+ transport and structural parameters
#     "beta": 0.185,
#     "V_SERCA" : 0.25, # microMolar seconds^(-1)
#     "K_SERCA" : 0.1, # microMolar
#     "V_pm" : 0.01, # microMolar seconds^(-1)
#     "K_pm" : 0.12, # microMolar
#     "v_0" : 0.0004, # microMolar seconds^(-1)
#     "phi" : 0.045, # seconds^(-1)
#     "epsilon": 0,
#     "c_tot" : 2, # microMolar ###### unused
#
#     ## IP3R parameters
#     "k_1" : 7.4, # seconds^(-1)
#     "k_2" : 0.00148, # seconds^(-1)
#     "K_a" : 0.2, # microMolar
#     "K_i" : 0.3, # microMolar
#     "K_p" : 0.13, # microMolar
#     "tau_r" : 6.6, # seconds
#
#     "V_PLC" : 0.00045 # microMolar seconds^(-1)
#
# }

''' positive feedback parameter values'''
params = {

    ## IP3 dynamics parameters
    "K_3K" : 0.4, # microMolar
    "k_3K" : 0, # seconds^(-1)
    "k_5P" : 0.66, # seconds^(-1)
    "K_PLC" : 0.20, # microMolar

    ## Ca2+ transport and structural parameters
    "beta": 0.185,
    "V_SERCA" : 0.9, # microMolar seconds^(-1)
    "K_SERCA" : 0.1, # microMolar
    "V_pm" : 0.01, # microMolar seconds^(-1)
    "K_pm" : 0.12, # microMolar
    "v_0" : 0.0004, # microMolar seconds^(-1)
    "phi" : 0.0047, # seconds^(-1)
    "epsilon": 0,
    "c_tot" : 2, # microMolar ###### unused

    ## IP3R parameters
    "k_1" : 1.11, # seconds^(-1)
    "k_2" : 0.0203, # seconds^(-1)
    "K_a" : 0.08, # microMolar
    "K_i" : 0.4, # microMolar
    "K_p" : 0.13, # microMolar
    "tau_r" : 12.5, # seconds
    
    # V_PLC is defined for each cell
    # ought to be defined in main
    "V_PLC" : {}, # microMolar seconds^(-1)
    # 0.779 is a lower bound on oscillatory behaviour.
    # It isn't the maximal lower bound.
    # It also only applies when epsilon = 0

    "D_IP3": 0.0 # seconds^(-1)
    # "D_IP3": 0.02 # seconds^(-1)
}