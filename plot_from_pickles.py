import pickle
import matplotlib.pyplot as plt
import random
import time

from lib import *
from functions import *
from plot import animate_var_by_color, animate_var_over_x_avg_y, plot_var_over_time_fixed_x_avg_y, plot_var_by_color, plot_colorbar


##################################################################################################
# results folder
results_dir = "results/"

save_dir = results_dir + 'V_PLC_gradient/0_787_to_1_100/'

##################################################################################################
# LOAD PICKLES

# var_strings = ['Ca_cyt', 'IP3']
# color_strings = ['Oranges', 'Blues']
var_strings = ['IP3']
color_strings = ['Blues']

pickle_dir = save_dir + "pickles/"

with open(pickle_dir + 'Ca_cyt' + '.pickle', 'rb') as handle:
    Ca_cyt_val_loc = pickle.load(handle)
with open(pickle_dir + 'ip3' + '.pickle', 'rb') as handle:
    ip3_val_loc = pickle.load(handle)
with open(pickle_dir + 'layout_dict.pickle', 'rb') as handle:
    layout_dict = pickle.load(handle)
pointy = create_layout_from_dict(layout_dict)

anim_param_str = ''

try:
    with open(pickle_dir + 'params.pickle', 'rb') as handle:
        params = pickle.load(handle)
    anim_param_str = r'$D_\mathrm{IP_3} $ = ' + str(params['D_IP3'])
    print('"D_IP3" loaded from pickle.')
except:
    params = {}
    params['D_IP3'] = 0.2
    anim_param_str = r'$D_\mathrm{IP_3} $ = ' + str(params['D_IP3'])
    print('"D_IP3" set manually.')

try:
    with open(pickle_dir + 'V_PLC.pickle', 'rb') as handle:
        V_PLC_val_loc_tuples = pickle.load(handle)
    V_PLC_values = [x[0] for x in V_PLC_val_loc_tuples]
    V_PLC_constant = all(i == V_PLC_values[0] for i in V_PLC_values)
    if V_PLC_constant:
        anim_param_str = r'$V_\mathrm{PLC} $ = ' + str(V_PLC_values[0])
    else:
        min_V_PLC = min(V_PLC_values)
        max_V_PLC = max(V_PLC_values)
        anim_param_str = r'$V_\mathrm{PLC} \in $ [' + str(min_V_PLC) + ', ' + str(max_V_PLC) + ']'
    print('"V_PLC" loaded from pickle.')
except:
    manual_V_PLC = 0.787
    anim_param_str = r'$V_\mathrm{PLC} $ = ' + str(manual_V_PLC)
    print('"V_PLC" set manually.')

hex_array_a, Ca_cyt = unpack_val_loc_tuple_std_layout(Ca_cyt_val_loc, pointy)
hex_array, ip3 = unpack_val_loc_tuple_std_layout(ip3_val_loc, pointy)

variables = [Ca_cyt, ip3]
variables = [ip3]

store_timepoint_N = 1200
store_timepoint_N = len(Ca_cyt[hex_array[0]])
store_dt = 0.5

hex_x_N,hex_y_N = 50, 30

##################################################################################################
# PLOT

start_time = time.time()

# plt.figure(figsize=(6, 6))
# plt.subplot(211)
# plt.plot(ip3[hex_array[0]])
# plt.subplot(212)
# plt.plot(ip3[hex_array[39]])
#
# plt.show()
# plot_time_indicies = [0, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210]
plot_time_indicies = [int(x/store_dt)for x in [0, 488, 494, 500, 506, 512, 518, 524, 530, 536, 542, 548]]
# plot_time_indicies = [0]

V_PLC_vary_min_max = [0.0598, 1.0702]
V_PLC_gradient_min_max = [0.0999, 1.492]

for var, var_str, color_str in zip(variables, var_strings, color_strings):
    # animate_var_by_color(var, store_timepoint_N, store_dt, hex_array, (hex_x_N,hex_y_N), pointy, 5, var_str, color_str, save_dir + var_str, show_time=True, show_colorbar=True, param_str=anim_param_str)
    # animate_var_over_x_avg_y(var, store_timepoint_N, hex_array, (hex_x_N,hex_y_N), pointy, 12, color_str, var_str, save_dir + var_str)
    # plot_var_over_time_fixed_x_avg_y(var, hex_array, pointy, 12, color_str, var_str, save_dir + var_str)
    plot_colorbar(var, 5, var_str, color_str, save_dir + var_str + '_colorbar', min_max=V_PLC_gradient_min_max)
    for plot_time_idx in plot_time_indicies:
        plot_var_by_color(var, plot_time_idx, store_dt, hex_array, (hex_x_N,hex_y_N), pointy, 12, color_str, save_dir + var_str + '_' + str(int(plot_time_idx*store_dt)), show_time=False, min_max=V_PLC_gradient_min_max)
    
end_time = time.time()

print("Time plotting", end_time - start_time)
    