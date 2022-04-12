import random

from lib import *
from plot import plot_hexes, plot_var_by_color

# set up hexagonal grid with (q,r,s) coordinates
hex_x_N = 40
hex_y_N = 6
hex_array = []
for x in range(hex_x_N):
    for y in range(hex_y_N):
        hex_array.append(Hex(-x-y,y,x))
        
var_dict = {}

for hexa in hex_array:
    var_dict[hexa] = random.randint(4,95)
        
#plot_hexes(hex_array, hex_grid_dim, pointy_radius, figsize_x)
plot_var_by_color(var_dict, hex_array, (hex_x_N,hex_y_N), 10.0, 12)




