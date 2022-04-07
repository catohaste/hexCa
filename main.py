from lib import *

from plot import plot_hexes

# set up hexagonal grid with (q,r,s) coordinates
hex_x_N = 40
hex_y_N = 6
hexes = []
for x in range(hex_x_N):
    for y in range(hex_y_N):
        hexes.append(Hex(-x-y,y,x))
        
plot_hexes(hexes, (hex_x_N,hex_y_N), 10.0, 12)
