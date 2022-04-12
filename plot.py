import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
from shapely.geometry import Polygon

from lib import *

def set_axes_lims_from_hexes(ax, hexes, pointy_layout):
    
    pointy_radius =  pointy_layout.size[0] # assumes unsquished haxagons
    
    hex_centers = [hex_to_pixel(pointy_layout, hexa) for hexa in hexes]
    hex_x_list = [center.x for center in hex_centers]
    hex_y_list = [center.y for center in hex_centers]

    min_x = np.min(hex_x_list)
    max_x = np.max(hex_x_list)
    min_y = np.min(hex_y_list)
    max_y = np.max(hex_y_list)
    ax.set_xlim([min_x - pointy_radius, max_x + pointy_radius])
    ax.set_ylim([min_y - pointy_radius, max_y + pointy_radius])

def plot_var_by_color(var_dict, hexes, hex_grid_dim, pointy_radius, figsize_x):
    
    pointy = Layout(layout_pointy, Point(pointy_radius, pointy_radius), Point(0, 0))
    
    value_loc_tuples = []
    for hexa in hexes:
        value_loc_tuples.append((var_dict[hexa], hex_to_pixel(pointy, hexa)))
    
    # get colormap
    min_val = min([val for (val, center) in value_loc_tuples])
    max_val = max([val for (val, center) in value_loc_tuples])
    var_norm = mpl.colors.Normalize(vmin=min_val, vmax=max_val)
    var_cmap = plt.get_cmap('Oranges')

    hex_centers = [hex_to_pixel(pointy, hexa) for hexa in hexes]

    grid_aspect_ratio = hex_grid_dim[1] / hex_grid_dim[0]
    fig = plt.figure(figsize=(figsize_x, figsize_x*grid_aspect_ratio))
    ax = fig.add_subplot(111)
    
    hex_patches = [RegularPolygon((center.x, center.y), facecolor=var_cmap(var_norm(val)), numVertices=6, radius=pointy_radius, alpha=0.2, edgecolor='k') for (val, center) in value_loc_tuples]
    for patch in hex_patches:
        ax.add_patch(patch)
    
    set_axes_lims_from_hexes(ax, hexes, pointy)

    ax.set_aspect('equal')
    fig.patch.set_visible(False)
    ax.axis('off')
    plt.tight_layout()
    plt.show()
    
def plot_hexes(hexes, hex_grid_dim, pointy_radius, figsize_x):
    
    pointy = Layout(layout_pointy, Point(pointy_radius, pointy_radius), Point(0, 0))

    hex_centers = [hex_to_pixel(pointy, hex) for hex in hexes]

    grid_aspect_ratio = hex_grid_dim[1] / hex_grid_dim[0]
    fig = plt.figure(figsize=(figsize_x, figsize_x*grid_aspect_ratio))
    ax = fig.add_subplot(111)

    hex_patches = [RegularPolygon((center.x, center.y), facecolor='C0', numVertices=6, radius=pointy_radius, alpha=0.2, edgecolor='k') for center in hex_centers]
    for patch in hex_patches:
        ax.add_patch(patch)
        
    set_axes_lims_from_hexes(ax, hexes, pointy)

    ax.set_aspect('equal')
    fig.patch.set_visible(False)
    ax.axis('off')
    plt.tight_layout()
    plt.show()

