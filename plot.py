import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
from shapely.geometry import Polygon

from lib import *
from functions import *

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
    
def animate_var_by_color(var_dict, timepoint_idx, hexes, hex_grid_dim, pointy_layout, figsize_x, save_dir):
    
    value_loc_tuples = create_val_loc_tuple_std_layout(var_dict, hexes, pointy_layout)
    
    # get colormap
    min_val = min([min(val_array) for (val_array, center) in value_loc_tuples])
    max_val = max([max(val_array) for (val_array, center) in value_loc_tuples])
    var_norm = mpl.colors.Normalize(vmin=min_val, vmax=max_val)
    var_cmap = plt.get_cmap('Oranges')

    hex_centers = [hex_to_pixel(pointy_layout, hexa) for hexa in hexes]
    
    grid_aspect_ratio = hex_grid_dim[1] / hex_grid_dim[0]
    fig = plt.figure(figsize=(figsize_x, figsize_x*grid_aspect_ratio))
    ax = fig.add_subplot(111)
    
    pointy_radius = pointy_layout.size[0]
    hex_patches = [RegularPolygon((center.x, center.y), facecolor=var_cmap(var_norm(val_array[timepoint_idx])), numVertices=6, radius=pointy_radius, alpha=0.2, edgecolor='k') for (val_array, center) in value_loc_tuples]
    for patch in hex_patches:
        ax.add_patch(patch)
    
    set_axes_lims_from_hexes(ax, hexes, pointy_layout)

    ax.set_aspect('equal')
    fig.patch.set_visible(False)
    ax.axis('off')
    plt.tight_layout()
    plt.savefig(save_dir + 'hexes_random.png')

def plot_var_by_color(var_dict, timepoint_idx, hexes, hex_grid_dim, pointy_layout, figsize_x, save_dir):
    
    value_loc_tuples = create_val_loc_tuple_std_layout(var_dict, hexes, pointy_layout)
    
    # get colormap
    min_val = min([min(val_array) for (val_array, center) in value_loc_tuples])
    max_val = max([max(val_array) for (val_array, center) in value_loc_tuples])
    var_norm = mpl.colors.Normalize(vmin=min_val, vmax=max_val)
    var_cmap = plt.get_cmap('Oranges')

    hex_centers = [hex_to_pixel(pointy_layout, hexa) for hexa in hexes]
    
    grid_aspect_ratio = hex_grid_dim[1] / hex_grid_dim[0]
    fig = plt.figure(figsize=(figsize_x, figsize_x*grid_aspect_ratio))
    ax = fig.add_subplot(111)
    
    pointy_radius = pointy_layout.size[0]
    hex_patches = [RegularPolygon((center.x, center.y), facecolor=var_cmap(var_norm(val_array[timepoint_idx])), numVertices=6, radius=pointy_radius, alpha=0.2, edgecolor='k') for (val_array, center) in value_loc_tuples]
    for patch in hex_patches:
        ax.add_patch(patch)
    
    set_axes_lims_from_hexes(ax, hexes, pointy_layout)

    ax.set_aspect('equal')
    fig.patch.set_visible(False)
    ax.axis('off')
    plt.tight_layout()
    plt.savefig(save_dir + 'hexes_random.png')
    
def plot_hexes(hexes, hex_grid_dim, pointy_layout, figsize_x, save_dir):
    
    pointy_radius = pointy_layout.size[0]

    hex_centers = [hex_to_pixel(pointy_layout, hex) for hex in hexes]

    grid_aspect_ratio = hex_grid_dim[1] / hex_grid_dim[0]
    fig = plt.figure(figsize=(figsize_x, figsize_x*grid_aspect_ratio))
    ax = fig.add_subplot(111)

    hex_patches = [RegularPolygon((center.x, center.y), facecolor='C0', numVertices=6, radius=pointy_radius, alpha=0.2, edgecolor='k') for center in hex_centers]
    for patch in hex_patches:
        ax.add_patch(patch)
        
    set_axes_lims_from_hexes(ax, hexes, pointy_layout)

    ax.set_aspect('equal')
    fig.patch.set_visible(False)
    ax.axis('off')
    plt.tight_layout()
    plt.savefig(save_dir + 'hexes_const.png')

