import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
from shapely.geometry import Polygon
from matplotlib.animation import FuncAnimation

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
    
def animate_var_by_color(var_dict, timepoint_N, hexes, hex_grid_dim, pointy_layout, figsize_x, color_str, file_str):
    
    value_loc_tuples = create_val_loc_tuple_std_layout(var_dict, hexes, pointy_layout)
    
    # get colormap
    min_val = min([min(val_array) for (val_array, center) in value_loc_tuples])
    max_val = max([max(val_array) for (val_array, center) in value_loc_tuples])
    var_norm = mpl.colors.Normalize(vmin=min_val, vmax=max_val)
    var_cmap = plt.get_cmap(color_str)
    
    grid_aspect_ratio = hex_grid_dim[1] / hex_grid_dim[0]
    fig = plt.figure(figsize=(figsize_x, figsize_x*grid_aspect_ratio))
    ax = fig.add_subplot(111)
    
    pointy_radius = pointy_layout.size[0]
    hex_patches = {}
    for hexa in hexes:
        val = var_dict[hexa][0]
        center = hex_to_pixel(pointy_layout, hexa)
        hex_patches[hexa] = RegularPolygon((center.x, center.y), facecolor=var_cmap(var_norm(val)), numVertices=6, radius=pointy_radius, edgecolor='k')
        ax.add_patch(hex_patches[hexa])
    
    set_axes_lims_from_hexes(ax, hexes, pointy_layout)
    
    ax.set_aspect('equal')
    fig.patch.set_visible(False)
    ax.axis('off')
    plt.tight_layout()
    
    video_length = 10 # seconds
    fps = 24
    interval_from_fps = 1000/fps
    frames_N = video_length * fps
    sample_rate = int(np.floor(timepoint_N / frames_N))
    # print("frames:" + str(frames_N) + ", timepoints:" + str(timepoint_N) + ", sample_rate:" + str(sample_rate))
    def animate(i):
        for hexa in hexes:
            val = var_dict[hexa][i*sample_rate]
            hex_patches[hexa].set_facecolor(var_cmap(var_norm(val)))
        return

    if file_str == 'show':
        print("animation not currently working")
        # OPTION 1
        from matplotlib import rc
        from IPython.display import HTML
        anim_jupyter = FuncAnimation(fig, animate, frames=frames_N, interval=interval_from_fps, blit=False)
        rc('animation', html='html5')
        HTML(anim_jupyter.to_html5_video())
        
        # OPTION 2
        # anim_jupyter = FuncAnimation(fig, animate, frames=frames_N, interval=interval_from_fps, blit=False)
        # plt.show()
        
        # OPTION 3
        # from IPython.display import HTML
        # anim_mp4 = FuncAnimation(fig, animate, frames=frames_N, blit=False)
        # anim_mp4.save('hex_anim.mp4', writer='ffmpeg', fps=fps)
        # HTML("""
        #     <video alt="test" controls>
        #         <source src="hex_anim.mp4" type="video/mp4">
        #     </video>
        # """)
    else:
        anim_mp4 = FuncAnimation(fig, animate, frames=frames_N, blit=False)
        anim_mp4.save(file_str + '.mp4', writer='ffmpeg', fps=fps)

def plot_var_by_color(var_dict, timepoint_idx, hexes, hex_grid_dim, pointy_layout, figsize_x, save_dir):
    
    value_loc_tuples = create_val_loc_tuple_std_layout(var_dict, hexes, pointy_layout)
    
    # get colormap
    min_val = min([min(val_array) for (val_array, center) in value_loc_tuples])
    max_val = max([max(val_array) for (val_array, center) in value_loc_tuples])
    var_norm = mpl.colors.Normalize(vmin=min_val, vmax=max_val)
    var_cmap = plt.get_cmap('Purples')

    hex_centers = [hex_to_pixel(pointy_layout, hexa) for hexa in hexes]
    
    grid_aspect_ratio = hex_grid_dim[1] / hex_grid_dim[0]
    fig = plt.figure(figsize=(figsize_x, figsize_x*grid_aspect_ratio))
    ax = fig.add_subplot(111)
    
    pointy_radius = pointy_layout.size[0]
    hex_patches = [RegularPolygon((center.x, center.y), facecolor=var_cmap(var_norm(val_array[timepoint_idx])), numVertices=6, radius=pointy_radius, edgecolor='k') for (val_array, center) in value_loc_tuples]
    for patch in hex_patches:
        ax.add_patch(patch)
    
    set_axes_lims_from_hexes(ax, hexes, pointy_layout)

    ax.set_aspect('equal')
    fig.patch.set_visible(False)
    ax.axis('off')
    plt.tight_layout()
    
    if save_dir == 'show':
        plt.show()
    else:
        plt.savefig(save_dir + 'hexes_random.png')
    
def plot_hexes(hexes, hex_grid_dim, pointy_layout, figsize_x, save_dir):
    
    pointy_radius = pointy_layout.size[0]

    hex_centers = [hex_to_pixel(pointy_layout, hex) for hex in hexes]

    grid_aspect_ratio = hex_grid_dim[1] / hex_grid_dim[0]
    fig = plt.figure(figsize=(figsize_x, figsize_x*grid_aspect_ratio))
    ax = fig.add_subplot(111)

    hex_patches = [RegularPolygon((center.x, center.y), facecolor='grey', numVertices=6, radius=pointy_radius, edgecolor='k') for center in hex_centers]
    for patch in hex_patches:
        ax.add_patch(patch)
        
    set_axes_lims_from_hexes(ax, hexes, pointy_layout)
    
    ax.set_aspect('equal')
    fig.patch.set_visible(False)
    ax.axis('off')
    plt.tight_layout()
    
    if save_dir == 'show':
        plt.show()
    else:
        plt.savefig(save_dir + 'hexes_const.png')

