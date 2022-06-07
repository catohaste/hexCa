import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
from shapely.geometry import Polygon
from matplotlib.animation import FuncAnimation
from copy import deepcopy

from lib import *
from functions import *
from plot import set_axes_lims_from_hexes

def animate_var_flag(flag_var_dict, timepoint_N, flag_hexes, stripe_dim, pointy_layout, figsize_x, file_str):
    
    all_hexes = []
    for stripe in flag_hexes:
        all_hexes = all_hexes + flag_hexes[stripe]
    
    colormap_strings = {
        'black': 'Greys',
        'brown': 'copper',
        'red': 'Reds',
        'orange': 'Oranges',
        'yellow': 'YlOrBr',
        'green': 'Greens',
        'blue': 'Blues',
        'purple': 'Purples'
    }
    norms = {}
    cmaps = {}
    
    
    number_of_stripes = len(flag_hexes)
    
    grid_aspect_ratio = (stripe_dim[1] * number_of_stripes) / stripe_dim[0]
    fig = plt.figure(figsize=(figsize_x, figsize_x*grid_aspect_ratio))
    ax = fig.add_subplot(111)
    
    hex_patches = {}
    pointy_radius = pointy_layout.size[0]
    
    for stripe in flag_hexes:
        
        var_dict = flag_var_dict[stripe]
        hexes = flag_hexes[stripe]
    
        value_loc_tuples = create_val_loc_tuple_std_layout(var_dict, hexes, pointy_layout)
        
        color_str = colormap_strings[stripe]
        if stripe == 'black':
            # get colormap
            min_val = 0
            max_val = 1
            var_norm = mpl.colors.Normalize(vmin=min_val, vmax=max_val)
            var_cmap = plt.get_cmap(color_str)
        elif stripe == 'brown':
            min_val = min([min(val_array) for (val_array, center) in value_loc_tuples])
            max_val = max([max(val_array) for (val_array, center) in value_loc_tuples])
            width = max_val - min_val
            # var_norm = mpl.colors.Normalize(vmin=min_val - 2*width, vmax=max_val ) # YlOrBr
            var_norm = mpl.colors.Normalize(vmin=min_val, vmax=max_val + width) # copper
            var_cmap = plt.get_cmap(color_str)
        elif stripe == 'yellow':
            min_val = min([min(val_array) for (val_array, center) in value_loc_tuples])
            max_val = max([max(val_array) for (val_array, center) in value_loc_tuples])
            width = max_val - min_val
            var_norm = mpl.colors.Normalize(vmin=min_val, vmax=max_val + width)
            var_cmap = plt.get_cmap(color_str)
        else:
            # get colormap
            min_val = min([min(val_array) for (val_array, center) in value_loc_tuples])
            max_val = max([max(val_array) for (val_array, center) in value_loc_tuples])
            var_norm = mpl.colors.Normalize(vmin=min_val, vmax=max_val)
            var_cmap = plt.get_cmap(color_str)
            
    
        norms[stripe] = var_norm
        cmaps[stripe] = var_cmap
        
        value_loc_tuples = create_val_loc_tuple_std_layout(var_dict, hexes, pointy_layout)    
        
        
        for hexa in hexes:
            val = var_dict[hexa][0]
            center = hex_to_pixel(pointy_layout, hexa)
            hex_patches[hexa] = RegularPolygon((center.x, center.y), facecolor=var_cmap(var_norm(val)), numVertices=6, radius=pointy_radius, edgecolor='k')
            ax.add_patch(hex_patches[hexa])
    
    set_axes_lims_from_hexes(ax, all_hexes, pointy_layout)
    
    ax.set_aspect('equal')
    fig.patch.set_visible(False)
    ax.axis('off')
    plt.tight_layout()
    
    video_length = 30 # seconds
    fps = 48
    interval_from_fps = 1000/fps
    frames_N = video_length * fps
    sample_rate = int(np.floor(timepoint_N / frames_N))
    if sample_rate == 0:
        sample_rate = 1
        frames_N = timepoint_N
    # print("frames:" + str(frames_N) + ", timepoints:" + str(timepoint_N) + ", sample_rate:" + str(sample_rate))
    
    def animate(i):
        for stripe in flag_hexes:
            var_dict = flag_var_dict[stripe]
            hexes = flag_hexes[stripe] 
            
            var_norm = norms[stripe]
            var_cmap = cmaps[stripe]
            
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

def plot_var_flag(flag_var_dict, timepoint_idx, flag_hexes, stripe_dim, pointy_layout, figsize_x, save_dir):
    
    all_hexes = []
    for stripe in flag_hexes:
        all_hexes = all_hexes + flag_hexes[stripe]
    
    colormap_strings = {
        'black': 'Greys',
        'brown': 'copper',
        'red': 'Reds',
        'orange': 'Oranges',
        'yellow': 'YlOrBr',
        'green': 'Greens',
        'blue': 'Blues',
        'purple': 'Purples'
    }
    
    number_of_stripes = len(flag_hexes)
    
    grid_aspect_ratio = (stripe_dim[1] * number_of_stripes) / stripe_dim[0]
    fig = plt.figure(figsize=(figsize_x, figsize_x*grid_aspect_ratio))
    ax = fig.add_subplot(111)
    
    for stripe in flag_hexes:
        
        var_dict = flag_var_dict[stripe]
        hexes = flag_hexes[stripe]
    
        value_loc_tuples = create_val_loc_tuple_std_layout(var_dict, hexes, pointy_layout)
        
        color_str = colormap_strings[stripe]
        if stripe == 'black':
            # get colormap
            min_val = 0
            max_val = 1
            var_norm = mpl.colors.Normalize(vmin=min_val, vmax=max_val)
            var_cmap = plt.get_cmap(color_str)
        elif stripe == 'brown':
            min_val = min([min(val_array) for (val_array, center) in value_loc_tuples])
            max_val = max([max(val_array) for (val_array, center) in value_loc_tuples])
            width = max_val - min_val
            # var_norm = mpl.colors.Normalize(vmin=min_val - 2*width, vmax=max_val ) # YlOrBr
            var_norm = mpl.colors.Normalize(vmin=min_val, vmax=max_val + width) # copper
            var_cmap = plt.get_cmap(color_str)
        elif stripe == 'yellow':
            min_val = min([min(val_array) for (val_array, center) in value_loc_tuples])
            max_val = max([max(val_array) for (val_array, center) in value_loc_tuples])
            width = max_val - min_val
            var_norm = mpl.colors.Normalize(vmin=min_val, vmax=max_val + width)
            var_cmap = plt.get_cmap(color_str)
        else:
            # get colormap
            min_val = min([min(val_array) for (val_array, center) in value_loc_tuples])
            max_val = max([max(val_array) for (val_array, center) in value_loc_tuples])
            var_norm = mpl.colors.Normalize(vmin=min_val, vmax=max_val)
            var_cmap = plt.get_cmap(color_str)

        hex_centers = [hex_to_pixel(pointy_layout, hexa) for hexa in hexes]
    
        pointy_radius = pointy_layout.size[0]
        hex_patches = [RegularPolygon((center.x, center.y), facecolor=var_cmap(var_norm(val_array[timepoint_idx])), numVertices=6, radius=pointy_radius, edgecolor='k') for (val_array, center) in value_loc_tuples]
        for patch in hex_patches:
            ax.add_patch(patch)
    
    set_axes_lims_from_hexes(ax, all_hexes, pointy_layout)

    ax.set_aspect('equal')
    fig.patch.set_visible(False)
    ax.axis('off')
    plt.tight_layout()
    
    if save_dir == 'show':
        plt.show()
    else:
        plt.savefig(save_dir + '.png')
    
def plot_hexes_flag(hexes, hex_grid_dim, pointy_layout, figsize_x, save_dir):
    
    pointy_radius = pointy_layout.size[0]

    hex_centers = [hex_to_pixel(pointy_layout, hex) for hex in hexes]

    grid_aspect_ratio = hex_grid_dim[1] / hex_grid_dim[0]
    fig = plt.figure(figsize=(figsize_x, figsize_x*grid_aspect_ratio))
    ax = fig.add_subplot(111)
    
    # hex_patches = [RegularPolygon((center.x, center.y), facecolor='grey', numVertices=6, radius=pointy_radius, edgecolor='k', orientation=np.pi/6) for center in hex_centers]
    hex_patches = [RegularPolygon((center.x, center.y), facecolor='grey', numVertices=6, radius=pointy_radius, edgecolor='k', orientation=0) for center in hex_centers]
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
        