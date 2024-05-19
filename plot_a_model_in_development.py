import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.patches import RegularPolygon
from shapely.geometry import Polygon
from matplotlib.animation import FuncAnimation
from copy import deepcopy

from lib import *
from functions import *
from plot import set_axes_lims_from_hexes

# def truncate_colormap(cmap, minval=0.0, maxval=1.0, trunc_min=0,trunc_min=1, n=100):
#     new_cmap = colors.LinearSegmentedColormap.from_list(
#         'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
#         cmap(np.linspace(minval, maxval, n)))
#     return new_cmap
    
def truncate_norm_for_cmap(data_min_val, data_max_val, trunc_min, trunc_max):
    """ trunc_min and trunc_max between 0 and 1 """
    
    data_width = data_max_val - data_min_val
    trunc_width = trunc_max - trunc_min
    
    scaling = data_width / trunc_width
    additive_const = data_min_val - trunc_min * scaling
    
    cmap_min = additive_const
    cmap_max = scaling + additive_const
    
    var_norm = mpl.colors.Normalize(vmin=cmap_min, vmax=cmap_max)

    return var_norm
    
def get_truncated_norm(stripe, value_loc_tuples):
    
    data_min_val = min([min(val_array) for (val_array, center) in value_loc_tuples])
    data_max_val = max([max(val_array) for (val_array, center) in value_loc_tuples])
    
    if stripe == 'black':
        min_val = 0
        max_val = 1
        var_norm = mpl.colors.Normalize(vmin=min_val, vmax=max_val)
    elif stripe == 'brown':
        var_norm = truncate_norm_for_cmap(data_min_val, data_max_val, 0, 0.5)
    elif stripe == 'red':
        var_norm = truncate_norm_for_cmap(data_min_val, data_max_val, 0.2, 1)
    elif stripe == 'orange':
        var_norm = truncate_norm_for_cmap(data_min_val, data_max_val, 0.3, 0.8)
    elif stripe == 'yellow':
        var_norm = truncate_norm_for_cmap(data_min_val, data_max_val, 0.2, 1)
    elif stripe == 'green':
        var_norm = truncate_norm_for_cmap(data_min_val, data_max_val, 0.2, 1)
    elif stripe == 'blue':
        var_norm = mpl.colors.Normalize(vmin=0, vmax=1)
    elif stripe == 'purple':
        var_norm = truncate_norm_for_cmap(data_min_val, data_max_val, 0.3, 0.9)
    else:
        var_norm = mpl.colors.Normalize(vmin=data_min_val, vmax=data_max_val)
        
    return var_norm

def animate_var_flag(run_selection, flag_var_dict, connect_var, timepoint_N, flag_hexes, animation_type, stripe_dim, pointy_layout, figsize_x, file_str):
    
    all_hexes = []
    for stripe in run_selection:
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
    
    number_of_stripes = len(run_selection)
    
    grid_aspect_ratio = (stripe_dim[1] * number_of_stripes) / stripe_dim[0]
    fig = plt.figure(figsize=(figsize_x, figsize_x*grid_aspect_ratio), layout='tight')
    ax = fig.add_subplot(111)
    
    hex_patches = {}
    pointy_radius = pointy_layout.size[0]
    
    # sort out the colourmap with truncated min and maxs for cmap norm
    # for connections type animation, only single value is needed
    for stripe in run_selection:
    
        var_dict = flag_var_dict[stripe]
        hexes = flag_hexes[stripe]        

        color_str = colormap_strings[stripe]
        var_cmap = plt.get_cmap(color_str)
        cmaps[stripe] = var_cmap

        value_loc_tuples = create_val_loc_tuple_std_layout(var_dict, hexes, pointy_layout)

        var_norm = get_truncated_norm(stripe, value_loc_tuples)
        norms[stripe] = var_norm

        for hexa in hexes:
            val = var_dict[hexa][0]
            center = hex_to_pixel(pointy_layout, hexa)
            hex_patches[hexa] = RegularPolygon((center.x, center.y), facecolor=var_cmap(var_norm(val)), numVertices=6, radius=pointy_radius, edgecolor='k')
            ax.add_patch(hex_patches[hexa])
            
        if animation_type[stripe] == 'connect':
            
            current_connect_var = connect_var[stripe]
            initial_connections = current_connect_var['initial_connections']
            birth_connections = current_connect_var['birth_connections']
            death_connections = current_connect_var['death_connections']
            
            connection_timepoints = list(set(list(birth_connections) + list(death_connections) + [0]))
            connection_timepoints.sort()
            
            edge_line_dict = {}
    
            # add initial connections
            for edge in initial_connections.edges:
        
                point1 = hex_to_pixel(pointy_layout, edge[0])
                point2 = hex_to_pixel(pointy_layout, edge[1])
        
                x_coords = [point[0] for point in (point1,point2)]
                y_coords = [point[1] for point in (point1,point2)]
        
                l, = ax.plot(x_coords, y_coords, color=var_cmap(0.95))
        
                edge_line_dict[edge] = l
                lines = [edge_line_dict[edge] for edge in edge_line_dict]
    
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
        video_length = frames_N / fps
    print("frames:" + str(frames_N) + ", timepoints:" + str(timepoint_N) + ", sample_rate:" + str(sample_rate),  ", video_length: ", video_length)
    
    def animate(i):
        for stripe in run_selection:
            var_dict = flag_var_dict[stripe]
            hexes = flag_hexes[stripe] 
        
            var_norm = norms[stripe]
            var_cmap = cmaps[stripe]
        
            for hexa in hexes:
                val = var_dict[hexa][i*sample_rate]
                hex_patches[hexa].set_facecolor(var_cmap(var_norm(val)))
                
            if  animation_type[stripe] == 'connect':

                current_t = i * sample_rate
                
                current_connect_var = connect_var[stripe]
                initial_connections = current_connect_var['initial_connections']
                birth_connections = current_connect_var['birth_connections']
                death_connections = current_connect_var['death_connections']
        
                try:
                    add_connections = birth_connections[current_t]
                    for connection in add_connections:
                        point1 = hex_to_pixel(pointy_layout, connection[0])
                        point2 = hex_to_pixel(pointy_layout, connection[1])
    
                        x_coords = [point[0] for point in (point1,point2)]
                        y_coords = [point[1] for point in (point1,point2)]
    
                        l, = ax.plot(x_coords, y_coords, color=var_cmap(0.95))
    
                        edge_line_dict[connection] = l
                except KeyError:
                    pass
        
                # remove death connections
                try:
                    remove_connections = death_connections[current_t]
                    for connection in remove_connections:
                        remove_line = edge_line_dict[connection]
                        remove_line.remove()
        
                        del edge_line_dict[connection]
                except KeyError:
                    pass
        
                lines = [edge_line_dict[edge] for edge in edge_line_dict]
                
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

def plot_var_flag(run_selection, flag_var_dict, connect_var, timepoint_idx, flag_hexes, animation_type, stripe_dim, pointy_layout, figsize_x, save_dir):
    
    all_hexes = []
    for stripe in run_selection:
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
    
    number_of_stripes = len(run_selection)
    
    grid_aspect_ratio = (stripe_dim[1] * number_of_stripes) / stripe_dim[0]
    fig = plt.figure(figsize=(figsize_x, figsize_x*grid_aspect_ratio), layout='tight')
    ax = fig.add_subplot(111)
    
    for stripe in run_selection:
        
        var_dict = flag_var_dict[stripe]
        hexes = flag_hexes[stripe] 
        
        color_str = colormap_strings[stripe]
        var_cmap = plt.get_cmap(color_str)
        
        value_loc_tuples = create_val_loc_tuple_std_layout(var_dict, hexes, pointy_layout)
        
        var_norm = get_truncated_norm(stripe, value_loc_tuples)

        hex_centers = [hex_to_pixel(pointy_layout, hexa) for hexa in hexes]
    
        pointy_radius = pointy_layout.size[0]
        hex_patches = [RegularPolygon((center.x, center.y), facecolor=var_cmap(var_norm(val_array[timepoint_idx])), numVertices=6, radius=pointy_radius, edgecolor='k') for (val_array, center) in value_loc_tuples]
        for patch in hex_patches:
            ax.add_patch(patch)
            
        if animation_type[stripe] == 'connect':
            
            current_connect_var = connect_var[stripe]
            initial_connections = current_connect_var['initial_connections']
            birth_connections = current_connect_var['birth_connections']
            death_connections = current_connect_var['death_connections']
            
            connection_timepoints = list(set(list(birth_connections) + list(death_connections) + [0]))
            connection_timepoints.sort()
            
            local_timepoint_N = len(connection_timepoints)
            
            edge_line_dict = {}
    
            for edge in initial_connections.edges:
        
                point1 = hex_to_pixel(pointy_layout, edge[0])
                point2 = hex_to_pixel(pointy_layout, edge[1])
        
                x_coords = [point[0] for point in (point1,point2)]
                y_coords = [point[1] for point in (point1,point2)]
        
                l, = ax.plot(x_coords, y_coords, color=var_cmap(0.95))
        
                edge_line_dict[edge] = l
                lines = [edge_line_dict[edge] for edge in edge_line_dict]
            
            for current_t in range(timepoint_idx):
        
                try:
                    add_connections = birth_connections[current_t]
                    for connection in add_connections:
                        point1 = hex_to_pixel(pointy_layout, connection[0])
                        point2 = hex_to_pixel(pointy_layout, connection[1])
        
                        x_coords = [point[0] for point in (point1,point2)]
                        y_coords = [point[1] for point in (point1,point2)]
        
                        l, = ax.plot(x_coords, y_coords, color=var_cmap(0.95))
        
                        edge_line_dict[connection] = l
                except KeyError:
                    pass
            
                # remove death connections
                try:
                    remove_connections = death_connections[current_t]
                    for connection in remove_connections:
                        remove_line = edge_line_dict[connection]
                        remove_line.remove()
            
                        del edge_line_dict[connection]
                except KeyError:
                    pass
            
                lines = [edge_line_dict[edge] for edge in edge_line_dict]
    
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
        