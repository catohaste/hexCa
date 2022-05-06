import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
from shapely.geometry import Polygon
from matplotlib.animation import FuncAnimation
from copy import deepcopy

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

def plot_var_over_time_fixed_x_avg_y(var_dict, hexes, pointy_layout, figsize_x, color_str, var_str, file_str):
    
    value_loc_tuples = create_val_loc_tuple_std_layout(var_dict, hexes, pointy_layout)
    
    value_loc_dict_averaged_over_y = create_val_loc_dict_average_over_y(value_loc_tuples)
    
    x_idx_max = len(list(value_loc_dict_averaged_over_y.keys())) - 1
    x_locs = list(value_loc_dict_averaged_over_y.keys())
    
    x_locs.sort()
    
    selected_x_fractions = [0.1,0.3,0.5,0.7,0.9]    
    selected_x_idxs = [int(np.round(frac * x_idx_max)) for frac in selected_x_fractions]
    
    fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(figsize_x, 3))
    
    axs[0].set_ylabel(var_str + '\n(Mean over y-axis)', fontsize=14)
    
    for idx, ax in enumerate(axs):
        
        ax.set_title('x = ' + str(selected_x_fractions[idx]), fontsize=16)
        ax.set_xlabel('Time', fontsize=14)
        var_cmap = plt.get_cmap(color_str)
        ax.plot(value_loc_dict_averaged_over_y[x_locs[selected_x_idxs[idx]]], color=var_cmap(0.8))
    
    fig.tight_layout()
    
    if file_str == 'show':
        plt.show()
    else:
        plt.savefig(file_str + '_fixed_x_avg_over_y.png')
    
def animate_var_over_x_avg_y(var_dict, timepoint_N, hexes, hex_grid_dim, pointy_layout, figsize_x, color_str, var_str, file_str):
    
    value_loc_tuples = create_val_loc_tuple_std_layout(var_dict, hexes, pointy_layout)
    
    value_loc_dict_averaged_over_y = create_val_loc_dict_average_over_y(value_loc_tuples)
    
    fig = plt.figure(figsize=(figsize_x, 3))
    ax = fig.add_subplot(111)
    
    min_val = min([min(val_array) for (val_array, center) in value_loc_tuples])
    max_val = max([max(val_array) for (val_array, center) in value_loc_tuples])
    
    var_cmap = plt.get_cmap(color_str)
    
    ax.set_ylim([0, max_val*1.05])
    ax.set_xticklabels([])
    ax.set_xlabel('Cell position along x-axis', fontsize=14)
    ax.set_ylabel(var_str + '\n(Mean over y-axis)', fontsize=14)
        
    x = list(value_loc_dict_averaged_over_y.keys())
    x.sort()
    
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
        for artist in plt.gca().lines + plt.gca().collections:
            artist.remove()
        y = []
        for loc in x:
            y.append(value_loc_dict_averaged_over_y[loc][i*sample_rate])
        line, = ax.plot(x, y, color=var_cmap(0.8))
        return line, 
        
    if file_str == 'show':
        print("animation not currently working")
        # OPTION 1
        # from matplotlib import rc
        # from IPython.display import HTML
        # anim_jupyter = FuncAnimation(fig, animate, frames=frames_N, interval=interval_from_fps, blit=False)
        # rc('animation', html='html5')
        # HTML(anim_jupyter.to_html5_video())
        
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
        anim_mp4 = FuncAnimation(fig, animate, frames=frames_N, blit=True)
        # plt.show()
        anim_mp4.save(file_str + '_avg_over_y.mp4', writer='ffmpeg', fps=fps)
    
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
    
    hex_patches = [RegularPolygon((center.x, center.y), facecolor='grey', numVertices=6, radius=pointy_radius, edgecolor='k', orientation=np.pi/6) for center in hex_centers]
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
        
def plot_hexes_highlight_cells(cell_locs, hexes, hex_grid_dim, pointy_layout, figsize_x, save_dir):

    hex_array_minus_chosen = deepcopy(hexes)
    
    pointy_radius = pointy_layout.size[0]
    
    for cell_loc in cell_locs:
        chosen_h = OffsetCoord(cell_loc[0], cell_loc[1])
        chosen_hexa = roffset_to_cube(-1, chosen_h)
        if chosen_hexa in hexes:
            hex_array_minus_chosen.remove(chosen_hexa)
        else:
            print('Chosen cell '+ str(cell_loc)+' not in hex array')
            cell_locs.remove(cell_loc)
            continue

    hex_centers = [hex_to_pixel(pointy_layout, hexa) for hexa in hex_array_minus_chosen]

    grid_aspect_ratio = hex_grid_dim[1] / hex_grid_dim[0]
    fig = plt.figure(figsize=(figsize_x, figsize_x*grid_aspect_ratio))
    ax = fig.add_subplot(111)
    
    hex_patches = [RegularPolygon((center.x, center.y), facecolor='grey', numVertices=6, radius=pointy_radius, edgecolor='k', orientation=0) for center in hex_centers]
    # hex_patches = [RegularPolygon((center.x, center.y), facecolor='grey', numVertices=6, radius=pointy_radius, edgecolor='k', orientation=np.pi/6) for center in hex_centers] # flat layout
    for patch in hex_patches:
        ax.add_patch(patch)
    
    for cell_loc in cell_locs:
        chosen_h = OffsetCoord(cell_loc[0], cell_loc[1])
        chosen_hexa = roffset_to_cube(-1, chosen_h)
        chosen_hex_center = hex_to_pixel(pointy_layout, chosen_hexa)
        chosen_patch = RegularPolygon((chosen_hex_center.x, chosen_hex_center.y), facecolor='C3', numVertices=6, radius=pointy_radius, edgecolor='k', orientation=0)
        ax.add_patch(chosen_patch)
    
    set_axes_lims_from_hexes(ax, hexes, pointy_layout)
    
    ax.set_aspect('equal')
    fig.patch.set_visible(False)
    ax.axis('off')
    plt.tight_layout()
    
    if save_dir == 'show':
        plt.show()
    else:
        plt.savefig(save_dir + 'chosen_cell.png')

def plot_all_vars_over_time_single_cells(cell_locs, variables, var_strings, col_strings, file_str):
    
    Ca_cyt, ip3, Ca_stored, ip3R_act, = variables
    
    chosen_cell_N = len(cell_locs)
    var_N = len(variables)
    
    col_N = chosen_cell_N
    row_N = var_N
    ax_N = col_N * row_N
    
    fig, axs = plt.subplots(nrows=row_N, ncols=col_N, figsize=(4*col_N, 2*row_N))
    
    for var_idx in range(var_N):
    
        for cell_idx in range(chosen_cell_N):
        
            cell_loc = cell_locs[cell_idx]
            chosen_h = OffsetCoord(cell_loc[0], cell_loc[1])
            chosen_hexa = roffset_to_cube(-1, chosen_h)
            
            axs[0][cell_idx].set_title('cell at ' + str(cell_loc), fontsize=16)
            axs[var_N - 1][cell_idx].set_xlabel('Time', fontsize=14)
            
            ax = axs[var_idx][cell_idx]
            
            var_cmap = plt.get_cmap(col_strings[var_idx])
            ax.plot(variables[var_idx][chosen_hexa], color = var_cmap(0.8) )
            if cell_idx == 0:
                ax.set_ylabel(var_strings[var_idx])
    
    fig.tight_layout()
    
    if file_str == 'show':
        plt.show()
    else:
        plt.savefig(file_str + '_all_vars_chosen_cells.png')

def plot_var_running_time_avg_single_cells(cell_locs, running_avg_N, var_dict, var_str, color_str, file_str):
    
    chosen_cell_N = len(cell_locs)
    first_loc = cell_locs[0]
    first_h = OffsetCoord(first_loc[0], first_loc[1])
    first_hexa = roffset_to_cube(-1, first_h)
    timepoint_N = len(var_dict[first_hexa])
    
    if chosen_cell_N == 1: 
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4, 3))
        axs = [ax]
    else: 
        fig, axs = plt.subplots(nrows=1, ncols=chosen_cell_N, figsize=(4*chosen_cell_N, 3))
        
    for ax_idx, ax in enumerate(axs):
        
        cell_loc = cell_locs[ax_idx]
        chosen_h = OffsetCoord(cell_loc[0], cell_loc[1])
        chosen_hexa = roffset_to_cube(-1, chosen_h)
        
        running_avg = np.convolve(var_dict[chosen_hexa], np.ones(running_avg_N)/running_avg_N, mode='valid')
        time_running_avg = range(running_avg_N, timepoint_N + 1)
        
        ax.set_title('cell at ' + str(cell_locs[ax_idx]))
        ax.set_xlabel('Time', fontsize=14)
        var_cmap = plt.get_cmap(color_str)
        ax.plot(time_running_avg, running_avg, color=var_cmap(0.8))
        ax.set_xlim([0,timepoint_N])
        
    axs[0].set_ylabel(var_str + '\n(Running average)', fontsize=14)
    
    fig.tight_layout()
    
    if file_str == 'show':
        plt.show()
    else:
        plt.savefig(file_str + '_running_avg_single_cells.png')

def plot_links(links, hexes, hex_grid_dim, pointy_layout, figsize_x, color_str, file_str): 
        
    pointy_radius = pointy_layout.size[0]
    
    timepoint_N = len(links)

    hex_centers = [hex_to_pixel(pointy_layout, hex) for hex in hexes]

    grid_aspect_ratio = hex_grid_dim[1] / hex_grid_dim[0]
    fig = plt.figure(figsize=(figsize_x, figsize_x*grid_aspect_ratio))
    ax = fig.add_subplot(111)
    
    var_cmap = plt.get_cmap(color_str)
    
    hex_patches = [RegularPolygon((center.x, center.y), edgecolor=var_cmap(0.8), facecolor=var_cmap(0.4), alpha=0.5, numVertices=6, radius=pointy_radius) for center in hex_centers]
    # hex_patches = [RegularPolygon((center.x, center.y), facecolor='grey', numVertices=6, radius=pointy_radius, edgecolor='k', orientation=np.pi/6) for center in hex_centers]
    for patch in hex_patches:
        ax.add_patch(patch)
        
    set_axes_lims_from_hexes(ax, hexes, pointy_layout)
    
    video_length = 30 # seconds
    fps = 48
    interval_from_fps = 1000/fps
    frames_N = video_length * fps
    sample_rate = int(np.floor(timepoint_N / frames_N))
    if sample_rate == 0:
        sample_rate = 1
        frames_N = timepoint_N
    # print("frames:" + str(frames_N) + ", timepoints:" + str(timepoint_N) + ", sample_rate:" + str(sample_rate))
    plot_memory = 5
    current_time_line_tuples = []
    def animate(i):
        t_links = links[i*sample_rate]
        
        for time_line_tuple in current_time_line_tuples:
            time = time_line_tuple[0]
            line = time_line_tuple[1]
            if time + plot_memory < i:
                line.remove()
                current_time_line_tuples.remove(time_line_tuple)
                
        for link in t_links:
            point1 = hex_to_pixel(pointy_layout, link["hexes"][0])
            point2 = hex_to_pixel(pointy_layout, link["hexes"][1])
            x_coords = [point[0] for point in (point1,point2)]
            y_coords = [point[1] for point in (point1,point2)]
            l, = ax.plot(x_coords, y_coords, color='k')
            current_time_line_tuples.append((i,l))
            
        return
    
    ax.set_aspect('equal')
    fig.patch.set_visible(False)
    ax.axis('off')
    plt.tight_layout()
    
    if file_str == 'show':
        plt.show()
    else:
        anim_mp4 = FuncAnimation(fig, animate, frames=frames_N, blit=False)
        anim_mp4.save(file_str + '.mp4', writer='ffmpeg', fps=fps)
        
        