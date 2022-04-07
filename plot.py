import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
from shapely.geometry import Polygon

from lib import *
    
def plot_hexes(hexes, hex_grid_N, pointy_radius, fig_x):
    
    pointy = Layout(layout_pointy, Point(pointy_radius, pointy_radius), Point(0, 0))

    hex_centers = [hex_to_pixel(pointy, hex) for hex in hexes]
    hex_x_list = [center.x for center in hex_centers]
    hex_y_list = [center.y for center in hex_centers]

    grid_aspect_ratio = hex_grid_N[1] / hex_grid_N[0]
    fig = plt.figure(figsize=(fig_x, fig_x*grid_aspect_ratio))
    ax = fig.add_subplot(111)

    hex_patches = [RegularPolygon((center.x, center.y), facecolor='C0', numVertices=6, radius=pointy_radius, alpha=0.2, edgecolor='k') for center in hex_centers]
    for patch in hex_patches:
        ax.add_patch(patch)
    
    min_x = np.min(hex_x_list)
    max_x = np.max(hex_x_list)
    min_y = np.min(hex_y_list)
    max_y = np.max(hex_y_list)

    ax.set_xlim([min_x - pointy_radius, max_x + pointy_radius])
    ax.set_ylim([min_y - pointy_radius, max_y + pointy_radius])

    ax.set_aspect('equal')
    fig.patch.set_visible(False)
    ax.axis('off')
    plt.tight_layout()
    plt.show()

