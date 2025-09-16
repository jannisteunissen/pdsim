#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
import argparse

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Plot pdsim results')
parser.add_argument('vtu_file', type=str, help='Input .vtu file')
parser.add_argument('-var', type=str, default='K_star',
                    help='Variable to plot')
parser.add_argument('-mesh', action='store_true',
                    help='Show the mesh')
parser.add_argument('-window_size', type=int, nargs=2, default=[1024, 768],
                    help='Size of the plot (in pixels)')
parser.add_argument('-savefig', type=str,
                    help='Save plot to this file')
parser.add_argument('-cmap', type=str, default='jet',
                    help='Colormap to use (e.g. viridis, jet, inferno, ...)')
parser.add_argument('-cbar_width', type=float, default=0.5,
                    help='Colorbar width')
parser.add_argument('-cbar_position', type=float, nargs=2,
                    default=[0.25, 0.85],
                    help='Colorbar position')
parser.add_argument('-font_size', type=int, default=20,
                    help='Font size')
parser.add_argument('-font_title_size', type=int, default=24,
                    help='Title font size')
parser.add_argument('-font_label_size', type=int, default=20,
                    help='Label font size')
parser.add_argument('-font_color', type=str, default='white',
                    help='Font color')
args = parser.parse_args()

mesh = pv.read(args.vtu_file)
mesh.set_active_scalars(args.var, preference='point')

off_screen = (args.savefig is not None)

pv.set_plot_theme(pv.themes.DocumentTheme())
pv.global_theme.cmap = args.cmap
pv.global_theme.font.size = args.font_size
pv.global_theme.font.title_size = args.font_title_size
pv.global_theme.font.label_size = args.font_label_size
pv.global_theme.font.color = args.font_color

if args.cbar_width is not None:
    pv.global_theme.colorbar_horizontal.width = args.cbar_width

if args.cbar_position is not None:
    pv.global_theme.colorbar_horizontal.position_x = args.cbar_position[0]
    pv.global_theme.colorbar_horizontal.position_y = args.cbar_position[1]

pl = pv.Plotter(off_screen=off_screen,
                window_size=args.window_size)
pl.add_mesh(mesh, show_edges=args.mesh, show_scalar_bar=False)
pl.camera_position = 'xy'
pl.camera.zoom('tight')
pl.add_scalar_bar(title=args.var)

if args.savefig is not None:
    pl.save_graphic(args.savefig)
else:
    pl.show()
