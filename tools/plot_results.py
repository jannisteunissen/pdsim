#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
# Importing vtk before pyvista can avoid a bug in latex rendering, see
# https://github.com/pyvista/pyvista/issues/5636
import vtk
import pyvista as pv
import argparse

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Plot pdsim results')
parser.add_argument('vtu_file', type=str, help='Input .vtu file')
parser.add_argument('-var', type=str,
                    help='Variable to plot')
parser.add_argument('-mesh', action='store_true',
                    help='Show the mesh')
parser.add_argument('-axes_triad', action='store_true',
                    help='Show axes triad')
parser.add_argument('-mesh_color', type=str, default='gray',
                    help='Color for the mesh')
parser.add_argument('-mesh_style', type=str, default='surface',
                    help='Style of the mesh')
parser.add_argument('-line_width', type=float, default=1.0,
                    help='Line width')
parser.add_argument('-window_size', type=int, nargs=2, default=[1024, 768],
                    help='Size of the plot (in pixels)')
parser.add_argument('-savefig', type=str,
                    help='Save plot to this file')
parser.add_argument('-title', type=str,
                    help='Title of the plot')
parser.add_argument('-title_position', type=str, default='lower_edge',
                    help='Position of the title')
parser.add_argument('-zoom', type=float,
                    help='Zoom factor')
parser.add_argument('-focal_point', type=float, nargs=3,
                    help='Focal point of the camera')
parser.add_argument('-view_vector', type=float, nargs=3, default=[0., 0., -1.],
                    help='Point the camera in this direction')
parser.add_argument('-view_up', type=float, nargs=3,
                    help='View up direction of the camera')
parser.add_argument('-reflect', type=float, nargs=3,
                    help='Reflect the mesh along this direction')
parser.add_argument('-slice_normal', type=float, nargs=3,
                    help='Slice using this normal direction')
parser.add_argument('-slice_orthogonal', action='store_true',
                    help='Apply an orthogonal slice')
parser.add_argument('-multi_samples', type=int,
                    help='Number of multisamples for MSAA anti-aliasing')
parser.add_argument('-show_boundary', action='store_true',
                    help='Show boundary of the mesh')
parser.add_argument('-ruler', type=float, nargs=6,
                    help='Start and end position of ruler')
parser.add_argument('-distance_unit', type=str, default='mm',
                    choices=['mm', 'cm', 'm'],
                    help='Distance unit to use')
parser.add_argument('-cmap', type=str, default='plasma',
                    help=r'''Colormap to use (e.g. viridis, jet, inferno, ...),
see https://docs.pyvista.org/api/utilities/named_colormaps''')
parser.add_argument('-cbar_hide', action='store_true',
                    help='Hide colorbar')
parser.add_argument('-cbar_width', type=float, default=0.7,
                    help='Colorbar width')
parser.add_argument('-cbar_height', type=float,
                    help='Colorbar height')
parser.add_argument('-cbar_position', type=float, nargs=2,
                    default=[0.15, 0.85],
                    help='Colorbar position')
parser.add_argument('-cbar_title', type=str,
                    help='Colorbar title')
parser.add_argument('-cbar_range', type=float, nargs=2,
                    help='Colorbar range')
parser.add_argument('-cbar_n_labels', type=int, default=5,
                    help='Number of labels in colorbar')
parser.add_argument('-font_size', type=int, default=20,
                    help='Font size')
parser.add_argument('-font_title_size', type=int, default=28,
                    help='Title font size')
parser.add_argument('-font_label_size', type=int, default=20,
                    help='Label font size')
parser.add_argument('-font_color', type=str, default='white',
                    help='Font color')
parser.add_argument('-print_info', action='store_true',
                    help='Print information about the mesh')
args = parser.parse_args()

pv.set_plot_theme(pv.themes.DocumentTheme())
pv.global_theme.cmap = args.cmap
pv.global_theme.font.size = args.font_size
pv.global_theme.font.title_size = args.font_title_size
pv.global_theme.font.label_size = args.font_label_size
pv.global_theme.font.color = args.font_color
pv.global_theme.line_width = args.line_width

if args.multi_samples is not None:
    pv.global_theme.multi_samples = args.multi_samples

if args.cbar_width is not None:
    pv.global_theme.colorbar_horizontal.width = args.cbar_width

if args.cbar_height is not None:
    pv.global_theme.colorbar_horizontal.height = args.cbar_height

if args.cbar_position is not None:
    pv.global_theme.colorbar_horizontal.position_x = args.cbar_position[0]
    pv.global_theme.colorbar_horizontal.position_y = args.cbar_position[1]

mesh = pv.read(args.vtu_file)
mesh.set_active_scalars(args.var, preference='point')

if args.reflect is not None:
    mesh.reflect(args.reflect, inplace=True)

if args.print_info:
    print(mesh)

off_screen = (args.savefig is not None)

pl = pv.Plotter(off_screen=off_screen,
                window_size=args.window_size)

if args.slice_normal is not None:
    mesh = mesh.slice(normal=args.slice_normal, generate_triangles=True)
elif args.slice_orthogonal:
    mesh = mesh.slice_orthogonal(generate_triangles=True)

pl.add_mesh(mesh, show_edges=args.mesh, show_scalar_bar=False,
            edge_color=args.mesh_color, style=args.mesh_style,
            line_width=args.line_width)

if args.var is not None and not args.cbar_hide:
    pl.add_scalar_bar(title=args.cbar_title, n_labels=args.cbar_n_labels)
    if args.cbar_range is not None:
        pl.update_scalar_bar_range(args.cbar_range)

pl.enable_parallel_projection()
pl.view_vector(args.view_vector, args.view_up)
pl.enable_anti_aliasing('msaa', multi_samples=64)

if args.show_boundary:
    bnd = mesh.extract_feature_edges(boundary_edges=True,
                                     feature_edges=False,
                                     non_manifold_edges=False,
                                     manifold_edges=False,
                                     clear_data=True)
    pl.add_mesh(bnd, line_width=args.line_width)

if args.focal_point is not None:
    pl.camera.focal_point = args.focal_point

if args.ruler is not None:
    scale_factors = {'mm': 1e3, 'cm': 1e2, 'm': 1.0}
    pl.add_ruler(args.ruler[0:3], args.ruler[3:],
                 title=args.distance_unit,
                 scale=scale_factors[args.distance_unit])

if args.title is not None:
    pl.add_text(args.title, position=args.title_position)

if args.zoom is not None:
    pl.camera.zoom(args.zoom)
else:
    pl.camera.zoom('tight')

if args.axes_triad:
    pl.show_axes()

if args.savefig is not None:
    pl.save_graphic(args.savefig, raster=True)
else:
    pl.show()
