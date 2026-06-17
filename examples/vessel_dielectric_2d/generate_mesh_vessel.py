#!/usr/bin/env python3

"""Create axisymmetric mesh for a vessel containing two disc electrodes
(with semi-circle sides) connected to rods, one on top and one on the bottom.
"""

import numpy as np
import gmsh
from skfem import Basis, BilinearForm, solve, condense
from skfem import ElementTriP0, ElementTriP1, ElementTriP3
from skfem.io.meshio import from_meshio, to_meshio
from skfem.helpers import dot, grad
import meshio
import argparse


def create_vessel_mesh(args):
    gmsh.initialize()
    model = gmsh.model
    geo = model.geo
    model.add("vessel")

    vessel_r = args.vessel_r
    vessel_h = args.vessel_h
    rod1_r = args.rod1_r
    rod1_h = args.rod1_h
    disc1_r = args.disc1_r
    disc1_h = args.disc1_h
    rod2_r = args.rod2_r
    rod2_h = args.rod2_h
    disc2_r = args.disc2_r
    disc2_h = args.disc2_h

    if disc1_r - 0.5 * disc1_h - rod1_r < 0:
        raise ValueError("disc1_r is smaller than supported")
    if disc2_r - 0.5 * disc2_h - rod2_r < 0:
        raise ValueError("disc2_r is smaller than supported")

    # --------------------------------------------------------------------
    # POINTS
    # --------------------------------------------------------------------

    # Axis points
    p_axis_disc1_top = geo.addPoint(0.0, rod1_h + disc1_h, 0.0, args.lc_disc)
    p_axis_disc2_bottom = geo.addPoint(0.0, vessel_h - (rod2_h + disc2_h),
                                       0.0, args.lc_disc)

    # Outer domain boundary along r = vessel_r
    p_outer_bottom = geo.addPoint(vessel_r, 0.0, 0.0, args.lc_vessel)
    p_outer_top = geo.addPoint(vessel_r, vessel_h, 0.0, args.lc_vessel)

    # Bottom rod (r = rod1_r, 0 <= z <= rod1_h)
    p_rod1_bottom = geo.addPoint(rod1_r, 0.0, 0.0, args.lc_rod)
    p_rod1_top = geo.addPoint(rod1_r, rod1_h, 0.0, args.lc_rod)

    # Bottom disc (semi-circle side), attached at z = rod1_h to rod1
    r_c1 = disc1_r - 0.5 * disc1_h
    z_c1 = rod1_h + 0.5 * disc1_h

    # Points on the semi-circle
    # bottom-right of the semi-circle (z = rod1_h)
    p_disc1_br = geo.addPoint(r_c1, rod1_h, 0.0, args.lc_disc)
    # rightmost point of semi-circle (r = disc1_r, z = z_c1)
    p_disc1_mid = geo.addPoint(disc1_r, z_c1, 0.0, args.lc_disc)
    # top-right of the semi-circle (z = rod1_h + disc1_h)
    p_disc1_tr = geo.addPoint(r_c1, rod1_h + disc1_h, 0.0, args.lc_disc)

    # top axis intersection for disc1
    p_disc1_axis_top = p_axis_disc1_top  # same as defined above

    # Top rod vertical extent: [vessel_h - rod2_h, vessel_h]
    p_rod2_bottom = geo.addPoint(rod2_r, vessel_h - rod2_h, 0.0, args.lc_rod)
    p_rod2_top = geo.addPoint(rod2_r, vessel_h, 0.0, args.lc_rod)

    # Disc2 vertical extent: [vessel_h - (rod2_h + disc2_h), vessel_h - rod2_h]
    z_disc2_bottom = vessel_h - (rod2_h + disc2_h)
    z_disc2_top = vessel_h - rod2_h
    r_c2 = disc2_r - 0.5 * disc2_h
    z_c2 = z_disc2_bottom + 0.5 * disc2_h

    # Points on the semi-circle of the top disc
    p_disc2_tr = geo.addPoint(r_c2, z_disc2_top, 0.0, args.lc_disc)
    p_disc2_mid = geo.addPoint(disc2_r, z_c2, 0.0, args.lc_disc)
    p_disc2_br = geo.addPoint(r_c2, z_disc2_bottom, 0.0, args.lc_disc)

    # bottom axis intersection for disc2
    p_disc2_axis_bottom = p_axis_disc2_bottom

    # --------------------------------------------------------------------
    # LINES AND ARCS
    # --------------------------------------------------------------------

    # Outer boundary (r = vessel_r)
    l_outer_bottom = geo.addLine(p_outer_bottom, p_outer_top)

    # Bottom boundary: from axis to rod1 to outer wall
    l_bottom_rod1_outer = geo.addLine(p_rod1_bottom, p_outer_bottom)

    # Top boundary: from axis to rod2 to outer wall
    l_top_rod2_outer = geo.addLine(p_rod2_top, p_outer_top)

    # Bottom rod (electrode) vertical side
    l_rod1_side = geo.addLine(p_rod1_bottom, p_rod1_top)

    # Connection from rod1 top to disc1 bottom-right point
    l_rod1_to_disc1 = geo.addLine(p_rod1_top, p_disc1_br)

    # Bottom disc semi-circle arcs
    p_disc1_center = geo.addPoint(r_c1, z_c1, 0.0, args.lc_disc)
    a_disc1_br_mid = geo.addCircleArc(p_disc1_br, p_disc1_center, p_disc1_mid)
    a_disc1_mid_tr = geo.addCircleArc(p_disc1_mid, p_disc1_center, p_disc1_tr)

    # Top of bottom disc: from tr to axis
    l_disc1_tr_axis = geo.addLine(p_disc1_tr, p_disc1_axis_top)

    # Top rod (electrode) vertical side
    l_rod2_side = geo.addLine(p_rod2_bottom, p_rod2_top)

    # Connection from disc2 top-right to rod2 bottom
    l_disc2_to_rod2 = geo.addLine(p_disc2_tr, p_rod2_bottom)

    # Top disc semi-circle arcs
    p_disc2_center = geo.addPoint(r_c2, z_c2, 0.0, args.lc_disc)
    a_disc2_tr_mid = geo.addCircleArc(p_disc2_tr, p_disc2_center, p_disc2_mid)
    a_disc2_mid_br = geo.addCircleArc(p_disc2_mid, p_disc2_center, p_disc2_br)

    # Bottom of top disc: from axis to br
    l_disc2_axis_br = geo.addLine(p_disc2_axis_bottom, p_disc2_br)

    # --------------------------------------------------------------------
    # Optional components
    # --------------------------------------------------------------------

    if args.use_ellipsoid:
        z = args.ellips_z
        h = args.ellips_h
        r = args.ellips_r
        ellips_cnt = geo.addPoint(0.0, z, 0.0, args.lc_ellips)
        ellips_bot = geo.addPoint(0.0, z-0.5*h, 0.0, args.lc_ellips)
        ellips_top = geo.addPoint(0.0, z+0.5*h, 0.0, args.lc_ellips)
        ellips_mid = geo.addPoint(r, z, 0.0, args.lc_ellips)
        ellips_bot_mid = geo.addEllipseArc(ellips_bot, ellips_cnt,
                                           ellips_top, ellips_mid)
        ellips_mid_top = geo.addEllipseArc(ellips_mid, ellips_cnt,
                                           ellips_top, ellips_top)
        ellips_axis = geo.addLine(ellips_bot, ellips_top)
        l_ellips = [ellips_mid_top, -ellips_axis, ellips_bot_mid]

        # On-axis line but around ellipsoid
        l_axis_disc1_disc2 = [
            geo.addLine(p_axis_disc2_bottom, ellips_top),
            ellips_mid_top,
            ellips_bot_mid,
            geo.addLine(ellips_bot, p_axis_disc1_top),
        ]

        cl_ellips = geo.addCurveLoop(l_ellips)
    else:
        # Axis line
        l_axis_disc1_disc2 = [
            geo.addLine(p_axis_disc2_bottom, p_axis_disc1_top)
        ]
        cl_ellips = None

    if args.use_ellipsoid:
        # Gas boundary loop, counter-clockwise
        cl_gas = geo.addCurveLoop([
            l_outer_bottom,
            -l_top_rod2_outer,
            -l_rod2_side,
            -l_disc2_to_rod2,
            a_disc2_tr_mid,
            a_disc2_mid_br,
            -l_disc2_axis_br,
            geo.addLine(p_axis_disc2_bottom, ellips_top),
            -ellips_mid_top,
            -ellips_bot_mid,
            geo.addLine(ellips_bot, p_axis_disc1_top),
            -l_disc1_tr_axis,
            -a_disc1_mid_tr,
            -a_disc1_br_mid,
            -l_rod1_to_disc1,
            -l_rod1_side,
            l_bottom_rod1_outer
        ])
        s_gas = geo.addPlaneSurface([cl_gas])
        s_ellips = geo.addPlaneSurface([cl_ellips])
    else:
        # Gas boundary loop, counter-clockwise
        cl_gas = geo.addCurveLoop([
            l_outer_bottom,
            -l_top_rod2_outer,
            -l_rod2_side,
            -l_disc2_to_rod2,
            a_disc2_tr_mid,
            a_disc2_mid_br,
            -l_disc2_axis_br,
            *l_axis_disc1_disc2,
            -l_disc1_tr_axis,
            -a_disc1_mid_tr,
            -a_disc1_br_mid,
            -l_rod1_to_disc1,
            -l_rod1_side,
            l_bottom_rod1_outer
        ])
        s_gas = geo.addPlaneSurface([cl_gas])

    # --------------------------------------------------------------------
    # PHYSICAL GROUPS
    # --------------------------------------------------------------------

    bottom_electrode_lines = [
        l_rod1_side,
        l_rod1_to_disc1,
        a_disc1_br_mid,
        a_disc1_mid_tr,
        l_disc1_tr_axis
    ]

    top_electrode_lines = [
        l_disc2_axis_br,
        a_disc2_mid_br,
        a_disc2_tr_mid,
        l_disc2_to_rod2,
        l_rod2_side
    ]

    # Domain boundaries
    domain_top_lines = [l_top_rod2_outer]
    domain_bottom_lines = [l_bottom_rod1_outer]
    domain_side_lines = [l_outer_bottom]

    geo.synchronize()

    model.addPhysicalGroup(1, bottom_electrode_lines, name="bottom_electrode")
    model.addPhysicalGroup(1, top_electrode_lines, name="top_electrode")
    model.addPhysicalGroup(1, domain_top_lines, name="domain_top")
    model.addPhysicalGroup(1, domain_bottom_lines, name="domain_bottom")
    model.addPhysicalGroup(1, domain_side_lines, name="domain_side")
    model.addPhysicalGroup(2, [s_gas], name="gas")

    if args.use_ellipsoid:
        model.addPhysicalGroup(2, [s_ellips], name="ellipsoid")

    # --------------------------------------------------------------------
    # MESH
    # --------------------------------------------------------------------
    model.mesh.generate(2)
    gmsh.write(args.gmshfile)
    gmsh.finalize()


def solve_field_fem(args):
    mio = meshio.read(args.gmshfile)
    mesh = from_meshio(mio)
    basis = Basis(mesh, ElementTriP3())
    grounded_bnd = basis.get_dofs("bottom_electrode") + \
        basis.get_dofs("domain_bottom") + \
        basis.get_dofs("domain_side")
    top_electrode_bnd = basis.get_dofs("top_electrode")

    # Set epsilon in domain
    eps = np.zeros(mesh.nelements)
    eps[mesh.subdomains["gas"]] = 1.0

    # Set material in domain (0: gas)
    material = np.zeros(mesh.nelements, dtype=int)

    if "ellipsoid" in mesh.subdomains:
        ellipsoid = mesh.subdomains["ellipsoid"]
        eps[ellipsoid] = args.ellips_eps
        material[ellipsoid] = 1

    # To use epsilon in Laplace operator
    basis0 = basis.with_element(ElementTriP0())
    eps_field = basis0.interpolate(eps)

    @BilinearForm
    def laplace(u, v, w):
        r = w.x[0]
        return w.eps * dot(grad(u), grad(v)) * r

    A = laplace.assemble(basis, eps=eps_field)

    u = basis.zeros()
    u[top_electrode_bnd] = args.V_top
    u[grounded_bnd] = 0.0

    u = solve(*condense(A, x=u, D=top_electrode_bnd+grounded_bnd))

    # Save output in P1 basis
    basis1 = basis.with_element(ElementTriP1())

    # Compute field only in the gas solution; it will be zero elsewhere
    gas_basis = Basis(mesh, ElementTriP3(), elements=mesh.subdomains["gas"])
    gas_basis1 = gas_basis.with_element(ElementTriP1())
    uh_gas = gas_basis.interpolate(u)
    Er_gas = gas_basis1.project(-uh_gas.grad[0])
    Ez_gas = gas_basis1.project(-uh_gas.grad[1])

    m_out = to_meshio(mesh, encode_cell_data=False)
    m_out.cell_data["material"] = [material]
    m_out.cell_data["eps"] = [eps]
    m_out.point_data["potential"] = basis1.project(basis.interpolate(u))
    m_out.point_data["E_r_gas"] = Er_gas
    m_out.point_data["E_z_gas"] = Ez_gas
    m_out.point_data["E_norm_gas"] = (Er_gas**2 + Ez_gas**2)**0.5

    meshio.write(args.output, m_out)
    print(f"Solution written to {args.output}")


if __name__ == '__main__':

    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Build axisymmetric sphere geometry')
    p.add_argument('--V_top', type=float, default=1.0,
                   help='Potential on top electrode [V]')
    p.add_argument('--rod1_r', type=float, default=5e-3,
                   help='Radius of bottom rod')
    p.add_argument('--rod1_h', type=float, default=20e-3,
                   help='Height of bottom rod')
    p.add_argument('--disc1_h', type=float, default=10e-3,
                   help='Height of bottom disc')
    p.add_argument('--disc1_r', type=float, default=40e-3,
                   help='Radius of bottom disc')
    p.add_argument('--rod2_r', type=float, default=5e-3,
                   help='Radius of top rod')
    p.add_argument('--rod2_h', type=float, default=20e-3,
                   help='Height of top rod')
    p.add_argument('--disc2_h', type=float, default=10e-3,
                   help='Height of top disc')
    p.add_argument('--disc2_r', type=float, default=40e-3,
                   help='Radius of top disc')
    p.add_argument('--vessel_r', type=float, default=100.0e-3,
                   help='Radius of vessel')
    p.add_argument('--vessel_h', type=float, default=150.0e-3,
                   help='Height of vessel')
    p.add_argument('--lc_disc', type=float, default=1e-3,
                   help='Mesh spacing on discs')
    p.add_argument('--lc_vessel', type=float, default=5e-3,
                   help='Mesh spacing on vessel')
    p.add_argument('--lc_rod', type=float, default=2e-3,
                   help='Mesh spacing on rod')
    p.add_argument('--lc_ellips', type=float, default=0.5e-3,
                   help='Mesh spacing on ellipsoid')
    p.add_argument('--use_ellipsoid', action='store_true',
                   help='Whether to include an dielectric ellipsoid')
    p.add_argument('--ellips_z', type=float, default=75.0e-3,
                   help='z-coordinate of ellipsoid center')
    p.add_argument('--ellips_h', type=float, default=20.0e-3,
                   help='Height of ellipsoid')
    p.add_argument('--ellips_r', type=float, default=5.0e-3,
                   help='Radius of ellipsoid')
    p.add_argument('--ellips_eps', type=float, default=2.0,
                   help='Relative permittivity of ellipsoid')
    p.add_argument('--gmshfile', type=str, default='vessel.msh',
                   help='File containing gmsh mesh')
    p.add_argument('--output', type=str, default='vessel.vtu',
                   help='File containing solution')
    args = p.parse_args()

    create_vessel_mesh(args)
    solve_field_fem(args)
