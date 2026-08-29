#!/usr/bin/env python3

"""Create axisymmetric mesh for a needle-to-plane geometry
"""

import numpy as np
import gmsh
from skfem import Basis, BilinearForm, solve, condense
from skfem import ElementTriP1, ElementTriP3
from skfem.io.meshio import from_meshio, to_meshio
from skfem.helpers import dot, grad
import meshio
import argparse


def create_needle_plate_mesh(args):
    if args.needle_tip_r >= args.needle_r:
        raise ValueError("Need needle_tip_r < needle_r")
    if args.d_gap <= 0:
        raise ValueError("Need d_gap > 0")
    if args.needle_cone_angle >= 90.0 or args.needle_cone_angle <= 0.0:
        raise ValueError("Need 0 < needle_cone_angle < 90")

    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 2)
    model = gmsh.model
    geo = model.geo
    model.add("needle_plate")

    cone_angle = args.needle_cone_angle * np.pi / 180.
    d_gap = args.d_gap
    needle_tip_r = args.needle_tip_r
    needle_r = args.needle_r
    domain_r = args.domain_r

    # Points
    p_tip_bottom = geo.addPoint(0.0, d_gap, 0.0, args.lc_tip)
    p_tip_center = geo.addPoint(0.0, d_gap + needle_tip_r, 0.0, args.lc_tip)

    # Get coordinate of endpoint of the tip
    r_tip = needle_tip_r * np.cos(-cone_angle)
    z_tip = d_gap + needle_tip_r - needle_tip_r * np.sin(cone_angle)
    p_tip_top = geo.addPoint(r_tip, z_tip, 0.0, args.lc_tip)

    # Cone height
    z_cone = z_tip + (needle_r - r_tip)/np.tan(cone_angle)
    p_cone_top = geo.addPoint(needle_r, z_cone, 0.0, args.lc_needle)

    z_top = z_cone + args.needle_h
    p_needle_top = geo.addPoint(needle_r, z_top, 0.0, args.lc_needle)

    p_top_right = geo.addPoint(domain_r, z_top, 0.0, args.lc_other)
    p_bottom_right = geo.addPoint(domain_r, 0.0, 0.0, args.lc_other)
    p_bottom_left = geo.addPoint(0.0, 0.0, 0.0, args.lc_tip)

    # Lines
    l_outer_bottom = geo.addLine(p_bottom_left, p_bottom_right)
    l_outer_right = geo.addLine(p_bottom_right, p_top_right)
    l_outer_top = geo.addLine(p_top_right, p_needle_top)
    l_needle_cyl = geo.addLine(p_needle_top, p_cone_top)
    l_needle_cone = geo.addLine(p_cone_top, p_tip_top)
    l_needle_tip = geo.addCircleArc(p_tip_top, p_tip_center, p_tip_bottom)
    l_axis = geo.addLine(p_tip_bottom, p_bottom_left)

    # Gas boundary loop, counter-clockwise
    cl_gas = geo.addCurveLoop([
            l_outer_bottom,
            l_outer_right,
            l_outer_top,
            l_needle_cyl,
            l_needle_cone,
            l_needle_tip,
            l_axis,
        ])
    s_gas = geo.addPlaneSurface([cl_gas])

    needle_lines = [
        l_needle_cyl,
        l_needle_cone,
        l_needle_tip,
    ]

    geo.synchronize()

    model.addPhysicalGroup(1, [l_outer_bottom], name="domain_bottom")
    model.addPhysicalGroup(1, needle_lines, name="needle")
    model.addPhysicalGroup(1, [l_outer_top], name="domain_top")
    model.addPhysicalGroup(1, [l_outer_right], name="domain_side")
    model.addPhysicalGroup(2, [s_gas], name="gas")

    model.mesh.generate(2)
    gmsh.write(args.gmshfile)
    gmsh.finalize()


def solve_field_fem(args):
    mio = meshio.read(args.gmshfile)
    mesh = from_meshio(mio)
    basis = Basis(mesh, ElementTriP3())
    grounded_bnd = basis.get_dofs("domain_bottom")
    needle_bnd = basis.get_dofs("needle")

    # Set epsilon in domain
    eps = np.zeros(mesh.nelements)
    eps[mesh.subdomains["gas"]] = 1.0

    # Set material in domain (0: gas)
    material = np.zeros(mesh.nelements, dtype=int)

    @BilinearForm
    def laplace(u, v, w):
        r = w.x[0]
        return dot(grad(u), grad(v)) * r

    A = laplace.assemble(basis)

    u = basis.zeros()
    u[needle_bnd] = args.V_needle
    u[grounded_bnd] = 0.0

    u = solve(*condense(A, x=u, D=needle_bnd+grounded_bnd))

    # Save output in P1 basis
    basis1 = basis.with_element(ElementTriP1())

    m_out = to_meshio(mesh, encode_cell_data=False)
    m_out.point_data["potential"] = basis1.project(basis.interpolate(u))
    m_out.point_data["E_r"] = basis1.project(-basis.interpolate(u).grad[0])
    m_out.point_data["E_z"] = basis1.project(-basis.interpolate(u).grad[1])
    m_out.point_data["E_norm"] = (m_out.point_data["E_r"]**2 +
                                  m_out.point_data["E_z"]**2)**0.5
    m_out.cell_data["material"] = [material]

    meshio.write(args.output, m_out)
    print(f"Solution written to {args.output}")


if __name__ == '__main__':

    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Build axisymmetric sphere geometry')
    p.add_argument('--V_needle', type=float, default=1.0,
                   help='Potential on the needle [V]')
    p.add_argument('--needle_r', type=float, default=0.25e-3,
                   help='Radius of needle (cylindrical part) [m]')
    p.add_argument('--needle_h', type=float, default=7e-3,
                   help='Height of needle (cylindrical part) [m]')
    p.add_argument('--needle_cone_angle', type=float, default=9.2,
                   help='Half of full angle of conical part (degree)')
    p.add_argument('--needle_tip_r', type=float, default=20e-6,
                   help='Needle tip radius of curvature [m]')
    p.add_argument('--domain_r', type=float, default=10e-3,
                   help='Domain/plate radius [m]')
    p.add_argument('--d_gap', type=float, default=100e-6,
                   help='Radius of top rod')
    p.add_argument('--lc_tip', type=float, default=5e-6,
                   help='Mesh spacing on needle tip [m]')
    p.add_argument('--lc_needle', type=float, default=0.1e-3,
                   help='Mesh spacing rest of needle [m]')
    p.add_argument('--lc_other', type=float, default=0.5e-3,
                   help='Mesh spacing on other domain parts [m]')
    p.add_argument('--gmshfile', type=str, default='needle_plate.msh',
                   help='File containing gmsh mesh')
    p.add_argument('--output', type=str, default='needle_plate.vtu',
                   help='File containing solution')
    args = p.parse_args()

    create_needle_plate_mesh(args)
    solve_field_fem(args)
