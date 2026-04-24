#!/usr/bin/env python3

"""Create axisymmetric mesh for inner sphere and outer sphere and solve
electrostatic problem
"""

import gmsh
from skfem import Basis, BilinearForm, solve, condense
from skfem import ElementTriP1, ElementTriP3
from skfem.io.meshio import from_meshio, to_meshio
import meshio
import argparse


def create_axi_sphere_mesh(args):
    gmsh.initialize()

    model = gmsh.model
    geo = model.geo

    model.add("sphere")

    # Center
    p0 = geo.addPoint(0, 0, 0, args.lc_sphere)

    # Sphere semicircle points
    p1 = geo.addPoint(0, args.r_sphere, 0, args.lc_sphere)
    p2 = geo.addPoint(args.r_sphere, 0, 0, args.lc_sphere)
    p3 = geo.addPoint(0, -args.r_sphere, 0, args.lc_sphere)

    # Outer boundary semicircle points
    p4 = geo.addPoint(0, args.r_outer, 0, args.lc_outer)
    p5 = geo.addPoint(args.r_outer, 0, 0, args.lc_outer)
    p6 = geo.addPoint(0, -args.r_outer, 0, args.lc_outer)

    # Sphere arcs (inner boundary)
    c1 = geo.addCircleArc(p1, p0, p2)
    c2 = geo.addCircleArc(p2, p0, p3)

    # Outer arcs (far-field boundary)
    c3 = geo.addCircleArc(p4, p0, p5)
    c4 = geo.addCircleArc(p5, p0, p6)

    # Axis segments (r = 0)
    l_axis_top = geo.addLine(p4, p1)
    l_axis_bot = geo.addLine(p3, p6)

    # Curve loop around the domain
    cl = geo.addCurveLoop([l_axis_top, c1, c2, l_axis_bot, -c4, -c3])
    surf = geo.addPlaneSurface([cl])

    geo.synchronize()

    # Physical groups
    model.addPhysicalGroup(1, [c1, c2], tag=1000, name="Sphere")
    model.addPhysicalGroup(1, [c3, c4], tag=2000, name="Outer")
    model.addPhysicalGroup(1, [l_axis_top, l_axis_bot], tag=3000,
                           name="Axis")
    model.addPhysicalGroup(2, [surf], tag=4000, name="Domain")

    # Build mesh
    mesh = model.mesh
    mesh.generate(2)

    gmsh.write(args.gmshfile)
    gmsh.finalize()


def solve_field_fem(args):
    mio = meshio.read(args.gmshfile)
    mesh = from_meshio(mio)
    basis = Basis(mesh, ElementTriP3())
    sphere_bnd = basis.get_dofs("Sphere")
    outer_bnd = basis.get_dofs("Outer")

    @BilinearForm
    def laplace_axi(u, v, w):
        r = w.x[0]
        return args.eps_r * (u.grad[0] * v.grad[0] + u.grad[1] * v.grad[1]) * r

    A = laplace_axi.assemble(basis)

    u = basis.zeros()
    u[sphere_bnd] = args.V_sphere
    u[outer_bnd] = args.V_outer

    u = solve(*condense(A, x=u, D=sphere_bnd+outer_bnd))

    # Save output in P1 basis
    basis1 = basis.with_element(ElementTriP1())

    m_out = to_meshio(mesh, encode_cell_data=False)
    m_out.point_data["potential"] = basis1.project(basis.interpolate(u))
    m_out.point_data["E_r"] = basis1.project(-basis.interpolate(u).grad[0])
    m_out.point_data["E_z"] = basis1.project(-basis.interpolate(u).grad[1])
    m_out.point_data["E_norm"] = (m_out.point_data["E_r"]**2 +
                                  m_out.point_data["E_z"]**2)**0.5

    meshio.write(args.output, m_out)
    print(f"Solution written to {args.output}")


if __name__ == '__main__':

    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Build axisymmetric sphere geometry')
    p.add_argument('--V_sphere', type=float, default=1.0,
                   help='Potential on inner electrode [V]')
    p.add_argument('--V_outer', type=float, default=0.1,
                   help='Potential on outer electrode [V]')
    p.add_argument('--eps_r', type=float, default=1.0,
                   help='Relative permittivity')
    p.add_argument('--r_sphere', type=float, default=1.0,
                   help='Radius of inner sphere')
    p.add_argument('--r_outer', type=float, default=10.0,
                   help='Radius of outer sphere')
    p.add_argument('--lc_sphere', type=float, default=0.05,
                   help='Mesh spacing on inner sphere')
    p.add_argument('--lc_outer', type=float, default=0.5,
                   help='Mesh spacing on outer sphere')
    p.add_argument('--gmshfile', type=str, default='sphere.msh',
                   help='File containing gmsh mesh')
    p.add_argument('--output', type=str, default='sphere.vtu',
                   help='File containing solution')
    args = p.parse_args()

    create_axi_sphere_mesh(args)
    solve_field_fem(args)
