#!/usr/bin/env python3

"""Create 2D Cartesian mesh for wire-to-wire geometry and solve electrostatic
problem
"""

import numpy as np
import gmsh
from skfem import Basis, BilinearForm, solve, condense
from skfem import ElementTriP0, ElementTriP1, ElementTriP3
from skfem.io.meshio import from_meshio, to_meshio
from skfem.helpers import dot, grad
import meshio
import argparse


def create_w2w_mesh(args):
    gmsh.initialize()

    model = gmsh.model
    geo = model.geo

    model.add("wire_to_wire")

    Lx = np.array(args.L_domain)
    R0 = args.r_wire
    R1 = args.r_wire + args.d_diel

    # Centers of wires
    c1 = np.array([-R1 - 0.5 * args.d_wires, 0., 0.])
    c2 = np.array([R1 + 0.5 * args.d_wires, 0., 0.])

    # Domain rectangle
    p1 = geo.addPoint(-0.5 * Lx[0], -0.5 * Lx[1], 0., args.lc_outer)
    p2 = geo.addPoint(0.5 * Lx[0], -0.5 * Lx[1], 0., args.lc_outer)
    p3 = geo.addPoint(0.5 * Lx[0], 0.5 * Lx[1], 0., args.lc_outer)
    p4 = geo.addPoint(-0.5 * Lx[0], 0.5 * Lx[1], 0., args.lc_outer)

    l1 = geo.addLine(p1, p2)
    l2 = geo.addLine(p2, p3)
    l3 = geo.addLine(p3, p4)
    l4 = geo.addLine(p4, p1)

    domain_loop = geo.addCurveLoop([l1, l2, l3, l4])

    # Offsets used to define circles as segments of 4 arcs
    offsets = np.array([[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]])

    # Wire 1
    pc1 = geo.addPoint(*c1)
    po1 = np.zeros(4, dtype=int)  # Outer points
    pi1 = np.zeros(4, dtype=int)  # Inner points
    ca1 = np.zeros(4, dtype=int)  # Outers arc loops
    ci1 = np.zeros(4, dtype=int)  # Inner arc loops

    # Wire 2
    pc2 = geo.addPoint(*c2)
    po2 = np.zeros(4, dtype=int)  # Outer points
    pi2 = np.zeros(4, dtype=int)  # Inner points
    ca2 = np.zeros(4, dtype=int)  # Outers arc loops
    ci2 = np.zeros(4, dtype=int)  # Inner arc loops

    for i, di in enumerate(offsets):
        po1[i] = geo.addPoint(*(c1 + R1 * di), args.lc_wire)
        pi1[i] = geo.addPoint(*(c1 + R0 * di), args.lc_wire)

        if i == 2 and args.d_wires <= 0:
            # Avoid duplicate point when two wires touch
            po2[i] = po1[0]
        else:
            po2[i] = geo.addPoint(*(c2 + R1 * di), args.lc_wire)
        pi2[i] = geo.addPoint(*(c2 + R0 * di), args.lc_wire)

    for i in range(4):
        ca1[i] = geo.addCircleArc(po1[i], pc1, po1[(i+1) % 4])
        ci1[i] = geo.addCircleArc(pi1[i], pc1, pi1[(i+1) % 4])

        ca2[i] = geo.addCircleArc(po2[i], pc2, po2[(i+1) % 4])
        ci2[i] = geo.addCircleArc(pi2[i], pc2, pi2[(i+1) % 4])

    outer_loop1 = geo.addCurveLoop(ca1)
    inner_loop1 = geo.addCurveLoop(ci1)

    outer_loop2 = geo.addCurveLoop(ca2)
    inner_loop2 = geo.addCurveLoop(ci2)

    # Dielectric rings (surface between loops)
    diel1 = geo.addPlaneSurface([outer_loop1, inner_loop1])
    diel2 = geo.addPlaneSurface([outer_loop2, inner_loop2])

    # Gas domain (rectangle with both wires cut out)
    gas = geo.addPlaneSurface([domain_loop, outer_loop1, outer_loop2])

    geo.synchronize()

    # Physical groups
    model.addPhysicalGroup(1, ci1, name="wire1_bnd")
    model.addPhysicalGroup(1, ci2, name="wire2_bnd")
    model.addPhysicalGroup(1, [l1, l2, l3, l4], name="outer_bnd")
    model.addPhysicalGroup(2, [diel1], name="diel1")
    model.addPhysicalGroup(2, [diel2], name="diel2")
    model.addPhysicalGroup(2, [gas], name="gas")

    # Build mesh
    mesh = model.mesh
    mesh.generate(2)

    gmsh.write(args.gmshfile)
    gmsh.finalize()


def solve_field_fem(args):
    mio = meshio.read(args.gmshfile)
    mesh = from_meshio(mio)
    basis = Basis(mesh, ElementTriP3())
    wire1_bnd = basis.get_dofs("wire1_bnd")
    wire2_bnd = basis.get_dofs("wire2_bnd")
    outer_bnd = basis.get_dofs("outer_bnd")

    # Set epsilon in domain
    eps = np.zeros(mesh.nelements)
    eps[mesh.subdomains["gas"]] = 1.0
    eps[mesh.subdomains["diel1"]] = args.eps_diel1
    eps[mesh.subdomains["diel2"]] = args.eps_diel2

    # Set material in domain (0: gas)
    material = np.zeros(mesh.nelements, dtype=int)
    material[mesh.subdomains["diel1"]] = 1
    material[mesh.subdomains["diel2"]] = 2

    # To use epsilon in Laplace operator
    basis0 = basis.with_element(ElementTriP0())
    eps_field = basis0.interpolate(eps)

    @BilinearForm
    def laplace(u, v, w):
        return w.eps * dot(grad(u), grad(v))

    A = laplace.assemble(basis, eps=eps_field)

    # Set boundary conditions
    u = basis.zeros()
    u[wire1_bnd] = args.V_wire1
    u[wire2_bnd] = args.V_wire2
    u[outer_bnd] = args.V_outer

    u = solve(*condense(A, x=u, D=wire1_bnd+wire2_bnd+outer_bnd))

    # Save output in P1 basis
    basis1 = basis.with_element(ElementTriP1())

    # Compute field only in the gas solution; it will be zero elsewhere
    gas_basis = Basis(mesh, ElementTriP3(), elements=mesh.subdomains["gas"])
    gas_basis1 = gas_basis.with_element(ElementTriP1())
    uh_gas = gas_basis.interpolate(u)
    Ex_gas = gas_basis1.project(-uh_gas.grad[0])
    Ey_gas = gas_basis1.project(-uh_gas.grad[1])

    m_out = to_meshio(mesh, encode_cell_data=False)
    m_out.cell_data["material"] = [material]
    m_out.cell_data["eps"] = [eps]
    m_out.point_data["potential"] = basis1.project(basis.interpolate(u))
    m_out.point_data["E_x_gas"] = Ex_gas
    m_out.point_data["E_y_gas"] = Ey_gas
    m_out.point_data["E_norm_gas"] = (Ex_gas**2 + Ey_gas**2)**0.5

    meshio.write(args.output, m_out)
    print(f"Solution written to {args.output}")


if __name__ == '__main__':

    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Build axisymmetric sphere geometry')
    p.add_argument('--V_wire1', type=float, default=1.0,
                   help='Potential on wire 1 [V]')
    p.add_argument('--V_wire2', type=float, default=-1.0,
                   help='Potential on wire 2 [V]')
    p.add_argument('--V_outer', type=float, default=0.1,
                   help='Potential on outer electrode [V]')
    p.add_argument('--eps_diel1', type=float, default=5.0,
                   help='Relative permittivity of wire1 dielectric')
    p.add_argument('--eps_diel2', type=float, default=5.0,
                   help='Relative permittivity of wire2 dielectric')
    p.add_argument('--eps_r', type=float, default=1.0,
                   help='Relative permittivity')
    p.add_argument('--L_domain', type=float, nargs=2, default=[1.0, 1.0],
                   help='Domain size (m)')
    p.add_argument('--r_wire', type=float, default=0.1,
                   help='Radius of wire')
    p.add_argument('--d_diel', type=float, default=0.05,
                   help='Thickness of dielectric')
    p.add_argument('--d_wires', type=float, default=0.0,
                   help='Spacing between wires')
    p.add_argument('--lc_wire', type=float, default=0.01,
                   help='Mesh spacing on wires')
    p.add_argument('--lc_outer', type=float, default=0.05,
                   help='Mesh spacing on outer boundary')
    p.add_argument('--gmshfile', type=str, default='wire_to_wire.msh',
                   help='File containing gmsh mesh')
    p.add_argument('--output', type=str, default='wire_to_wire.vtu',
                   help='File containing solution')
    args = p.parse_args()

    create_w2w_mesh(args)
    solve_field_fem(args)
