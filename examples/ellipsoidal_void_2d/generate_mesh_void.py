#!/usr/bin/env python3

"""Create axisymmetric mesh for a dielectric containing an ellipsoidal void
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
    model.add("ellipsoidal_void")

    domain_r = args.domain_r
    domain_h = args.domain_h
    z = args.ellips_z
    h = args.ellips_h
    r = args.ellips_r

    # --------------------------------------------------------------------
    # POINTS
    # --------------------------------------------------------------------
    p_bot_left = geo.addPoint(0.0, 0.0, 0.0, args.lc_domain)
    p_bot_right = geo.addPoint(domain_r, 0.0, 0.0, args.lc_domain)
    p_top_left = geo.addPoint(0.0, domain_h, 0.0, args.lc_domain)
    p_top_right = geo.addPoint(domain_r, domain_h, 0.0, args.lc_domain)

    ellips_cnt = geo.addPoint(0.0, z, 0.0, args.lc_ellips)
    ellips_bot = geo.addPoint(0.0, z-0.5*h, 0.0, args.lc_ellips)
    ellips_top = geo.addPoint(0.0, z+0.5*h, 0.0, args.lc_ellips)
    ellips_mid = geo.addPoint(r, z, 0.0, args.lc_ellips)

    # --------------------------------------------------------------------
    # LINES AND ARCS
    # --------------------------------------------------------------------
    l_ellips_bot_mid = geo.addEllipseArc(
        ellips_bot, ellips_cnt, ellips_top, ellips_mid)
    l_ellips_mid_top = geo.addEllipseArc(
        ellips_mid, ellips_cnt, ellips_top, ellips_top)
    l_ellips_axis = geo.addLine(ellips_top, ellips_bot)
    l_ellips = [l_ellips_axis, l_ellips_bot_mid, l_ellips_mid_top]

    cl_ellips = geo.addCurveLoop(l_ellips)

    l_domain_bot = geo.addLine(p_bot_left, p_bot_right)
    l_domain_side = geo.addLine(p_bot_right, p_top_right)
    l_domain_top = geo.addLine(p_top_right, p_top_left)

    # Outer boundary loop, counter-clockwise
    cl_outer = geo.addCurveLoop([
        l_domain_bot,
        l_domain_side,
        l_domain_top,
        geo.addLine(p_top_left, ellips_top),
        l_ellips_axis,
        geo.addLine(ellips_bot, p_bot_left),
    ])

    s_gas = geo.addPlaneSurface([cl_ellips])
    s_dielectric = geo.addPlaneSurface([cl_outer, cl_ellips])

    geo.synchronize()

    # --------------------------------------------------------------------
    # PHYSICAL GROUPS
    # --------------------------------------------------------------------
    model.addPhysicalGroup(1, [l_domain_top], name="domain_top")
    model.addPhysicalGroup(1, [l_domain_bot], name="domain_bottom")
    model.addPhysicalGroup(1, [l_domain_side], name="domain_side")
    model.addPhysicalGroup(2, [s_gas], name="gas")
    model.addPhysicalGroup(2, [s_dielectric], name="dielectric")

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
    bot_bnd = basis.get_dofs("domain_bottom")
    top_bnd = basis.get_dofs("domain_top")

    # Set epsilon in domain
    eps = np.ones(mesh.nelements) * args.eps_r
    eps[mesh.subdomains["gas"]] = 1.0

    # Set material in domain (0: gas)
    material = np.zeros(mesh.nelements, dtype=int)
    material[mesh.subdomains["dielectric"]] = 1

    # To use epsilon in Laplace operator
    basis0 = basis.with_element(ElementTriP0())
    eps_field = basis0.interpolate(eps)

    @BilinearForm
    def laplace(u, v, w):
        r = w.x[0]
        return w.eps * dot(grad(u), grad(v)) * r

    A = laplace.assemble(basis, eps=eps_field)

    u = basis.zeros()
    u[top_bnd] = args.V_top
    u[bot_bnd] = 0.0

    u = solve(*condense(A, x=u, D=top_bnd+bot_bnd))

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
    p.add_argument('--eps_r', type=float, default=5.0,
                   help='Relative permittivity of dielectric')
    p.add_argument('--domain_r', type=float, default=10.0e-3,
                   help='Radius of domain')
    p.add_argument('--domain_h', type=float, default=20.0e-3,
                   help='Height of domain')
    p.add_argument('--lc_domain', type=float, default=1e-3,
                   help='Mesh spacing on domain')
    p.add_argument('--lc_ellips', type=float, default=0.1e-3,
                   help='Mesh spacing on ellipsoid')
    p.add_argument('--ellips_z', type=float, default=10.0e-3,
                   help='z-coordinate of ellipsoid center')
    p.add_argument('--ellips_h', type=float, default=2.0e-3,
                   help='Height of ellipsoid')
    p.add_argument('--ellips_r', type=float, default=0.5e-3,
                   help='Radius of ellipsoid')
    p.add_argument('--gmshfile', type=str, default='void.msh',
                   help='File containing gmsh mesh')
    p.add_argument('--output', type=str, default='void.vtu',
                   help='File containing solution')
    args = p.parse_args()

    create_vessel_mesh(args)
    solve_field_fem(args)
