"""
Common Finite Element definitions and utilities for Cantilever Steel Beam.
Implements 3D elasticity forms, system matrices assembly, and von Mises stress evaluation.
"""

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import numpy as np
from mpi4py import MPI
import dolfinx
from dolfinx import fem, mesh
from dolfinx.fem import petsc
import ufl
from petsc4py import PETSc
from mesh.generate_mesh import create_cantilever_mesh


class BeamModel:
    """Minimal Cantilever Steel Beam FEA Model."""

    def __init__(self, L=0.3, H=0.04, B=0.01, E=210e9, nu=0.3, rho=7850.0, damping_ratio=0.03):
        self.L, self.H, self.B = L, H, B
        self.E, self.nu, self.rho = E, nu, rho
        self.damping_ratio = damping_ratio
        self.comm = MPI.COMM_WORLD

        # Mesh and Function Spaces
        self.domain = create_cantilever_mesh(self.comm, L, H, B)
        self.V = fem.functionspace(self.domain, ("Lagrange", 1, (3,)))
        self.W = fem.functionspace(self.domain, ("Discontinuous Lagrange", 0))

        # Tip surface area and component mass
        self.A_tip = H * B
        self.mass = rho * (L * H * B)

        # Boundary Conditions: clamped at x = 0
        fdim = self.domain.topology.dim - 1
        fixed_facets = mesh.locate_entities_boundary(self.domain, fdim, lambda x: np.isclose(x[0], 0.0))
        self.bc = fem.dirichletbc(np.zeros(3, dtype=np.float64), fem.locate_dofs_topological(self.V, fdim, fixed_facets), self.V)

        # Tip boundary at x = L
        tip_facets = mesh.locate_entities_boundary(self.domain, fdim, lambda x: np.isclose(x[0], self.L))
        facet_tags = mesh.meshtags(self.domain, fdim, tip_facets, np.full_like(tip_facets, 1, dtype=np.int32))
        self.ds = ufl.Measure("ds", domain=self.domain, subdomain_data=facet_tags)
        self.dx = ufl.Measure("dx", domain=self.domain)

        # Bilinear Forms
        u = ufl.TrialFunction(self.V)
        v = ufl.TestFunction(self.V)
        mu = E / (2.0 * (1.0 + nu))
        lmbda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
        eps = lambda w: ufl.sym(ufl.grad(w))
        sig = lambda w: 2.0 * mu * eps(w) + lmbda * ufl.tr(eps(w)) * ufl.Identity(3)
        self.sigma = sig

        self.k_form = fem.form(ufl.inner(sig(u), eps(v)) * self.dx)
        self.m_form = fem.form(rho * ufl.inner(u, v) * self.dx)

        # Assemble Matrices
        self.K = petsc.assemble_matrix(self.k_form, bcs=[self.bc])
        self.K.assemble()
        self.M = petsc.assemble_matrix(self.m_form)
        self.M.assemble()

        # Rayleigh Damping: C = alpha * M + beta * K
        I_strong = (B * (H ** 3)) / 12.0
        omega_1 = 3.516 * np.sqrt((E * I_strong) / (rho * (H * B) * (L ** 4)))
        self.omega_1 = omega_1
        self.freq_1 = omega_1 / (2.0 * np.pi)

        alpha_m = damping_ratio * omega_1
        beta_k = damping_ratio / omega_1
        self.C = self.M.copy()
        self.C.scale(alpha_m)
        self.C.axpy(beta_k, self.K)
        self.C.assemble()

        # Pre-assemble unit load vectors: tip force (downward Y) and body force (X direction)
        t_unit = fem.Constant(self.domain, np.array([0.0, -1.0 / self.A_tip, 0.0], dtype=np.float64))
        self.vec_tip_y = petsc.assemble_vector(fem.form(ufl.inner(t_unit, v) * self.ds(1)))
        self.vec_tip_y.ghostUpdate()

        b_x = fem.Constant(self.domain, np.array([1.0, 0.0, 0.0], dtype=np.float64))
        self.vec_body_x = petsc.assemble_vector(fem.form(ufl.inner(b_x, v) * self.dx))
        self.vec_body_x.ghostUpdate()

    def compute_von_mises(self, uh: fem.Function) -> float:
        """Compute max von Mises stress in Pa."""
        sig = self.sigma(uh)
        s_dev = sig - (1.0 / 3.0) * ufl.tr(sig) * ufl.Identity(3)
        vm = ufl.sqrt(ufl.inner(s_dev, s_dev) * 1.5)
        expr = fem.Expression(vm, self.W.element.interpolation_points)
        f_vm = fem.Function(self.W)
        f_vm.interpolate(expr)
        return float(np.max(f_vm.x.array))
