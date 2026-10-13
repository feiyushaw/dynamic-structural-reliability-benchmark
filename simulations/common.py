"""
Common Finite Element definitions and rate-dependent constitutive utilities
for Cantilever Steel Beam under dynamic impact loading.
Implements:
  - 3D elastodynamics forms and system matrices (Newmark-beta)
  - Strain rate tensor and equivalent strain rate calculation
  - Cowper-Symonds rate-dependent dynamic yield strength:
      sigma_y(edot) = sigma_0 * [1 + ((edot - edot_0) / D)^(1/p)] (for edot > edot_0)
  - Dynamic ductility reduction / embrittlement evaluation
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


# Cowper-Symonds material constants for structural steel
COWPER_SYMONDS_SIGMA0 = 250e6  # 250 MPa quasi-static yield stress
COWPER_SYMONDS_D = 40.4        # D parameter: 40.4 s^-1
COWPER_SYMONDS_P = 5.0         # p parameter: 5.0
STATIC_STRAIN_RATE = 1.0e-4    # Quasi-static reference rate: 10^-4 s^-1


def cowper_symonds_yield(
    edot: float,
    sigma0: float = COWPER_SYMONDS_SIGMA0,
    D: float = COWPER_SYMONDS_D,
    p: float = COWPER_SYMONDS_P,
    edot_0: float = STATIC_STRAIN_RATE,
) -> float:
    """
    Compute dynamic yield stress (Pa) via Cowper-Symonds rate-dependent constitutive law.
    Normalized such that at quasi-static rate edot <= edot_0, sigma_y = sigma0 exactly.
    """
    edot_val = float(edot)
    if edot_val <= edot_0:
        return sigma0
    return sigma0 * (1.0 + ((edot_val - edot_0) / D) ** (1.0 / p))


def dynamic_ductility_ratio(edot: float, C_eps: float = 0.08, edot_0: float = STATIC_STRAIN_RATE) -> float:
    """
    Compute relative dynamic failure ductility ratio eps_f(edot) / eps_f0.
    Reflects dynamic embrittlement under high strain rates.
    """
    edot_val = max(float(edot), edot_0)
    rate_factor = 1.0 + C_eps * np.log(edot_val / edot_0)
    return 1.0 / rate_factor


class BeamModel:
    """Cantilever Steel Beam FEA Model with Rate-Dependent Constitutive Evaluation."""

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
        self.epsilon = eps
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

    def compute_strain_rate(self, v_func: fem.Function) -> float:
        """
        Compute maximum equivalent strain rate (s^-1) from velocity field v_func:
        edot_eq = sqrt(2/3 * eps_dot : eps_dot)
        """
        eps_dot = self.epsilon(v_func)
        edot_dev = eps_dot - (1.0 / 3.0) * ufl.tr(eps_dot) * ufl.Identity(3)
        edot_eq = ufl.sqrt(ufl.inner(edot_dev, edot_dev) * (2.0 / 3.0) + 1e-16)
        expr = fem.Expression(edot_eq, self.W.element.interpolation_points)
        f_edot = fem.Function(self.W)
        f_edot.interpolate(expr)
        return float(np.max(f_edot.x.array))
