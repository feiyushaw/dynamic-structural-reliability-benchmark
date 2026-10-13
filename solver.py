"""
Brake Pedal Structural Dynamics Finite Element Solver
Based on DOLFINx 0.11 and PETSc.

Solves 3D solid linear elasticity and transient elastodynamics (Newmark-beta)
for a brake pedal / bracket cantilever structure under static, rapid pedal loading,
vehicle deceleration inertia, and transient road impacts.
"""

from typing import Dict, Any, Tuple, Optional
import numpy as np
from mpi4py import MPI
import dolfinx
from dolfinx import fem, mesh
from dolfinx.fem import petsc
import ufl
from petsc4py import PETSc


class BrakePedalFEA:
    """
    Finite element analyzer for automotive structural components (brake pedal/bracket).
    """

    def __init__(
        self,
        length: float = 0.3,       # 300 mm
        height: float = 0.04,      # 40 mm (strong axis height)
        thickness: float = 0.01,   # 10 mm (thickness / width)
        nx: int = 60,              # elements along length (5 mm)
        ny: int = 8,               # elements along height (5 mm)
        nz: int = 2,               # elements along thickness (5 mm)
        E: float = 210e9,          # Young's modulus: 210 GPa
        nu: float = 0.3,           # Poisson's ratio: 0.3
        rho: float = 7850.0,       # Density: 7850 kg/m^3
        damping_ratio: float = 0.03,  # Damping ratio: 3%
    ):
        self.L = length
        self.H = height
        self.B = thickness
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.E = E
        self.nu = nu
        self.rho = rho
        self.damping_ratio = damping_ratio

        self.comm = MPI.COMM_WORLD
        self._build_mesh_and_spaces()
        self._setup_boundary_conditions()
        self._setup_variational_forms()
        self._assemble_matrices()
        self._preassemble_load_vectors()

    def _build_mesh_and_spaces(self):
        """Create 3D structured Hex8 mesh and function spaces."""
        self.domain = mesh.create_box(
            self.comm,
            [np.array([0.0, 0.0, 0.0]), np.array([self.L, self.H, self.B])],
            [self.nx, self.ny, self.nz],
            cell_type=dolfinx.mesh.CellType.hexahedron,
        )

        # Displacement vector space: Lagrange P1 (3 components)
        self.V = fem.functionspace(self.domain, ("Lagrange", 1, (3,)))

        # Scalar discontinuous space for von Mises stress evaluation
        self.W = fem.functionspace(self.domain, ("Discontinuous Lagrange", 0))

        # Tip area and component mass
        self.A_tip = self.H * self.B
        self.volume = self.L * self.H * self.B
        self.mass = self.rho * self.volume

    def _setup_boundary_conditions(self):
        """Set up Dirichlet BC: fully clamped at x = 0 (mounting bracket)."""
        fdim = self.domain.topology.dim - 1

        def left_boundary(x):
            return np.isclose(x[0], 0.0)

        fixed_facets = mesh.locate_entities_boundary(self.domain, fdim, left_boundary)
        bc_dofs = fem.locate_dofs_topological(self.V, fdim, fixed_facets)
        self.bc = fem.dirichletbc(np.zeros(3, dtype=np.float64), bc_dofs, self.V)

        # Mark free end facet at x = L (tip surface for pedal load)
        def tip_boundary(x):
            return np.isclose(x[0], self.L)

        tip_facets = mesh.locate_entities_boundary(self.domain, fdim, tip_boundary)
        facet_tags = mesh.meshtags(
            self.domain, fdim, tip_facets, np.full_like(tip_facets, 1, dtype=np.int32)
        )
        self.ds = ufl.Measure("ds", domain=self.domain, subdomain_data=facet_tags)
        self.dx = ufl.Measure("dx", domain=self.domain)

    def _setup_variational_forms(self):
        """Define constitutive relations and bilinear forms."""
        u = ufl.TrialFunction(self.V)
        v = ufl.TestFunction(self.V)

        # Lame parameters
        self.mu = self.E / (2.0 * (1.0 + self.nu))
        self.lmbda = self.E * self.nu / ((1.0 + self.nu) * (1.0 - 2.0 * self.nu))

        def epsilon(w):
            return ufl.sym(ufl.grad(w))

        def sigma(w):
            return 2.0 * self.mu * epsilon(w) + self.lmbda * ufl.tr(epsilon(w)) * ufl.Identity(3)

        self.epsilon = epsilon
        self.sigma = sigma

        # Stiffness and Mass forms
        self.k_form = fem.form(ufl.inner(sigma(u), epsilon(v)) * self.dx)
        self.m_form = fem.form(self.rho * ufl.inner(u, v) * self.dx)

    def _assemble_matrices(self):
        """Assemble K and M matrices, and compute Rayleigh damping matrix C."""
        self.K = petsc.assemble_matrix(self.k_form, bcs=[self.bc])
        self.K.assemble()

        self.M = petsc.assemble_matrix(self.m_form)
        self.M.assemble()

        # Fundamental bending frequency estimate for cantilever beam
        # omega_1 = 3.516 * sqrt(E * I / (rho * A * L^4))
        # I = B * H^3 / 12
        I_strong = (self.B * (self.H ** 3)) / 12.0
        A_sec = self.H * self.B
        omega_1 = 3.516 * np.sqrt((self.E * I_strong) / (self.rho * A_sec * (self.L ** 4)))
        self.omega_1 = omega_1
        self.freq_1 = omega_1 / (2.0 * np.pi)

        # Rayleigh damping: C = alpha * M + beta * K
        alpha_m = self.damping_ratio * omega_1
        beta_k = self.damping_ratio / omega_1

        self.C = self.M.copy()
        self.C.scale(alpha_m)
        self.C.axpy(beta_k, self.K)
        self.C.assemble()

    def _preassemble_load_vectors(self):
        """Pre-assemble unit load vectors for rapid scaling during simulation."""
        v = ufl.TestFunction(self.V)

        # 1. Unit downward force at tip (F_y = -1.0 N applied as surface traction)
        t_unit = fem.Constant(self.domain, np.array([0.0, -1.0 / self.A_tip, 0.0], dtype=np.float64))
        form_tip = fem.form(ufl.inner(t_unit, v) * self.ds(1))
        self.vec_unit_tip = petsc.assemble_vector(form_tip)
        self.vec_unit_tip.ghostUpdate()

        # 2. Unit body force in X direction (1.0 N/m^3)
        b_x = fem.Constant(self.domain, np.array([1.0, 0.0, 0.0], dtype=np.float64))
        form_bx = fem.form(ufl.inner(b_x, v) * self.dx)
        self.vec_unit_bx = petsc.assemble_vector(form_bx)
        self.vec_unit_bx.ghostUpdate()

        # 3. Unit body force in Y direction (1.0 N/m^3)
        b_y = fem.Constant(self.domain, np.array([0.0, 1.0, 0.0], dtype=np.float64))
        form_by = fem.form(ufl.inner(b_y, v) * self.dx)
        self.vec_unit_by = petsc.assemble_vector(form_by)
        self.vec_unit_by.ghostUpdate()

    def compute_von_mises_stress(self, u_func: fem.Function) -> Tuple[fem.Function, float]:
        """Compute von Mises stress field and return its global maximum (in Pa)."""
        sig = self.sigma(u_func)
        s_dev = sig - (1.0 / 3.0) * ufl.tr(sig) * ufl.Identity(3)
        von_Mises = ufl.sqrt(ufl.inner(s_dev, s_dev) * 1.5)

        sigma_expr = fem.Expression(von_Mises, self.W.element.interpolation_points)
        sigma_func = fem.Function(self.W)
        sigma_func.interpolate(sigma_expr)
        max_vm = float(np.max(sigma_func.x.array))
        return sigma_func, max_vm

    def solve_static(
        self,
        F_load: float = 1500.0,
        accel_x: float = 0.0,
        accel_y: float = 0.0,
    ) -> Dict[str, Any]:
        """
        Solve linear static elasticity problem.
        Ku = F_load * vec_tip - rho * accel_x * vec_bx - rho * accel_y * vec_by
        """
        rhs = self.K.createVecRight()
        rhs.zeroEntries()

        # Tip downward force
        if abs(F_load) > 1e-12:
            rhs.axpy(F_load, self.vec_unit_tip)

        # Base acceleration equivalent body force: f_body = -rho * a
        if abs(accel_x) > 1e-12:
            f_bx = -self.rho * accel_x
            rhs.axpy(f_bx, self.vec_unit_bx)
        if abs(accel_y) > 1e-12:
            f_by = -self.rho * accel_y
            rhs.axpy(f_by, self.vec_unit_by)

        petsc.apply_lifting(rhs, [self.k_form], [[self.bc]])
        rhs.ghostUpdate()
        petsc.set_bc(rhs, [self.bc])

        uh = fem.Function(self.V)
        ksp = PETSc.KSP().create(self.comm)
        ksp.setOperators(self.K)
        ksp.setType("preonly")
        ksp.getPC().setType("lu")
        ksp.solve(rhs, uh.x.petsc_vec)
        uh.x.scatter_forward()

        # Displacements
        uy = uh.x.array[1::3]
        max_deflection = float(np.min(uy))  # downward deflection is negative

        # Stress
        _, max_stress = self.compute_von_mises_stress(uh)

        return {
            "u_func": uh,
            "max_deflection_mm": abs(max_deflection) * 1e3,
            "max_stress_MPa": max_stress / 1e6,
            "F_load": F_load,
            "accel_x": accel_x,
            "accel_y": accel_y,
        }

    def solve_transient(
        self,
        t_span: np.ndarray,
        F_t: np.ndarray,
        accel_x_t: Optional[np.ndarray] = None,
        accel_y_t: Optional[np.ndarray] = None,
        record_stress_interval: int = 1,
    ) -> Dict[str, Any]:
        """
        Solve transient elastodynamics using the unconditionally stable implicit Newmark-beta method:
        M u'' + C u' + K u = F_ext(t)
        beta = 0.25, gamma = 0.5 (average acceleration method)
        """
        n_steps = len(t_span)
        dt = t_span[1] - t_span[0]

        if accel_x_t is None:
            accel_x_t = np.zeros(n_steps, dtype=np.float64)
        if accel_y_t is None:
            accel_y_t = np.zeros(n_steps, dtype=np.float64)

        # Newmark coefficients
        beta = 0.25
        gamma = 0.5

        a0 = 1.0 / (beta * dt * dt)
        a1 = gamma / (beta * dt)
        a2 = 1.0 / (beta * dt)
        a3 = 1.0 / (2.0 * beta) - 1.0
        a4 = gamma / beta - 1.0
        a5 = (dt / 2.0) * (gamma / beta - 2.0)
        a6 = dt * (1.0 - gamma)
        a7 = gamma * dt

        # Effective dynamic stiffness: K_hat = K + a0 * M + a1 * C
        K_hat = self.K.copy()
        K_hat.axpy(a0, self.M)
        K_hat.axpy(a1, self.C)
        bc_dofs, _ = self.bc.dof_indices()
        K_hat.zeroRowsColumns(bc_dofs, diag=1.0)
        K_hat.assemble()
        ksp = PETSc.KSP().create(self.comm)
        ksp.setOperators(K_hat)
        ksp.setType("preonly")
        ksp.getPC().setType("lu")
        ksp.setUp()

        # Vectors for state variables
        u_curr = self.K.createVecRight()
        v_curr = self.K.createVecRight()
        a_curr = self.K.createVecRight()
        u_next = self.K.createVecRight()
        v_next = self.K.createVecRight()
        a_next = self.K.createVecRight()

        temp_m = self.K.createVecRight()
        temp_c = self.K.createVecRight()
        rhs = self.K.createVecRight()

        # Zero initial state
        u_curr.set(0.0)
        v_curr.set(0.0)
        a_curr.set(0.0)

        # Function container for stress calculation
        uh = fem.Function(self.V)

        history_time = []
        history_tip_disp = []
        history_max_stress = []

        for step in range(n_steps):
            t = t_span[step]
            F_val = F_t[step]
            ax_val = accel_x_t[step]
            ay_val = accel_y_t[step]

            # 1. Effective dynamic load: F_eff = F_ext(t) + M(...) + C(...)
            # temp_m = a0 * u + a2 * v + a3 * a
            temp_m.zeroEntries()
            temp_m.axpy(a0, u_curr)
            temp_m.axpy(a2, v_curr)
            temp_m.axpy(a3, a_curr)

            # temp_c = a1 * u + a4 * v + a5 * a
            temp_c.zeroEntries()
            temp_c.axpy(a1, u_curr)
            temp_c.axpy(a4, v_curr)
            temp_c.axpy(a5, a_curr)

            # Assemble external force
            rhs.zeroEntries()
            if abs(F_val) > 1e-12:
                rhs.axpy(F_val, self.vec_unit_tip)
            if abs(ax_val) > 1e-12:
                rhs.axpy(-self.rho * ax_val, self.vec_unit_bx)
            if abs(ay_val) > 1e-12:
                rhs.axpy(-self.rho * ay_val, self.vec_unit_by)

            # Add inertia and damping history forces
            self.M.multAdd(temp_m, rhs, rhs)
            self.C.multAdd(temp_c, rhs, rhs)

            # Apply homogeneous Dirichlet BC
            petsc.set_bc(rhs, [self.bc])

            # Solve linear system
            ksp.solve(rhs, u_next)

            # 2. Update acceleration and velocity
            # a_next = a0 * (u_next - u_curr) - a2 * v_curr - a3 * a_curr
            a_next.zeroEntries()
            a_next.axpy(a0, u_next)
            a_next.axpy(-a0, u_curr)
            a_next.axpy(-a2, v_curr)
            a_next.axpy(-a3, a_curr)

            # v_next = v_curr + a6 * a_curr + a7 * a_next
            v_next.zeroEntries()
            v_next.axpy(1.0, v_curr)
            v_next.axpy(a6, a_curr)
            v_next.axpy(a7, a_next)

            # State propagation
            u_next.copy(u_curr)
            v_next.copy(v_curr)
            a_next.copy(a_curr)

            # Record tip vertical displacement
            u_arr = u_curr.getArray()
            # Tip deflection: minimum y-displacement
            tip_disp_mm = abs(float(np.min(u_arr[1::3]))) * 1e3
            history_time.append(t)
            history_tip_disp.append(tip_disp_mm)

            # Calculate stress field periodically
            if step % record_stress_interval == 0 or step == n_steps - 1:
                u_curr.copy(uh.x.petsc_vec)
                uh.x.scatter_forward()
                _, max_vm = self.compute_von_mises_stress(uh)
                history_max_stress.append(max_vm / 1e6)
            else:
                # Approximate with last recorded stress
                history_max_stress.append(history_max_stress[-1])

        return {
            "time": np.array(history_time),
            "tip_deflection_mm": np.array(history_tip_disp),
            "max_stress_MPa": np.array(history_max_stress),
            "peak_stress_MPa": float(np.max(history_max_stress)),
            "peak_deflection_mm": float(np.max(history_tip_disp)),
        }
