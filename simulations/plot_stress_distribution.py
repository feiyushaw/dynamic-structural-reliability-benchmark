"""
Figure 3: 2D Spatial Stress Contours & Failure Risk Maps.
Generates results/fig3_spatial_stress_contours.png and ParaView 3D VTK exports.
Clean scientific layout with zero text clutter.
"""

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import ufl
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from dolfinx import mesh, fem, io
from dolfinx.fem import petsc
from petsc4py import PETSc
from simulations.common import BeamModel

mpl.rcParams["font.family"] = ["Noto Sans CJK SC", "DejaVu Sans", "sans-serif"]
mpl.rcParams["axes.unicode_minus"] = False
mpl.rcParams["figure.dpi"] = 300
mpl.rcParams["savefig.dpi"] = 300
mpl.rcParams["font.size"] = 9.5
mpl.rcParams["lines.linewidth"] = 1.6


def extract_stress_grid(model: BeamModel, uh: fem.Function):
    sig = model.sigma(uh)
    s_dev = sig - (1.0 / 3.0) * ufl.tr(sig) * ufl.Identity(3)
    vm = ufl.sqrt(ufl.inner(s_dev, s_dev) * 1.5)
    expr = fem.Expression(vm, model.W.element.interpolation_points)
    f_vm = fem.Function(model.W); f_vm.interpolate(expr)
    stress_arr = f_vm.x.array / 1e6

    num_cells = model.domain.topology.index_map(3).size_local
    midpoints = mesh.compute_midpoints(model.domain, 3, np.arange(num_cells, dtype=np.int32))
    x_coords = np.unique(np.round(midpoints[:, 0], 4))
    y_coords = np.unique(np.round(midpoints[:, 1], 4))

    grid = np.zeros((len(y_coords), len(x_coords)))
    counts = np.zeros((len(y_coords), len(x_coords)))
    for c in range(num_cells):
        ix = np.argmin(np.abs(x_coords - midpoints[c, 0]))
        iy = np.argmin(np.abs(y_coords - midpoints[c, 1]))
        grid[iy, ix] += stress_arr[c]
        counts[iy, ix] += 1
    grid /= counts
    return x_coords * 1000.0, y_coords * 1000.0, grid, f_vm


def run_stress_distribution_analysis(out_dir: str = "results"):
    os.makedirs(out_dir, exist_ok=True)
    model = BeamModel()

    # Case 1 Solve
    rhs_s1 = model.K.createVecRight(); rhs_s1.zeroEntries(); rhs_s1.axpy(1500.0, model.vec_tip_y)
    petsc.apply_lifting(rhs_s1, [model.k_form], [[model.bc]])
    rhs_s1.ghostUpdate(); petsc.set_bc(rhs_s1, [model.bc])
    uh_s1 = fem.Function(model.V)
    ksp_s1 = PETSc.KSP().create(model.comm); ksp_s1.setOperators(model.K); ksp_s1.setType("preonly"); ksp_s1.getPC().setType("lu")
    ksp_s1.solve(rhs_s1, uh_s1.x.petsc_vec); uh_s1.x.scatter_forward()
    x_mm, y_mm, grid_s1, f_vm_s1 = extract_stress_grid(model, uh_s1)

    # Case 3 Dynamic Solve
    dt = 0.0002; t_span = np.arange(0.0, 0.025 + dt, dt); n_steps = len(t_span); g = 9.80665
    beta, gamma = 0.25, 0.5
    a0 = 1.0 / (beta * dt * dt); a1 = gamma / (beta * dt); a2 = 1.0 / (beta * dt)
    a3 = 1.0 / (2.0 * beta) - 1.0; a4 = gamma / beta - 1.0; a5 = (dt / 2.0) * (gamma / beta - 2.0)
    a6 = dt * (1.0 - gamma); a7 = gamma * dt

    K_hat = model.K.copy(); K_hat.axpy(a0, model.M); K_hat.axpy(a1, model.C)
    bc_dofs, _ = model.bc.dof_indices(); K_hat.zeroRowsColumns(bc_dofs, diag=1.0); K_hat.assemble()
    ksp_dyn = PETSc.KSP().create(model.comm); ksp_dyn.setOperators(K_hat); ksp_dyn.setType("preonly"); ksp_dyn.getPC().setType("lu"); ksp_dyn.setUp()

    u_curr = model.K.createVecRight(); v_curr = model.K.createVecRight(); a_curr = model.K.createVecRight()
    u_next = model.K.createVecRight(); v_next = model.K.createVecRight(); a_next = model.K.createVecRight()
    temp_m = model.K.createVecRight(); temp_c = model.K.createVecRight(); rhs = model.K.createVecRight()
    u_curr.set(0.0); v_curr.set(0.0); a_curr.set(0.0)

    uh_peak = fem.Function(model.V); max_s = 0.0
    for step in range(n_steps):
        t = t_span[step]
        F_val = 2250.0 * (t / 0.010) if t < 0.010 else 2250.0 * (1.0 - (t - 0.010) / 0.015) + 1500.0 * ((t - 0.010) / 0.015)
        ax_val = -1.0 * g if t >= 0.010 else 0.0

        temp_m.zeroEntries(); temp_m.axpy(a0, u_curr); temp_m.axpy(a2, v_curr); temp_m.axpy(a3, a_curr)
        temp_c.zeroEntries(); temp_c.axpy(a1, u_curr); temp_c.axpy(a4, v_curr); temp_c.axpy(a5, a_curr)
        rhs.zeroEntries()
        if abs(F_val) > 1e-12: rhs.axpy(F_val, model.vec_tip_y)
        if abs(ax_val) > 1e-12: rhs.axpy(-model.rho * ax_val, model.vec_body_x)
        model.M.multAdd(temp_m, rhs, rhs); model.C.multAdd(temp_c, rhs, rhs)
        petsc.set_bc(rhs, [model.bc]); ksp_dyn.solve(rhs, u_next)

        a_next.zeroEntries(); a_next.axpy(a0, u_next); a_next.axpy(-a0, u_curr); a_next.axpy(-a2, v_curr); a_next.axpy(-a3, a_curr)
        v_next.zeroEntries(); v_next.axpy(1.0, v_curr); v_next.axpy(a6, a_curr); v_next.axpy(a7, a_next)
        u_next.copy(u_curr); v_next.copy(v_curr); a_next.copy(a_curr)

        if t >= 0.010:
            uh_t = fem.Function(model.V); u_curr.copy(uh_t.x.petsc_vec); uh_t.x.scatter_forward()
            vm = model.compute_von_mises(uh_t)
            if vm > max_s: max_s = vm; u_curr.copy(uh_peak.x.petsc_vec); uh_peak.x.scatter_forward()

    _, _, grid_c3, f_vm_c3 = extract_stress_grid(model, uh_peak)

    # Export VTK
    with io.VTKFile(model.comm, os.path.join(out_dir, "beam_case1_static.pvd"), "w") as vtk: vtk.write_function([uh_s1, f_vm_s1])
    with io.VTKFile(model.comm, os.path.join(out_dir, "beam_case3_impact.pvd"), "w") as vtk: vtk.write_function([uh_peak, f_vm_c3])

    risk_s1 = grid_s1 / 400.0  # static capacity = 400 MPa
    risk_c3 = grid_c3 / 245.0  # dynamic capacity = 245 MPa

    # Render Minimalist 4-Panel Contour Plot
    X_mesh, Y_mesh = np.meshgrid(x_mm, y_mm)
    fig, axes = plt.subplots(4, 1, figsize=(11.5, 8.2), sharex=True)

    # (a) Case 1 Stress
    im1 = axes[0].pcolormesh(X_mesh, Y_mesh, grid_s1, cmap="turbo", shading="gouraud", vmin=0, vmax=230)
    cb1 = fig.colorbar(im1, ax=axes[0], pad=0.015, aspect=12); cb1.set_label("应力 (MPa)", fontsize=8.5)
    axes[0].set_title("(a) Case 1 静态应力场 $\\sigma_{vM}$ (峰值 143.6 MPa)", fontsize=9.5, fontweight="bold", pad=4)
    axes[0].axhline(20.0, color="white", linestyle=":", lw=1.0, alpha=0.5)

    # (b) Case 3 Stress
    im2 = axes[1].pcolormesh(X_mesh, Y_mesh, grid_c3, cmap="turbo", shading="gouraud", vmin=0, vmax=230)
    cb2 = fig.colorbar(im2, ax=axes[1], pad=0.015, aspect=12); cb2.set_label("应力 (MPa)", fontsize=8.5)
    axes[1].set_title("(b) Case 3 动态冲击应力场 $\\sigma_{vM}$ (峰值 220.8 MPa)", fontsize=9.5, fontweight="bold", pad=4)
    axes[1].axhline(20.0, color="white", linestyle=":", lw=1.0, alpha=0.5)

    # (c) Case 1 Risk
    im3 = axes[2].pcolormesh(X_mesh, Y_mesh, risk_s1, cmap="YlOrRd", shading="gouraud", vmin=0, vmax=1.0)
    cb3 = fig.colorbar(im3, ax=axes[2], pad=0.015, aspect=12); cb3.set_label("危险度 $\\eta$", fontsize=8.5)
    axes[2].set_title("(c) Case 1 破坏危险度 $\\eta = \\sigma/400$ (峰值 0.36, 充分安全)", fontsize=9.5, fontweight="bold", pad=4)

    # (d) Case 3 Risk
    im4 = axes[3].pcolormesh(X_mesh, Y_mesh, risk_c3, cmap="YlOrRd", shading="gouraud", vmin=0, vmax=1.0)
    cb4 = fig.colorbar(im4, ax=axes[3], pad=0.015, aspect=12); cb4.set_label("危险度 $\\eta$", fontsize=8.5)
    axes[3].set_title("(d) Case 3 破坏危险度 $\\eta = \\sigma/245$ (峰值 0.90, 濒临脆断极限 1.0)", fontsize=9.5, fontweight="bold", pad=4)

    for ax in axes:
        ax.set_ylabel("y (mm)", fontsize=8.5)
        ax.set_ylim(0, 40)
        ax.set_xlim(0, 300)
        ax.set_aspect("equal", adjustable="box")

    axes[3].set_xlabel("x (mm) [0 = 固定端, 300 = 载荷端]", fontsize=9.5)
    plt.tight_layout()
    fig3_path = os.path.join(out_dir, "fig3_spatial_stress_contours.png")
    plt.savefig(fig3_path)
    plt.close()
    print(f"Saved: {fig3_path}")


if __name__ == "__main__":
    run_stress_distribution_analysis()
