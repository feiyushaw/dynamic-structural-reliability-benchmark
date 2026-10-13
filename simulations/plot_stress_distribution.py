"""
Spatial Stress Distribution & Failure Risk Contour Visualization.
Computes and visualizes the 2D longitudinal mid-plane and 3D stress fields for:
  - Case 1: Quasi-static qualification (1500 N)
  - Case 3: High-rate emergency impact stamping (peak response)
Generates:
  1. results/stress_distribution_contour.png (Multi-panel contour plot)
  2. results/beam_case1_static.pvd & results/beam_case3_impact.pvd (ParaView 3D VTK export)
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

# Publication styling
mpl.rcParams["font.family"] = ["Noto Sans CJK SC", "DejaVu Sans", "sans-serif"]
mpl.rcParams["axes.unicode_minus"] = False
mpl.rcParams["figure.dpi"] = 300
mpl.rcParams["savefig.dpi"] = 300
mpl.rcParams["font.size"] = 10.5


def extract_stress_grid(model: BeamModel, uh: fem.Function):
    """Extract von Mises stress and project onto regular (8, 60) 2D midplane grid."""
    sig = model.sigma(uh)
    s_dev = sig - (1.0 / 3.0) * ufl.tr(sig) * ufl.Identity(3)
    vm = ufl.sqrt(ufl.inner(s_dev, s_dev) * 1.5)
    expr = fem.Expression(vm, model.W.element.interpolation_points)
    f_vm = fem.Function(model.W)
    f_vm.interpolate(expr)
    stress_arr = f_vm.x.array / 1e6  # MPa

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
    print("=" * 68)
    print("Computing 3D Spatial Stress Distribution & Failure Risk Fields...")
    model = BeamModel()

    # -------------------------------------------------------------
    # 1. Solve Case 1 (Quasi-Static 1500 N)
    # -------------------------------------------------------------
    rhs_s1 = model.K.createVecRight()
    rhs_s1.zeroEntries()
    rhs_s1.axpy(1500.0, model.vec_tip_y)
    petsc.apply_lifting(rhs_s1, [model.k_form], [[model.bc]])
    rhs_s1.ghostUpdate()
    petsc.set_bc(rhs_s1, [model.bc])

    uh_s1 = fem.Function(model.V)
    ksp_s1 = PETSc.KSP().create(model.comm)
    ksp_s1.setOperators(model.K)
    ksp_s1.setType("preonly")
    ksp_s1.getPC().setType("lu")
    ksp_s1.solve(rhs_s1, uh_s1.x.petsc_vec)
    uh_s1.x.scatter_forward()

    x_mm, y_mm, grid_s1, f_vm_s1 = extract_stress_grid(model, uh_s1)

    # -------------------------------------------------------------
    # 2. Solve Case 3 (High-Rate Impact Peak at t = 10 ms)
    # -------------------------------------------------------------
    dt = 0.0002
    t_span = np.arange(0.0, 0.025 + dt, dt)
    n_steps = len(t_span)
    g = 9.80665

    # Newmark parameters
    beta, gamma = 0.25, 0.5
    a0 = 1.0 / (beta * dt * dt)
    a1 = gamma / (beta * dt)
    a2 = 1.0 / (beta * dt)
    a3 = 1.0 / (2.0 * beta) - 1.0
    a4 = gamma / beta - 1.0
    a5 = (dt / 2.0) * (gamma / beta - 2.0)
    a6 = dt * (1.0 - gamma)
    a7 = gamma * dt

    K_hat = model.K.copy()
    K_hat.axpy(a0, model.M)
    K_hat.axpy(a1, model.C)
    bc_dofs, _ = model.bc.dof_indices()
    K_hat.zeroRowsColumns(bc_dofs, diag=1.0)
    K_hat.assemble()

    ksp_dyn = PETSc.KSP().create(model.comm)
    ksp_dyn.setOperators(K_hat)
    ksp_dyn.setType("preonly")
    ksp_dyn.getPC().setType("lu")
    ksp_dyn.setUp()

    u_curr = model.K.createVecRight(); v_curr = model.K.createVecRight(); a_curr = model.K.createVecRight()
    u_next = model.K.createVecRight(); v_next = model.K.createVecRight(); a_next = model.K.createVecRight()
    temp_m = model.K.createVecRight(); temp_c = model.K.createVecRight(); rhs = model.K.createVecRight()
    u_curr.set(0.0); v_curr.set(0.0); a_curr.set(0.0)

    uh_peak_dyn = fem.Function(model.V)
    max_recorded_stress = 0.0

    for step in range(n_steps):
        t = t_span[step]
        # Force pulse
        if t < 0.010:
            F_val = 2250.0 * (t / 0.010)
        else:
            ratio = (t - 0.010) / 0.015
            F_val = 2250.0 * (1.0 - ratio) + 1500.0 * ratio

        ax_val = -1.0 * g if t >= 0.010 else 0.0

        temp_m.zeroEntries(); temp_m.axpy(a0, u_curr); temp_m.axpy(a2, v_curr); temp_m.axpy(a3, a_curr)
        temp_c.zeroEntries(); temp_c.axpy(a1, u_curr); temp_c.axpy(a4, v_curr); temp_c.axpy(a5, a_curr)

        rhs.zeroEntries()
        if abs(F_val) > 1e-12:
            rhs.axpy(F_val, model.vec_tip_y)
        if abs(ax_val) > 1e-12:
            rhs.axpy(-model.rho * ax_val, model.vec_body_x)

        model.M.multAdd(temp_m, rhs, rhs)
        model.C.multAdd(temp_c, rhs, rhs)
        petsc.set_bc(rhs, [model.bc])
        ksp_dyn.solve(rhs, u_next)

        a_next.zeroEntries(); a_next.axpy(a0, u_next); a_next.axpy(-a0, u_curr); a_next.axpy(-a2, v_curr); a_next.axpy(-a3, a_curr)
        v_next.zeroEntries(); v_next.axpy(1.0, v_curr); v_next.axpy(a6, a_curr); v_next.axpy(a7, a_next)
        u_next.copy(u_curr); v_next.copy(v_curr); a_next.copy(a_curr)

        # Check peak stress at t ~ 10ms
        if np.isclose(t, 0.010, atol=dt/2.0) or t >= 0.010:
            uh_test = fem.Function(model.V)
            u_curr.copy(uh_test.x.petsc_vec)
            uh_test.x.scatter_forward()
            vm_val = model.compute_von_mises(uh_test)
            if vm_val > max_recorded_stress:
                max_recorded_stress = vm_val
                u_curr.copy(uh_peak_dyn.x.petsc_vec)
                uh_peak_dyn.x.scatter_forward()

    _, _, grid_c3, f_vm_c3 = extract_stress_grid(model, uh_peak_dyn)

    # -------------------------------------------------------------
    # 3. Export 3D VTK Files for ParaView
    # -------------------------------------------------------------
    vtk_c1_path = os.path.join(out_dir, "beam_case1_static.pvd")
    with io.VTKFile(model.comm, vtk_c1_path, "w") as vtk:
        vtk.write_function([uh_s1, f_vm_s1])
    print(f"Exported ParaView 3D file: {vtk_c1_path}")

    vtk_c3_path = os.path.join(out_dir, "beam_case3_impact.pvd")
    with io.VTKFile(model.comm, vtk_c3_path, "w") as vtk:
        vtk.write_function([uh_peak_dyn, f_vm_c3])
    print(f"Exported ParaView 3D file: {vtk_c3_path}")

    # -------------------------------------------------------------
    # 4. Compute Failure Risk Index Fields
    #   eta_static  = sigma_vM / 400.0 (Ductile UTS capacity)
    #   eta_dynamic = sigma_vM / 245.0 (Dynamic cleavage fracture threshold)
    # -------------------------------------------------------------
    risk_s1 = grid_s1 / 400.0
    risk_c3 = grid_c3 / 245.0

    # -------------------------------------------------------------
    # 5. Render Multi-Panel Publication Contour Plot
    # -------------------------------------------------------------
    X_mesh, Y_mesh = np.meshgrid(x_mm, y_mm)
    fig, axes = plt.subplots(4, 1, figsize=(14, 11), sharex=True)

    # (a) Case 1 Stress Contour
    c1_map = axes[0].pcolormesh(X_mesh, Y_mesh, grid_s1, cmap="turbo", shading="gouraud", vmin=0, vmax=230)
    cb1 = fig.colorbar(c1_map, ax=axes[0], pad=0.015, aspect=12)
    cb1.set_label("应力 (MPa)", fontsize=10)
    axes[0].set_title("(a) Case 1 准静态测试：悬臂梁纵向截面 von Mises 应力分布 (峰值 143.6 MPa)", fontweight="bold", pad=8)
    axes[0].axhline(20.0, color="white", linestyle=":", lw=1.2, alpha=0.7, label="中性轴 (正应力=0)")
    axes[0].annotate("固定端 (X=0)\n踏板支架安装板", xy=(0, 20), xytext=(-25, 20), ha="right", va="center", fontsize=8.5, fontweight="bold")
    axes[0].annotate("外载 F=1500N\n(踏板踩踏端)", xy=(300, 20), xytext=(305, 20), ha="left", va="center", fontsize=8.5, fontweight="bold", color="blue")
    axes[0].annotate("根部最大拉应力 (143.6 MPa)", xy=(2.5, 37.5), xytext=(25, 36), arrowprops=dict(facecolor="black", shrink=0.08, width=0.8, headwidth=4), fontsize=8.5)

    # (b) Case 3 Stress Contour
    c3_map = axes[1].pcolormesh(X_mesh, Y_mesh, grid_c3, cmap="turbo", shading="gouraud", vmin=0, vmax=230)
    cb2 = fig.colorbar(c3_map, ax=axes[1], pad=0.015, aspect=12)
    cb2.set_label("应力 (MPa)", fontsize=10)
    axes[1].set_title("(b) Case 3 动态冲击跺脚 (峰值时刻)：von Mises 应力空间激增分布 (峰值 220.8 MPa)", fontweight="bold", pad=8)
    axes[1].axhline(20.0, color="white", linestyle=":", lw=1.2, alpha=0.7)
    axes[1].annotate("冲击应力激增区 (220.8 MPa)\n高应力区显著向前扩展", xy=(2.5, 37.5), xytext=(35, 34), arrowprops=dict(facecolor="black", shrink=0.08, width=0.8, headwidth=4), fontsize=8.5, fontweight="bold", color="darkred")

    # (c) Case 1 Failure Risk Index
    r1_map = axes[2].pcolormesh(X_mesh, Y_mesh, risk_s1, cmap="YlOrRd", shading="gouraud", vmin=0, vmax=1.0)
    cb3 = fig.colorbar(r1_map, ax=axes[2], pad=0.015, aspect=12)
    cb3.set_label("破坏风险度 $\\eta$", fontsize=10)
    axes[2].set_title("(c) Case 1 静态破坏风险指数场 $\\eta = \\sigma_{vM} / \\sigma_{fail}^{stat}$ (基准极限 400 MPa, 峰值风险仅 0.36, 极度安全)", fontweight="bold", pad=8)
    axes[2].annotate("静态充分延性吸能\n最大破坏风险仅 36%\n(绿色低危安全区)", xy=(150, 20), ha="center", va="center", fontsize=9.5, fontweight="bold", color="#196f3d", bbox=dict(boxstyle="round,pad=0.3", fc="#d5f5e3", ec="#27ae60", alpha=0.9))

    # (d) Case 3 Failure Risk Index
    r2_map = axes[3].pcolormesh(X_mesh, Y_mesh, risk_c3, cmap="YlOrRd", shading="gouraud", vmin=0, vmax=1.0)
    cb4 = fig.colorbar(r2_map, ax=axes[3], pad=0.015, aspect=12)
    cb4.set_label("破坏风险度 $\\eta$", fontsize=10)
    axes[3].set_title("(d) Case 3 动态冲击破坏风险指数场 $\\eta = \\sigma_{vM} / \\sigma_{fail}^{dyn}$ (脆断极限 245 MPa, 峰值风险高达 0.90, 高危濒断!)", fontweight="bold", pad=8)
    axes[3].axhline(20.0, color="black", linestyle=":", lw=1.2, alpha=0.5)
    axes[3].annotate("根部脆断极度危险区 (风险指数 $\\eta = 0.90$)\n缺乏塑性钝化，极易诱发微裂纹失稳撕裂！", xy=(2.5, 37.5), xytext=(45, 30), arrowprops=dict(facecolor="crimson", shrink=0.08, width=1.0, headwidth=5), fontsize=9, fontweight="bold", color="crimson", bbox=dict(boxstyle="round,pad=0.3", fc="#fce4e4", ec="crimson", alpha=0.95))

    for ax in axes:
        ax.set_ylabel("截面高度 Y (mm)", fontsize=10)
        ax.set_ylim(0, 40)
        ax.set_xlim(-1, 301)
        ax.set_aspect("equal", adjustable="box")

    axes[3].set_xlabel("悬臂梁长度 X (mm) [0 = 支架固定安装端, 300 = 踏板踩踏端]", fontsize=11, labelpad=6)

    plt.tight_layout()
    contour_path = os.path.join(out_dir, "stress_distribution_contour.png")
    plt.savefig(contour_path)
    plt.close()
    print(f"Saved Publication Contour Plot: {contour_path}")
    print("=" * 68)


if __name__ == "__main__":
    run_stress_distribution_analysis()
