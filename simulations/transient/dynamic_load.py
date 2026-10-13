"""
Case 2: Transient Dynamic Force Loading
Simulates rapid driver pedal application (0 -> 1500 N in 50 ms, hold 50 ms, release).
Solves: M * u'' + C * u' + K * u = F(t) using Newmark-beta method.
"""

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))

import json
import numpy as np
from dolfinx import fem
from dolfinx.fem import petsc
from petsc4py import PETSc
from simulations.common import BeamModel


def run_dynamic_case(output_file: str = "results/case2_dynamic.json"):
    print("=" * 60)
    print("[Case 2] Transient Dynamic Loading (50 ms Ramp)")
    model = BeamModel()

    # Time discretization
    dt = 0.0002
    t_span = np.arange(0.0, 0.20 + dt, dt)
    n_steps = len(t_span)

    # Force profile: 0-50ms linear ramp to 1500N, 50-100ms hold, 100-200ms unload
    F_t = np.zeros(n_steps, dtype=np.float64)
    for i, t in enumerate(t_span):
        if t <= 0.05:
            F_t[i] = 1500.0 * (t / 0.05)
        elif t <= 0.10:
            F_t[i] = 1500.0
        else:
            F_t[i] = 1500.0 * max(0.0, 1.0 - (t - 0.10) / 0.10)

    # Newmark parameters (beta=0.25, gamma=0.5)
    beta, gamma = 0.25, 0.5
    a0 = 1.0 / (beta * dt * dt)
    a1 = gamma / (beta * dt)
    a2 = 1.0 / (beta * dt)
    a3 = 1.0 / (2.0 * beta) - 1.0
    a4 = gamma / beta - 1.0
    a5 = (dt / 2.0) * (gamma / beta - 2.0)
    a6 = dt * (1.0 - gamma)
    a7 = gamma * dt

    # Dynamic effective stiffness: K_hat = K + a0*M + a1*C
    K_hat = model.K.copy()
    K_hat.axpy(a0, model.M)
    K_hat.axpy(a1, model.C)
    bc_dofs, _ = model.bc.dof_indices()
    K_hat.zeroRowsColumns(bc_dofs, diag=1.0)
    K_hat.assemble()

    ksp = PETSc.KSP().create(model.comm)
    ksp.setOperators(K_hat)
    ksp.setType("preonly")
    ksp.getPC().setType("lu")
    ksp.setUp()

    u_curr = model.K.createVecRight()
    v_curr = model.K.createVecRight()
    a_curr = model.K.createVecRight()
    u_next = model.K.createVecRight()
    v_next = model.K.createVecRight()
    a_next = model.K.createVecRight()
    temp_m = model.K.createVecRight()
    temp_c = model.K.createVecRight()
    rhs = model.K.createVecRight()

    u_curr.set(0.0); v_curr.set(0.0); a_curr.set(0.0)

    history_t, history_disp, history_stress = [], [], []
    uh = fem.Function(model.V)

    for step in range(n_steps):
        t = t_span[step]
        F_val = F_t[step]

        temp_m.zeroEntries()
        temp_m.axpy(a0, u_curr); temp_m.axpy(a2, v_curr); temp_m.axpy(a3, a_curr)

        temp_c.zeroEntries()
        temp_c.axpy(a1, u_curr); temp_c.axpy(a4, v_curr); temp_c.axpy(a5, a_curr)

        rhs.zeroEntries()
        if abs(F_val) > 1e-12:
            rhs.axpy(F_val, model.vec_tip_y)
        model.M.multAdd(temp_m, rhs, rhs)
        model.C.multAdd(temp_c, rhs, rhs)

        petsc.set_bc(rhs, [model.bc])
        ksp.solve(rhs, u_next)

        a_next.zeroEntries()
        a_next.axpy(a0, u_next); a_next.axpy(-a0, u_curr); a_next.axpy(-a2, v_curr); a_next.axpy(-a3, a_curr)

        v_next.zeroEntries()
        v_next.axpy(1.0, v_curr); v_next.axpy(a6, a_curr); v_next.axpy(a7, a_next)

        u_next.copy(u_curr); v_next.copy(v_curr); a_next.copy(a_curr)

        tip_disp = abs(float(np.min(u_curr.getArray()[1::3]))) * 1e3
        history_t.append(t)
        history_disp.append(tip_disp)

        if step % 5 == 0 or step == n_steps - 1:
            u_curr.copy(uh.x.petsc_vec)
            uh.x.scatter_forward()
            vm = model.compute_von_mises(uh) / 1e6
            history_stress.append(vm)
        else:
            history_stress.append(history_stress[-1])

    peak_disp = float(np.max(history_disp))
    peak_stress = float(np.max(history_stress))
    print(f"  Peak Tip Deflection : {peak_disp:.3f} mm")
    print(f"  Peak von Mises Stress: {peak_stress:.2f} MPa")
    print("=" * 60)

    result = {
        "case": "Case 2: Dynamic Force (50ms Ramp)",
        "time_s": history_t,
        "force_N": F_t.tolist(),
        "deflection_mm": history_disp,
        "stress_MPa": history_stress,
        "peak_deflection_mm": peak_disp,
        "peak_stress_MPa": peak_stress,
    }

    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2, ensure_ascii=False)
    return result


if __name__ == "__main__":
    run_dynamic_case()
