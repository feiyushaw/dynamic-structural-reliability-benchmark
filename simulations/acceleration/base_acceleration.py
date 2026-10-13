"""
Case 3: High-Rate Impact Stamping & Deceleration Coupling
Simulates emergency panic foot impact (10 ms initial impact spike up to 2250 N decaying to 1500 N)
coupled with vehicle base deceleration (-1.0g).
Evaluates high strain-rate constitutive escalation (Cowper-Symonds) and dynamic embrittlement risk.
"""

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))

import json
import numpy as np
from dolfinx import fem
from dolfinx.fem import petsc
from petsc4py import PETSc
from simulations.common import (
    BeamModel,
    cowper_symonds_yield,
    dynamic_ductility_ratio,
)


def run_acceleration_case(output_file: str = "results/case3_acceleration.json"):
    print("=" * 65)
    print("[Case 3] High-Rate Impact Stamping (10 ms Impact Spike + 1g Deceleration)")
    model = BeamModel()

    dt = 0.0002
    t_span = np.arange(0.0, 0.20 + dt, dt)
    n_steps = len(t_span)
    g = 9.80665

    # Impact force profile: 0-10ms initial impact spike to 2250N (1.5x),
    # 10-25ms transition to 1500N sustained, 100-200ms unload
    F_t = np.zeros(n_steps, dtype=np.float64)
    # Vehicle longitudinal deceleration profile (-1g)
    a_x_t = np.zeros(n_steps, dtype=np.float64)

    for i, t in enumerate(t_span):
        if t < 0.010:  # Rapid initial impact
            F_t[i] = 2250.0 * (t / 0.010)
        elif t < 0.025:  # Decay to sustained pedal force
            ratio = (t - 0.010) / 0.015
            F_t[i] = 2250.0 * (1.0 - ratio) + 1500.0 * ratio
        elif t <= 0.10:
            F_t[i] = 1500.0
        else:
            F_t[i] = 1500.0 * max(0.0, 1.0 - (t - 0.10) / 0.10)

        # Vehicle deceleration ramps up as hydraulic brake builds up
        if t < 0.03:
            a_x_t[i] = 0.0
        elif t <= 0.06:
            a_x_t[i] = -1.0 * g * ((t - 0.03) / 0.03)
        elif t <= 0.12:
            a_x_t[i] = -1.0 * g
        else:
            a_x_t[i] = -1.0 * g * max(0.0, 1.0 - (t - 0.12) / 0.08)

    # Newmark setup
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
    history_edot, history_sig_y, history_duct = [], [], []
    uh = fem.Function(model.V)
    vh = fem.Function(model.V)

    for step in range(n_steps):
        t = t_span[step]
        F_val = F_t[step]
        ax_val = a_x_t[step]

        temp_m.zeroEntries()
        temp_m.axpy(a0, u_curr); temp_m.axpy(a2, v_curr); temp_m.axpy(a3, a_curr)

        temp_c.zeroEntries()
        temp_c.axpy(a1, u_curr); temp_c.axpy(a4, v_curr); temp_c.axpy(a5, a_curr)

        rhs.zeroEntries()
        if abs(F_val) > 1e-12:
            rhs.axpy(F_val, model.vec_tip_y)
        if abs(ax_val) > 1e-12:
            rhs.axpy(-model.rho * ax_val, model.vec_body_x)

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

            v_curr.copy(vh.x.petsc_vec)
            vh.x.scatter_forward()
            edot = model.compute_strain_rate(vh)
            sig_y = cowper_symonds_yield(edot) / 1e6
            duct = dynamic_ductility_ratio(edot)

            history_edot.append(edot)
            history_sig_y.append(sig_y)
            history_duct.append(duct)
        else:
            history_stress.append(history_stress[-1])
            history_edot.append(history_edot[-1])
            history_sig_y.append(history_sig_y[-1])
            history_duct.append(history_duct[-1])

    peak_disp = float(np.max(history_disp))
    peak_stress = float(np.max(history_stress))
    peak_edot = float(np.max(history_edot))
    peak_sig_y = float(np.max(history_sig_y))
    min_duct = float(np.min(history_duct)) * 100.0

    print(f"  Peak Tip Deflection     : {peak_disp:.3f} mm")
    print(f"  Peak von Mises Stress   : {peak_stress:.2f} MPa")
    print(f"  Peak Transient Strain Rate: {peak_edot:.2e} s^-1 (High Strain-Rate Regime)")
    print(f"  Dynamic Yield Strength  : {peak_sig_y:.2f} MPa (Rate Hardening: +{(peak_sig_y/250-1)*100:.1f}%)")
    print(f"  Min Ductility Retention : {min_duct:.1f}% (Severe Embrittlement Risk!)")
    print("=" * 65)

    result = {
        "case": "Case 3: High-Rate Impact Stamping (10ms Spike + 1g Accel)",
        "time_s": history_t,
        "force_N": F_t.tolist(),
        "accel_x_mps2": a_x_t.tolist(),
        "deflection_mm": history_disp,
        "stress_MPa": history_stress,
        "strain_rate_s1": history_edot,
        "dynamic_yield_MPa": history_sig_y,
        "ductility_ratio": history_duct,
        "peak_deflection_mm": peak_disp,
        "peak_stress_MPa": peak_stress,
        "peak_strain_rate_s1": peak_edot,
        "peak_dynamic_yield_MPa": peak_sig_y,
        "min_ductility_pct": min_duct,
    }

    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2, ensure_ascii=False)
    return result


if __name__ == "__main__":
    run_acceleration_case()
