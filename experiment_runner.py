"""
Experiment Runner for Brake Pedal Structural Dynamic Reliability Study.
Executes Case S1, S2, D1, D2, D3, Step Benchmark, and Parametric Sweep.
Saves all simulation curves and key indices to JSON for scientific plotting.
"""

import os
import json
import numpy as np
from solver import BrakePedalFEA


def run_all_experiments(output_dir: str = "results"):
    os.makedirs(output_dir, exist_ok=True)
    print("=" * 75)
    print("Initializing Brake Pedal Finite Element Model (DOLFINx 0.11)...")
    fea = BrakePedalFEA()
    print(f"Mesh: {fea.nx}x{fea.ny}x{fea.nz} Hex8 elements.")
    print(f"Component dimensions: L={fea.L*1e3:.0f}mm, H={fea.H*1e3:.0f}mm, B={fea.B*1e3:.0f}mm")
    print(f"Mass: {fea.mass:.3f} kg | Fundamental bending frequency: {fea.freq_1:.1f} Hz (Period: {1000/fea.freq_1:.2f} ms)")
    print("=" * 75)

    results = {}
    g = 9.80665

    # ---------------------------------------------------------
    # Case S1: 1500 N Standard Quasi-Static Test (QC/T 788-2018)
    # ---------------------------------------------------------
    print("\n[Running Case S1: 1500 N Static Standard Benchmark (QC/T 788-2018)]")
    s1 = fea.solve_static(F_load=1500.0)
    print(f"  -> Case S1 Max Deflection: {s1['max_deflection_mm']:.3f} mm")
    print(f"  -> Case S1 Max von Mises Stress: {s1['max_stress_MPa']:.2f} MPa")
    results["Case_S1"] = {
        "name": "Case S1: 1500N 静态标准测试 (国标基准)",
        "type": "static",
        "F_load": 1500.0,
        "max_deflection_mm": s1["max_deflection_mm"],
        "max_stress_MPa": s1["max_stress_MPa"],
    }

    # ---------------------------------------------------------
    # Case S2: 3000 N Strict OEM Quasi-Static Test
    # ---------------------------------------------------------
    print("\n[Running Case S2: 3000 N Static Strict Benchmark]")
    s2 = fea.solve_static(F_load=3000.0)
    print(f"  -> Case S2 Max Deflection: {s2['max_deflection_mm']:.3f} mm")
    print(f"  -> Case S2 Max von Mises Stress: {s2['max_stress_MPa']:.2f} MPa")
    results["Case_S2"] = {
        "name": "Case S2: 3000N 静态极限加严测试 (2倍冗余)",
        "type": "static",
        "F_load": 3000.0,
        "max_deflection_mm": s2["max_deflection_mm"],
        "max_stress_MPa": s2["max_stress_MPa"],
    }

    # Common time grid for 0.2s duration
    t_total = 0.2
    dt = 0.0002
    t_span = np.arange(0.0, t_total + dt, dt)
    n_steps = len(t_span)

    # ---------------------------------------------------------
    # Case D1: Rapid Driver Braking (0-50ms linear ramp to 1500N, 50-100ms hold, 100-200ms unload)
    # ---------------------------------------------------------
    print("\n[Running Case D1: Rapid Pedal Application (50ms Rise Time)]")
    F_d1 = np.zeros(n_steps, dtype=np.float64)
    for i, t in enumerate(t_span):
        if t <= 0.05:
            F_d1[i] = 1500.0 * (t / 0.05)
        elif t <= 0.10:
            F_d1[i] = 1500.0
        else:
            F_d1[i] = 1500.0 * max(0.0, 1.0 - (t - 0.10) / 0.10)

    d1 = fea.solve_transient(t_span, F_d1, record_stress_interval=5)
    daf_d1 = d1["peak_stress_MPa"] / s1["max_stress_MPa"]
    print(f"  -> Case D1 Peak Deflection: {d1['peak_deflection_mm']:.3f} mm")
    print(f"  -> Case D1 Peak von Mises Stress: {d1['peak_stress_MPa']:.2f} MPa (DAF = {daf_d1:.3f})")

    results["Case_D1"] = {
        "name": "Case D1: 快速平滑踩踏 (50ms上升沿)",
        "type": "dynamic",
        "time": d1["time"].tolist(),
        "force": F_d1.tolist(),
        "deflection_mm": d1["tip_deflection_mm"].tolist(),
        "stress_MPa": d1["max_stress_MPa"].tolist(),
        "peak_deflection_mm": d1["peak_deflection_mm"],
        "peak_stress_MPa": d1["peak_stress_MPa"],
        "daf": daf_d1,
    }

    # ---------------------------------------------------------
    # Case D1_Impact: Panic Stamping (Driver foot impact spike: 2250 N spike at 10ms decaying to 1500 N)
    # ---------------------------------------------------------
    print("\n[Running Case D1_Impact: Emergency Panic Stamping Impact]")
    F_d1_imp = np.zeros(n_steps, dtype=np.float64)
    for i, t in enumerate(t_span):
        if t < 0.010:  # 0-10ms initial impact spike up to 2250 N (1.5x nominal)
            F_d1_imp[i] = 2250.0 * (t / 0.010)
        elif t < 0.025:  # 10-25ms transition back to sustained 1500 N
            ratio = (t - 0.010) / 0.015
            F_d1_imp[i] = 2250.0 * (1.0 - ratio) + 1500.0 * ratio
        elif t <= 0.10:  # 25-100ms sustained 1500 N
            F_d1_imp[i] = 1500.0
        else:  # unload
            F_d1_imp[i] = 1500.0 * max(0.0, 1.0 - (t - 0.10) / 0.10)

    d1_imp = fea.solve_transient(t_span, F_d1_imp, record_stress_interval=5)
    daf_d1_imp = d1_imp["peak_stress_MPa"] / s1["max_stress_MPa"]
    print(f"  -> Case D1_Impact Peak Deflection: {d1_imp['peak_deflection_mm']:.3f} mm")
    print(f"  -> Case D1_Impact Peak von Mises Stress: {d1_imp['peak_stress_MPa']:.2f} MPa (DAF = {daf_d1_imp:.3f})")

    results["Case_D1_Impact"] = {
        "name": "Case D1-Imp: 紧急跺脚冲击踩踏 (2250N初峰值)",
        "type": "dynamic",
        "time": d1_imp["time"].tolist(),
        "force": F_d1_imp.tolist(),
        "deflection_mm": d1_imp["tip_deflection_mm"].tolist(),
        "stress_MPa": d1_imp["max_stress_MPa"].tolist(),
        "peak_deflection_mm": d1_imp["peak_deflection_mm"],
        "peak_stress_MPa": d1_imp["peak_stress_MPa"],
        "daf": daf_d1_imp,
    }

    # ---------------------------------------------------------
    # Case D2: Coupled Vehicle Braking Deceleration Inertia (-1g base acceleration)
    # ---------------------------------------------------------
    print("\n[Running Case D2: Driver Braking Coupled with -1g Vehicle Deceleration]")
    a_x_d2 = np.zeros(n_steps, dtype=np.float64)
    for i, t in enumerate(t_span):
        if t < 0.03:
            a_x_d2[i] = 0.0
        elif t <= 0.06:
            a_x_d2[i] = -1.0 * g * ((t - 0.03) / 0.03)
        elif t <= 0.12:
            a_x_d2[i] = -1.0 * g
        else:
            a_x_d2[i] = -1.0 * g * max(0.0, 1.0 - (t - 0.12) / 0.08)

    d2 = fea.solve_transient(t_span, F_d1, accel_x_t=a_x_d2, record_stress_interval=5)
    daf_d2 = d2["peak_stress_MPa"] / s1["max_stress_MPa"]
    print(f"  -> Case D2 Peak Deflection: {d2['peak_deflection_mm']:.3f} mm")
    print(f"  -> Case D2 Peak von Mises Stress: {d2['peak_stress_MPa']:.2f} MPa (DAF = {daf_d2:.3f})")

    results["Case_D2"] = {
        "name": "Case D2: 车辆制动惯性耦合 (-1g 基座减速)",
        "type": "dynamic",
        "time": d2["time"].tolist(),
        "force": F_d1.tolist(),
        "accel_x": a_x_d2.tolist(),
        "deflection_mm": d2["tip_deflection_mm"].tolist(),
        "stress_MPa": d2["max_stress_MPa"].tolist(),
        "peak_deflection_mm": d2["peak_deflection_mm"],
        "peak_stress_MPa": d2["peak_stress_MPa"],
        "daf": daf_d2,
    }

    # ---------------------------------------------------------
    # Case D3: Road Impact Shock (3g 30ms Half-Sine vertical shock)
    # ---------------------------------------------------------
    print("\n[Running Case D3: Road Impact Shock (3g 30ms Half-Sine)]")
    a_y_d3 = np.zeros(n_steps, dtype=np.float64)
    t_shock_start = 0.05
    T_shock = 0.030
    for i, t in enumerate(t_span):
        if t_shock_start <= t <= t_shock_start + T_shock:
            a_y_d3[i] = 3.0 * g * np.sin(np.pi * (t - t_shock_start) / T_shock)

    d3 = fea.solve_transient(t_span, F_d1, accel_x_t=a_x_d2, accel_y_t=a_y_d3, record_stress_interval=5)
    daf_d3 = d3["peak_stress_MPa"] / s1["max_stress_MPa"]
    print(f"  -> Case D3 Peak Deflection: {d3['peak_deflection_mm']:.3f} mm")
    print(f"  -> Case D3 Peak von Mises Stress: {d3['peak_stress_MPa']:.2f} MPa (DAF = {daf_d3:.3f})")

    results["Case_D3"] = {
        "name": "Case D3: 极端道路冲击 (3g 半正弦垂向冲击)",
        "type": "dynamic",
        "time": d3["time"].tolist(),
        "force": F_d1.tolist(),
        "accel_y": a_y_d3.tolist(),
        "deflection_mm": d3["tip_deflection_mm"].tolist(),
        "stress_MPa": d3["max_stress_MPa"].tolist(),
        "peak_deflection_mm": d3["peak_deflection_mm"],
        "peak_stress_MPa": d3["peak_stress_MPa"],
        "daf": daf_d3,
    }

    # ---------------------------------------------------------
    # Case Step: Instantaneous Step Load Benchmark (Theoretical Upper Bound)
    # ---------------------------------------------------------
    print("\n[Running Case Step: Instantaneous Step Load (1500 N Theoretical Limit)]")
    t_step_span = np.arange(0.0, 0.05 + dt, dt)
    F_step = np.full_like(t_step_span, 1500.0)
    d_step = fea.solve_transient(t_step_span, F_step, record_stress_interval=1)
    daf_step = d_step["peak_stress_MPa"] / s1["max_stress_MPa"]
    print(f"  -> Case Step Peak Deflection: {d_step['peak_deflection_mm']:.3f} mm")
    print(f"  -> Case Step Peak von Mises Stress: {d_step['peak_stress_MPa']:.2f} MPa (DAF = {daf_step:.3f})")

    results["Case_Step"] = {
        "name": "Case Step: 突加阶跃载荷 (理论上限基准)",
        "type": "dynamic",
        "time": d_step["time"].tolist(),
        "deflection_mm": d_step["tip_deflection_mm"].tolist(),
        "stress_MPa": d_step["max_stress_MPa"].tolist(),
        "peak_deflection_mm": d_step["peak_deflection_mm"],
        "peak_stress_MPa": d_step["peak_stress_MPa"],
        "daf": daf_step,
    }

    # ---------------------------------------------------------
    # Parametric Study: Loading Rate (Rise Time) Sensitivity Sweep
    # ---------------------------------------------------------
    print("\n[Running Parametric Sweep: Rise Time Sensitivity Analysis]")
    rise_times = [1.0, 0.20, 0.10, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001]
    sweep_results = []

    for t_rise in rise_times:
        t_sim_end = max(0.05, t_rise + 0.03)
        t_vec = np.arange(0.0, t_sim_end + dt, dt)
        f_vec = np.array([1500.0 * min(t / t_rise, 1.0) for t in t_vec])

        res_sw = fea.solve_transient(t_vec, f_vec, record_stress_interval=2)
        daf_sw = res_sw["peak_stress_MPa"] / s1["max_stress_MPa"]
        print(f"  Rise time = {t_rise*1000:6.1f} ms | Peak Stress = {res_sw['peak_stress_MPa']:.2f} MPa | DAF = {daf_sw:.3f}")
        sweep_results.append({
            "rise_time_ms": t_rise * 1000.0,
            "peak_stress_MPa": res_sw["peak_stress_MPa"],
            "daf": daf_sw,
        })

    results["Parametric_Sweep"] = sweep_results

    # ---------------------------------------------------------
    # Theoretical Decoupling Analysis: 2-ton vs 3-ton Vehicle
    # ---------------------------------------------------------
    print("\n[Computing Vehicle Weight & Inertial Decoupling Analysis]")
    m_pedal = fea.mass
    f_ground_2t = 2000.0 * g
    f_ground_3t = 3000.0 * g
    f_pedal_inertia = m_pedal * (1.0 * g)

    decoupling_analysis = {
        "component_mass_kg": m_pedal,
        "driver_foot_force_N": 1500.0,
        "vehicle_2t": {
            "mass_kg": 2000.0,
            "ground_adhesion_kN": f_ground_2t / 1000.0,
            "deceleration_g": 1.0,
            "pedal_inertia_N": f_pedal_inertia,
            "pedal_inertia_ratio_pct": (f_pedal_inertia / 1500.0) * 100.0,
        },
        "vehicle_3t": {
            "mass_kg": 3000.0,
            "ground_adhesion_kN": f_ground_3t / 1000.0,
            "deceleration_g": 1.0,
            "pedal_inertia_N": f_pedal_inertia,
            "pedal_inertia_ratio_pct": (f_pedal_inertia / 1500.0) * 100.0,
        },
        "inertia_difference_N": 0.0,
        "key_takeaway": (
            "地面制动力通过轮胎-悬架-副车架-白车身闭环平衡整车动能，不会传递给驾驶舱踏板。"
            "踏板自身仅承受微弱的自身惯性力 (~9.2 N)，仅占踩踏力的 0.61%，2吨与3吨车在此项上毫无本质差异。"
        )
    }
    results["Decoupling_Analysis"] = decoupling_analysis

    # Save to JSON
    json_path = os.path.join(output_dir, "simulation_results.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(results, f, ensure_ascii=False, indent=2)
    print(f"\nAll simulation results successfully exported to: {json_path}")
    print("=" * 75)


if __name__ == "__main__":
    run_all_experiments()
