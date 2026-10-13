"""
Case 1: Static Qualification Test (1500 N Benchmark)
Solves linear elasticity problem: K * u = F
Evaluates quasi-static strain rate (~10^-4 s^-1) and baseline yield stress.
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
    STATIC_STRAIN_RATE,
)


def run_static_case(F_load: float = 1500.0, output_file: str = "results/case1_static.json"):
    print("=" * 65)
    print(f"[Case 1] Linear Static Qualification: F = {F_load:.1f} N")
    model = BeamModel()

    rhs = model.K.createVecRight()
    rhs.zeroEntries()
    rhs.axpy(F_load, model.vec_tip_y)

    petsc.apply_lifting(rhs, [model.k_form], [[model.bc]])
    rhs.ghostUpdate()
    petsc.set_bc(rhs, [model.bc])

    uh = fem.Function(model.V)
    ksp = PETSc.KSP().create(model.comm)
    ksp.setOperators(model.K)
    ksp.setType("preonly")
    ksp.getPC().setType("lu")
    ksp.solve(rhs, uh.x.petsc_vec)
    uh.x.scatter_forward()

    max_deflection_mm = abs(float(np.min(uh.x.array[1::3]))) * 1e3
    max_stress_MPa = model.compute_von_mises(uh) / 1e6

    # Quasi-static strain rate at 20 mm/min loading velocity (QC/T 788 standard)
    # delta_t ~ 3.6 s to reach 1.2 mm tip deflection
    edot_static = STATIC_STRAIN_RATE
    sigma_y_MPa = cowper_symonds_yield(edot_static) / 1e6
    ductility_pct = dynamic_ductility_ratio(edot_static) * 100.0

    print(f"  Max Tip Deflection     : {max_deflection_mm:.3f} mm")
    print(f"  Max von Mises Stress   : {max_stress_MPa:.2f} MPa")
    print(f"  Quasi-Static Strain Rate: {edot_static:.1e} s^-1")
    print(f"  Baseline Yield Strength: {sigma_y_MPa:.2f} MPa")
    print(f"  Ductility Retention    : {ductility_pct:.1f}% (Fully ductile regime)")
    print("=" * 65)

    result = {
        "case": "Case 1: Static Qualification (QC/T 788 Benchmark)",
        "load_N": F_load,
        "max_deflection_mm": max_deflection_mm,
        "max_stress_MPa": max_stress_MPa,
        "strain_rate_s1": edot_static,
        "dynamic_yield_MPa": sigma_y_MPa,
        "ductility_ratio": ductility_pct / 100.0,
    }

    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2, ensure_ascii=False)
    return result


if __name__ == "__main__":
    run_static_case()
