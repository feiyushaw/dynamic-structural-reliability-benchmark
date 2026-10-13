"""
Case 1: Static Qualification Test (1500 N)
Solves linear elasticity problem: K * u = F
Outputs maximum deflection and peak von Mises stress.
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


def run_static_case(F_load: float = 1500.0, output_file: str = "results/case1_static.json"):
    print("=" * 60)
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

    print(f"  Max Tip Deflection : {max_deflection_mm:.3f} mm")
    print(f"  Max von Mises Stress: {max_stress_MPa:.2f} MPa")
    print("=" * 60)

    result = {
        "case": "Case 1: Static Qualification",
        "load_N": F_load,
        "max_deflection_mm": max_deflection_mm,
        "max_stress_MPa": max_stress_MPa,
    }

    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2, ensure_ascii=False)
    return result


if __name__ == "__main__":
    run_static_case()
