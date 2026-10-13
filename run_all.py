"""
Top-level execution script for Cantilever Beam FEM Validation Study.
Sequentially runs:
  1. simulations/static/static_load.py (Case 1)
  2. simulations/transient/dynamic_load.py (Case 2)
  3. simulations/acceleration/base_acceleration.py (Case 3)
  4. simulations/postprocess.py (generates rate-dependent constitutive curves & summary table)
  5. simulations/plot_stress_distribution.py (generates 2D/3D spatial stress distribution contour & VTK files)
"""

import sys
import os

from simulations.static.static_load import run_static_case
from simulations.transient.dynamic_load import run_dynamic_case
from simulations.acceleration.base_acceleration import run_acceleration_case
from simulations.postprocess import postprocess_all
from simulations.plot_stress_distribution import run_stress_distribution_analysis


def main():
    print("=" * 72)
    print("  Minimal FEM Validation: Static vs Dynamic Structural Reliability")
    print("=" * 72)

    # 1. Static qualification
    run_static_case()

    # 2. Dynamic force
    run_dynamic_case()

    # 3. Dynamic force + base deceleration
    run_acceleration_case()

    # 4. Post-processing & constitutive figure generation
    postprocess_all()

    # 5. Spatial stress distribution contour & VTK generation
    run_stress_distribution_analysis()


if __name__ == "__main__":
    main()
