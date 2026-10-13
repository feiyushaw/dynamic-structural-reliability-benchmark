"""
Top-level execution script for Cantilever Beam FEM Validation Study.
Sequentially runs:
  1. simulations/static/static_load.py (Case 1: 1500 N Static)
  2. simulations/transient/dynamic_load.py (Case 2: 50 ms Dynamic Ramp)
  3. simulations/acceleration/base_acceleration.py (Case 3: 10 ms Impact Stamping)
  4. simulations/plot_schematic.py (Figure 1: 2D Model & Load Waveforms)
  5. simulations/postprocess.py (Figure 2: Three Cases Systematic Comparison Suite)
  6. simulations/plot_stress_distribution.py (Figure 3: 2D Stress Contours & 3D VTK)
"""

import sys
import os

from simulations.static.static_load import run_static_case
from simulations.transient.dynamic_load import run_dynamic_case
from simulations.acceleration.base_acceleration import run_acceleration_case
from simulations.plot_schematic import plot_model_and_loading
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

    # 4. Figure 1: 2D Model & Load Waveforms
    plot_model_and_loading()

    # 5. Figure 2: Three Cases Systematic Comparison Suite
    postprocess_all()

    # 6. Figure 3: 2D Stress Contours & 3D VTK
    run_stress_distribution_analysis()


if __name__ == "__main__":
    main()
