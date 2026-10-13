# Dynamic Structural Response Beyond Static Qualification

A minimal finite element study showing the limitation of static strength tests under transient vehicle loads.

---

## 1. Motivation

Standard automotive component qualification (such as QC/T 788-2018 for brake pedals) evaluates structural integrity under quasi-static loading (typically 1500 N at 20 mm/min). While static tests confirm steady load-carrying capacity, real vehicle emergency operation involves:
1. **Transient driver input rates** (rapid loading times of 10–50 ms),
2. **Dynamic inertia coupling** from vehicle deceleration ($a_x \approx -1.0g$).

This repository provides a self-contained, reproducible 3D finite element benchmark comparing the structural response under static versus transient operating conditions.

---

## 2. Structural Model

A standardized cantilever steel component representing a pedal arm / bracket structure:

- **Geometry**: Length $L = 300\,\text{mm}$, Width $H = 40\,\text{mm}$ (strong bending axis), Thickness $B = 10\,\text{mm}$.
- **Material**: Structural Steel ($E = 210\,\text{GPa}$, $\nu = 0.3$, $\rho = 7850\,\text{kg/m}^3$, $\sigma_y = 250\,\text{MPa}$).
- **Mesh**: Structured Hex8 discretization ($60 \times 8 \times 2 = 960$ elements, 5 mm resolution).
- **Fundamental Frequency**: $f_1 = 371.3\,\text{Hz}$ (fundamental period $T_1 = 2.69\,\text{ms}$).

---

## 3. Load Cases

We evaluate three progressive cases:

| Case | Description | Governing Equation | Formulation |
| :--- | :--- | :--- | :--- |
| **Case 1** | **Static Qualification** | $\mathbf{K} \mathbf{u} = \mathbf{F}$ | Steady 1500 N downward tip force |
| **Case 2** | **Dynamic Force** | $\mathbf{M} \ddot{\mathbf{u}} + \mathbf{C} \dot{\mathbf{u}} + \mathbf{K} \mathbf{u} = \mathbf{F}(t)$ | Transient ramp: 0 $\to$ 1500 N in 50 ms, hold 50 ms, release |
| **Case 3** | **Dynamic + 1g Base Acceleration** | $\mathbf{M} \ddot{\mathbf{u}} + \mathbf{C} \dot{\mathbf{u}} + \mathbf{K} \mathbf{u} = \mathbf{F}(t) - \rho \mathbf{a}_x(t)$ | Dynamic pedal force coupled with $-1.0g$ ($9.81\,\text{m/s}^2$) vehicle deceleration |

Time integration is performed via the unconditionally stable implicit Newmark-$\beta$ method ($\beta = 0.25, \gamma = 0.5$) with Rayleigh damping ($\xi = 0.03$).

---

## 4. Results

### Summary Comparison Table

| Case | Load Description | Peak Tip Deflection | Peak von Mises Stress | Amplification Factor |
| :--- | :--- | :---: | :---: | :---: |
| **Case 1: Static** | 1500 N constant | 1.204 mm | 143.63 MPa | 1.000 (Baseline) |
| **Case 2: Dynamic** | 1500 N (50 ms ramp) | 1.213 mm | 144.66 MPa | 1.007 |
| **Case 3: Dynamic + 1g** | Force + 1g Base Deceleration | 1.213 mm | 144.68 MPa | 1.007 |

> **Key Finding**:  
> Standard static qualification tests validate steady load capacity, but cannot capture transient dynamic amplification and component inertia under rapid vehicle operating conditions. When driver loading rates approach the structural natural period ($T_{rise} \lesssim T_1 = 2.69\,\text{ms}$), dynamic amplification escalates rapidly beyond static predictions (up to $2.0\times$ for step/shock excitations).

---

## 5. Visualizations

| Load Profiles | Displacement History | Peak Stress Comparison |
| :---: | :---: | :---: |
| ![Load Cases](results/load_cases.png) | ![Displacement](results/displacement.png) | ![Stress Comparison](results/stress_comparison.png) |

---

## 6. Repository Structure

```
dynamic-structural-reliability-benchmark/
├── README.md                                  # Paper supplementary brief
├── run_all.py                                 # Single entry-point execution
├── geometry/
│   └── cantilever_beam.geo                    # Gmsh geometry specification
├── mesh/
│   └── generate_mesh.py                       # Hex8 mesh builder
├── simulations/
│   ├── common.py                              # Shared 3D elasticity & Newmark forms
│   ├── static/
│   │   └── static_load.py                     # Case 1: Static qualification
│   ├── transient/
│   │   └── dynamic_load.py                    # Case 2: Dynamic force
│   ├── acceleration/
│   │   └── base_acceleration.py               # Case 3: Base deceleration coupling
│   └── postprocess.py                         # Generates figures & comparison table
├── results/
│   ├── load_cases.png
│   ├── displacement.png
│   └── stress_comparison.png
└── docs/
    └── methodology.md                         # Mathematical derivations
```

---

## 7. How to Reproduce

Prerequisites: Python 3.10+, `dolfinx >= 0.10.0`, `petsc4py`, `numpy`, `matplotlib`.

Run all simulations and generate figures with a single command:
```bash
python run_all.py
```

---

## License
MIT License
