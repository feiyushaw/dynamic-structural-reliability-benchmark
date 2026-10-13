# Dynamic Structural Response Beyond Static Qualification: Rate-Dependent Impact Dynamics

A minimal 3D finite element benchmark demonstrating how dynamic loading impact rates alter material constitutive behavior (Cowper-Symonds strain-rate hardening & dynamic embrittlement), exposing the fundamental limitations of quasi-static strength qualification.

---

## 1. Motivation: Rate-Dependent Constitutive Shift

Standard automotive component qualification (e.g., QC/T 788-2018 for brake pedals) specifies quasi-static testing at $20\,\text{mm/min}$ ($\dot{\varepsilon} \sim 10^{-4}\,\text{s}^{-1}$). Under slow loading, structural steel exhibits high plastic ductility ($\varepsilon_f \approx 25\% - 30\%$), allowing components to redistribute local stresses without cracking.

In real-world emergency panic braking, driver foot impact occurs over $10 - 50\,\text{ms}$ ($\dot{\varepsilon} \sim 10^{-1} - 10^{2}\,\text{s}^{-1}$). This **multi-order-of-magnitude surge in strain rate** fundamentally alters the material constitutive law:
1. **Strain-Rate Hardening**: Flow stress and yield threshold escalate via dislocation drag (modeled by the Cowper-Symonds law);
2. **Dynamic Embrittlement**: Plastic strain capacity and fracture toughness $K_{IC}$ drop significantly ($< 63\%$ ductility retention), transforming ductile steel into a notch-sensitive brittle regime prone to abrupt fracture.

---

## 2. Structural Model & Constitutive Formulation

- **Geometry**: Cantilever beam ($L = 300\,\text{mm}$, strong-axis height $H = 40\,\text{mm}$, thickness $B = 10\,\text{mm}$).
- **Material Baseline**: Structural steel ($E = 210\,\text{GPa}, \nu = 0.3, \rho = 7850\,\text{kg/m}^3, \sigma_0 = 250\,\text{MPa}$).
- **Rate-Dependent Law (Cowper-Symonds)**:
  $$\sigma_y(\dot{\varepsilon}_{eq}) = \sigma_0 \left[ 1 + \left( \frac{\max(0, \dot{\varepsilon}_{eq} - \dot{\varepsilon}_0)}{D} \right)^{1/p} \right] \quad (D = 40.4\,\text{s}^{-1}, p = 5.0, \dot{\varepsilon}_0 = 10^{-4}\,\text{s}^{-1})$$
- **Mesh**: 3D Hex8 structured mesh ($60 \times 8 \times 2 = 960$ elements, 5 mm resolution).
- **Time Integration**: Implicit Newmark-$\beta$ ($\beta = 0.25, \gamma = 0.5$) with Rayleigh damping ($\xi = 0.03$).

---

## 3. Load Cases

1. **Case 1: Quasi-Static Benchmark**: 1500 N steady force at $\dot{\varepsilon} \approx 1.0 \times 10^{-4}\,\text{s}^{-1}$ (laboratory test standard).
2. **Case 2: Dynamic Ramp Loading**: 1500 N applied over 50 ms ($\dot{\varepsilon} \approx 2.2 \times 10^{-2}\,\text{s}^{-1}$, moderate rate hardening).
3. **Case 3: High-Rate Impact Stamping**: 10 ms impact spike (2250 N initial peak decaying to 1500 N) coupled with $-1.0g$ base deceleration ($\dot{\varepsilon} \approx 1.6 \times 10^{-1}\,\text{s}^{-1}$, high strain-rate impact regime).

---

## 4. Summary Validation Table

| Case | Strain Rate $\dot{\varepsilon}_{eq}$ | Peak Stress $\sigma_{vM}$ | Cowper-Symonds $\sigma_y^{dyn}$ | Ductility Retention | Failure Regime |
| :--- | :---: | :---: | :---: | :---: | :--- |
| **Case 1: Static** | $1.0 \times 10^{-4}\,\text{s}^{-1}$ | 143.63 MPa | 250.00 MPa | **100.0%** | Fully ductile, high safety margin |
| **Case 2: Dynamic** | $2.2 \times 10^{-2}\,\text{s}^{-1}$ | 144.66 MPa | 305.47 MPa | **69.9%** | Moderate rate hardening (+22%) |
| **Case 3: Impact** | $1.6 \times 10^{-1}\,\text{s}^{-1}$ | **220.84 MPa** | **333.07 MPa** | **62.8%** | **Dynamic embrittlement risk!** |

> **Key Takeaway**:  
> Static qualification cannot detect dynamic impact vulnerability. While static tests show ample safety margin, high-rate impact elevates operating stresses toward 220 MPa while simultaneously reducing material fracture ductility by over 37%, creating severe risk of unheralded brittle cleavage under emergency vehicle service.

---

## 5. Visualizations

| Constitutive Rate Sensitivity | Transient Rate Evolution | Stress vs Yield Threshold |
| :---: | :---: | :---: |
| ![Constitutive](results/strain_rate_constitutive.png) | ![Evolution](results/transient_rate_evolution.png) | ![Stress Comparison](results/stress_comparison.png) |

---

## 6. How to Run

Prerequisites: Python 3.10+, `dolfinx >= 0.10.0`, `petsc4py`, `numpy`, `matplotlib`.

Run all simulations and generate figures:
```bash
python run_all.py
```

---

## 7. Structure

```
├── README.md                                  # Paper supplementary brief
├── run_all.py                                 # Single entry-point execution
├── geometry/cantilever_beam.geo               # Gmsh geometry
├── mesh/generate_mesh.py                      # Hex8 mesh builder
├── simulations/
│   ├── common.py                              # Elastodynamics & Cowper-Symonds constitutive law
│   ├── static/static_load.py                  # Case 1: Quasi-static benchmark
│   ├── transient/dynamic_load.py              # Case 2: Dynamic ramp loading
│   ├── acceleration/base_acceleration.py      # Case 3: High-rate impact stamping
│   └── postprocess.py                         # Generates figures & summary table
├── results/                                   # Generated figures and JSON data
└── docs/methodology.md                        # Theoretical derivations
```

---

## License
MIT License
