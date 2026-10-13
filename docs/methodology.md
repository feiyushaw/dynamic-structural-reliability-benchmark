# Methodology: Rate-Dependent Impact Dynamics & Constitutive Degradation

## 1. Governing Elastodynamic Equations

Under dynamic vehicle impact loading, the structural response of the 3D cantilever component is governed by the continuum elastodynamics balance equation:

$$\rho \ddot{\mathbf{u}} - \nabla \cdot \boldsymbol{\sigma}(\mathbf{u}, \dot{\mathbf{u}}) = \mathbf{f}_b$$

where:
- $\mathbf{u}(\mathbf{x}, t)$ is the 3D displacement vector field,
- $\rho = 7850\,\text{kg/m}^3$ is the material density,
- $\mathbf{f}_b(\mathbf{x}, t) = -\rho \mathbf{a}(t)$ is the body force density from vehicle base deceleration,
- $\boldsymbol{\sigma}$ is the stress tensor.

---

## 2. Rate-Dependent Constitutive Law (Cowper-Symonds Model)

In high-rate impact events (such as panic foot stamping), loading rates span 5 to 6 orders of magnitude above quasi-static laboratory standards. Metal plasticity becomes highly rate-sensitive. We model the rate-dependent flow stress via the classic **Cowper-Symonds empirical constitutive equation**:

$$\sigma_y(\dot{\varepsilon}_{eq}) = \sigma_0 \left[ 1 + \left( \frac{\dot{\varepsilon}_{eq}}{D} \right)^{1/p} \right]$$

where:
- $\sigma_0 = 250\,\text{MPa}$ is the quasi-static yield strength at reference rate $\dot{\varepsilon}_0 \sim 10^{-4}\,\text{s}^{-1}$,
- $D = 40.4\,\text{s}^{-1}$ is the material rate sensitivity parameter for structural steel,
- $p = 5.0$ is the strain-rate hardening exponent,
- $\dot{\varepsilon}_{eq}$ is the equivalent von Mises strain rate computed from the strain rate tensor $\dot{\boldsymbol{\varepsilon}}$:

$$\dot{\boldsymbol{\varepsilon}} = \frac{1}{2}\left(\nabla \dot{\mathbf{u}} + (\nabla \dot{\mathbf{u}})^T\right), \quad \dot{\varepsilon}_{eq} = \sqrt{\frac{2}{3} \dot{\boldsymbol{\varepsilon}}' : \dot{\boldsymbol{\varepsilon}}'}$$

where $\dot{\boldsymbol{\varepsilon}}' = \dot{\boldsymbol{\varepsilon}} - \frac{1}{3}\operatorname{tr}(\dot{\boldsymbol{\varepsilon}})\mathbf{I}$ is the deviatoric strain rate tensor.

---

## 3. High Strain-Rate Dynamic Embrittlement & Ductility Loss

Under quasi-static tensile tests ($v = 20\,\text{mm/min}, \dot{\varepsilon} \sim 10^{-4}\,\text{s}^{-1}$), crystal dislocations have ample time to thermally activate and slip, providing superior plastic ductility ($\varepsilon_f \approx 25\% - 30\%$). 

Under high-rate impact ($\dot{\varepsilon} \sim 10^0 - 10^2\,\text{s}^{-1}$), dislocation motion is inhibited by phonon drag, shifting deformation towards twinning and localized adiabatic shearing:
1. **Flow Stress Escalation**: Material yield strength increases significantly ($\sigma_y^{dyn} > 380\,\text{MPa}$);
2. **Dynamic Embrittlement**: The dynamic failure strain $\varepsilon_f(\dot{\varepsilon})$ degrades according to the Johnson-Cook rate-ductility degradation relation:

$$\frac{\varepsilon_f(\dot{\varepsilon})}{\varepsilon_{f0}} = \left[ 1 + C_\varepsilon \ln\left( \max\left(1, \frac{\dot{\varepsilon}_{eq}}{\dot{\varepsilon}_0}\right) \right) \right]^{-1}$$

with $C_\varepsilon \approx 0.08$. The material loses up to $40\% - 50\%$ of its plastic strain capacity.

---

## 4. Why Quasi-Static Tests Fail to Predict Dynamic Reliability

| Metric | Quasi-Static Qualification (QC/T 788) | Emergency Panic Stamping Impact | Physical Mechanism |
| :--- | :---: | :---: | :--- |
| **Loading Velocity** | $20\,\text{mm/min}$ ($0.33\,\text{mm/s}$) | $1.5 - 2.5\,\text{m/s}$ | $\sim 5000\times$ faster impact |
| **Strain Rate $\dot{\varepsilon}_{eq}$** | $\sim 10^{-4}\,\text{s}^{-1}$ | $\sim 1.0 - 5.0\,\text{s}^{-1}$ | Rate jumps across 5 orders of magnitude |
| **Yield Limit $\sigma_y^{dyn}$** | $250.0\,\text{MPa}$ (nominal) | **$380 - 410\,\text{MPa}$** | Dislocation pinned by rate hardening |
| **Ductility Retention** | $100\%$ (fully ductile) | **$< 60\%$ (embrittled)** | Significant drop in fracture toughness $K_{IC}$ |
| **Failure Mode** | Ductile plastic yielding | **Brittle micro-crack initiation** | Zero warning abrupt brittle failure |

A structural component that withstands 3000 N under slow quasi-static testing by virtue of plastic redistribution may undergo sudden brittle cleavage under rapid impact loading because the material constitutive law has been fundamentally altered by the strain rate.
