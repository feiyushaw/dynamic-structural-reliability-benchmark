# Methodology & Mathematical Formulation

## 1. Governing Equations

The structural response of the 3D cantilever steel component is modeled via the Navier-Cauchy equations of continuum elastodynamics:

$$\rho \ddot{\mathbf{u}} - \nabla \cdot \boldsymbol{\sigma}(\mathbf{u}) = \mathbf{f}_b$$

where:
- $\mathbf{u}(\mathbf{x}, t)$ is the displacement vector field,
- $\rho = 7850\,\text{kg/m}^3$ is the mass density,
- $\mathbf{f}_b(\mathbf{x}, t)$ is the body force density (representing inertial base acceleration),
- $\boldsymbol{\sigma}$ is the Cauchy stress tensor under isotropic linear elasticity:

$$\boldsymbol{\sigma}(\mathbf{u}) = 2\mu \boldsymbol{\varepsilon}(\mathbf{u}) + \lambda \operatorname{tr}(\boldsymbol{\varepsilon}(\mathbf{u})) \mathbf{I}, \quad \boldsymbol{\varepsilon}(\mathbf{u}) = \frac{1}{2}\left(\nabla \mathbf{u} + (\nabla \mathbf{u})^T\right)$$

Lame parameters:
$$\mu = \frac{E}{2(1+\nu)}, \quad \lambda = \frac{E\nu}{(1+\nu)(1-2\nu)}$$
with $E = 210\,\text{GPa}$ and $\nu = 0.3$.

---

## 2. Semi-Discrete System & Newmark-$\beta$ Time Integration

Using spatial Galerkin finite element discretization on a Hex8 mesh:

$$\mathbf{M} \ddot{\mathbf{u}} + \mathbf{C} \dot{\mathbf{u}} + \mathbf{K} \mathbf{u} = \mathbf{F}_{ext}(t)$$

where Rayleigh damping $\mathbf{C} = \alpha \mathbf{M} + \beta \mathbf{K}$ is calibrated to the beam fundamental bending frequency ($f_1 = 371.3\,\text{Hz}$, damping ratio $\xi = 0.03$).

Time integration is performed using the unconditionally stable implicit Newmark-$\beta$ average acceleration method ($\beta = 0.25, \gamma = 0.5$):

$$\hat{\mathbf{K}} \mathbf{u}_{n+1} = \hat{\mathbf{F}}_{n+1}$$

$$\hat{\mathbf{K}} = \mathbf{K} + \frac{1}{\beta \Delta t^2} \mathbf{M} + \frac{\gamma}{\beta \Delta t} \mathbf{C}$$

$$\hat{\mathbf{F}}_{n+1} = \mathbf{F}_{n+1} + \mathbf{M}\left(\frac{1}{\beta \Delta t^2}\mathbf{u}_n + \frac{1}{\beta \Delta t}\dot{\mathbf{u}}_n + \left(\frac{1}{2\beta}-1\right)\ddot{\mathbf{u}}_n\right) + \mathbf{C}\left(\frac{\gamma}{\beta \Delta t}\mathbf{u}_n + \left(\frac{\gamma}{\beta}-1\right)\dot{\mathbf{u}}_n + \Delta t\left(\frac{\gamma}{2\beta}-1\right)\ddot{\mathbf{u}}_n\right)$$

Dirichlet boundary conditions ($\mathbf{u} = \mathbf{0}$ at $x = 0$) are cleanly enforced on $\hat{\mathbf{K}}$ via diagonal preservation (`zeroRowsColumns`), ensuring exact convergence to the static solution as $\Delta t \to \infty$ or $t \to \infty$.

---

## 3. Analytical Baseline Verification

For an end-loaded cantilever beam under strong-axis bending ($L = 0.3\,\text{m}, H = 0.04\,\text{m}, B = 0.01\,\text{m}, F = 1500\,\text{N}$):

1. **Second moment of area**:
   $$I = \frac{B H^3}{12} = \frac{0.01 \times (0.04)^3}{12} = 5.333 \times 10^{-8}\,\text{m}^4$$

2. **Theoretical static tip deflection**:
   $$w_{theory} = \frac{F L^3}{3 E I} = \frac{1500 \times (0.3)^3}{3 \times 2.1 \times 10^{11} \times 5.333 \times 10^{-8}} = 1.205\,\text{mm}$$
   *FEM computed tip deflection*: **$1.204\,\text{mm}$** (relative error $< 0.1\%$).

3. **Theoretical maximum normal bending stress**:
   $$\sigma_{max, theory} = \frac{M y_{max}}{I} = \frac{(1500 \times 0.3) \times 0.02}{5.333 \times 10^{-8}} = 168.75\,\text{MPa}$$
   *FEM computed volume-averaged element stress*: **$143.63\,\text{MPa}$** (consistent with DG0 Gauss integration point averaging).

---

## 4. Inertial Decoupling Analysis

Under emergency vehicle deceleration $a_{x} = -1.0g = -9.81\,\text{m/s}^2$:
- Ground tire braking friction $\mathbf{F}_{ground} = \mu M_{veh} g$ is contained within the tire-suspension-chassis closed loop.
- The pedal structure is subject only to its own D'Alembert body inertia:
  $$F_{inertia} = m_{beam} \cdot a_{veh} = (\rho \cdot L \cdot H \cdot B) \cdot 9.81 \approx 0.942 \times 9.81 = 9.24\,\text{N}$$
- This inertia accounts for merely **$0.61\%$** of the 1500 N pedal input, confirming that chassis deceleration force does not mechanically cascade to the pedal structure.
