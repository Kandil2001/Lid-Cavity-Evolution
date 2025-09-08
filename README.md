# Lid Cavity Evolution
_A benchmark suite for unsteady incompressible CFD: From MATLAB fundamentals to industrial applications_

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Table of Contents
- [🌟 Introduction](#🌟-introduction)
- [💡 Why the Lid-Driven Cavity?](#💡-why-the-lid-driven-cavity)
- [⚡ Why the SIMPLE Algorithm?](#⚡-why-the-simple-algorithm)
- [🧮 Governing Equations](#🧮-governing-equations)
- [🛠 SIMPLE Algorithm Steps](#🛠-simple-algorithm-steps)
- [🧑‍💻 Numerical Methods & Boundary Conditions](#🧑‍💻-numerical-methods--boundary-conditions)
- [🏁 Project Roadmap](#🏁-project-roadmap)
- [📊 Benchmark Table](#📊-benchmark-table)
- [🚀 Getting Started](#🚀-getting-started)
- [🤝 Contributing](#🤝-contributing)
- [📜 License](#📜-license)
- [📚 References](#📚-references)

## 🌟 Introduction

**Lid Cavity Evolution** is an open-source benchmark suite that chronicles the development of the classic lid-driven cavity CFD problem—from foundational MATLAB scripts to industrial-grade solvers. The project emphasizes accuracy, performance, and reproducibility for unsteady incompressible flow simulation.

## 💡 Why the Lid-Driven Cavity?

- **Standard Test Case:** Simple geometry, well-defined boundaries, and established reference solutions make it ideal for CFD code verification.
- **Rich Physics:** Captures vortex formation, boundary layers, and evolving flow structures.
- **Unsteady Simulation:** Tracks time evolution of flow fields, including transient and nonlinear phenomena.

## ⚡ Why the SIMPLE Algorithm?

Simulating incompressible flows is numerically challenging due to tight coupling of velocity and pressure. The **SIMPLE (Semi-Implicit Method for Pressure-Linked Equations)** algorithm is widely used because it:
- Decouples momentum and continuity equations for robust convergence.
- Uses a staggered grid to prevent checkerboard pressure artifacts.
- Applies under-relaxation for improved stability and efficiency.

## 🧮 Governing Equations

The solver models unsteady, incompressible, two-dimensional flow in a square cavity with a moving lid.

### 1. Continuity Equation (Incompressibility)
$\nabla \cdot \mathbf{u} = 0$
- $\mathbf{u}$: Velocity vector, $\mathbf{u} = (u, v)$
  - $u$: velocity in $x$-direction
  - $v$: velocity in $y$-direction
This equation enforces conservation of mass for incompressible flow: the net flow into any control volume is zero.

### 2. Momentum Equations (Navier-Stokes, Unsteady, 2D)
$\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla)\mathbf{u} = -\nabla p + \frac{1}{Re}\nabla^2 \mathbf{u}$
Where each term means:
- $\frac{\partial \mathbf{u}}{\partial t}$: **Unsteady term** — Time rate of change of velocity.
- $(\mathbf{u} \cdot \nabla)\mathbf{u}$: **Convection term** — Transport of momentum by fluid motion.
- $- \nabla p$: **Pressure gradient term** — Acceleration caused by pressure differences.
- $\frac{1}{Re}\nabla^2 \mathbf{u}$: **Diffusion term** — Viscous spreading of momentum.
- $Re = \frac{UL}{\nu}$: **Reynolds number**  
    - $U$: Lid velocity  
    - $L$: Cavity length  
    - $\nu$: Kinematic viscosity

#### Expanded in 2D Components
- **x-Momentum:**
  $\frac{\partial u}{\partial t} + u \frac{\partial u}{\partial x} + v \frac{\partial u}{\partial y} = -\frac{\partial p}{\partial x} + \frac{1}{Re} \left( \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} \right )$
  - $\frac{\partial u}{\partial t}$: Time derivative of $u$
  - $u \frac{\partial u}{\partial x}$: Convection of $u$ in $x$
  - $v \frac{\partial u}{\partial y}$: Convection of $u$ in $y$
  - $-\frac{\partial p}{\partial x}$: Pressure gradient in $x$
  - $\frac{1}{Re} \frac{\partial^2 u}{\partial x^2}$: Viscous diffusion in $x$
  - $\frac{1}{Re} \frac{\partial^2 u}{\partial y^2}$: Viscous diffusion in $y$
- **y-Momentum:**
  $\frac{\partial v}{\partial t} + u \frac{\partial v}{\partial x} + v \frac{\partial v}{\partial y} = -\frac{\partial p}{\partial y} + \frac{1}{Re} \left( \frac{\partial^2 v}{\partial x^2} + \frac{\partial^2 v}{\partial y^2} \right )$
  - $\frac{\partial v}{\partial t}$: Time derivative of $v$
  - $u \frac{\partial v}{\partial x}$: Convection of $v$ in $x$
  - $v \frac{\partial v}{\partial y}$: Convection of $v$ in $y$
  - $-\frac{\partial p}{\partial y}$: Pressure gradient in $y$
  - $\frac{1}{Re} \frac{\partial^2 v}{\partial x^2}$: Viscous diffusion in $x$
  - $\frac{1}{Re} \frac{\partial^2 v}{\partial y^2}$: Viscous diffusion in $y`

## 🛠 SIMPLE Algorithm Steps

The SIMPLE algorithm solves these equations with the following procedure:
1. **Predictor Step:**  
   - Solve momentum equations for an intermediate velocity $\mathbf{u}^*$ using the current pressure estimate.
2. **Pressure Correction:**  
   - Solve the pressure correction Poisson equation:  
     $\nabla^2 p' = \frac{1}{\Delta t} \nabla \cdot \mathbf{u}^*$  
     where $p'$ is the pressure correction and $\Delta t$ is the time step.
3. **Corrector Step:**  
   - Update velocities and pressure:  
     $\mathbf{u}^{n+1} = \mathbf{u}^* - \Delta t \nabla p'$  
     $p^{n+1} = p^{n} + \alpha p'$  
     where $\alpha$ is the under-relaxation factor ($0 < \alpha \leq 1$).

## 🧑‍💻 Numerical Methods & Boundary Conditions

- **Spatial Discretization:** Second-order central differencing
- **Time Discretization:** First-order implicit Euler
- **Grid:** Staggered (pressure at cell centers, velocities at faces)
- **Boundary Conditions:**
    - **Top lid:** $u = 1$, $v = 0$ (moving wall)
    - **Other walls:** $u = v = 0$ (no-slip)
    - **Pressure:** Neumann ($\frac{\partial p}{\partial n} = 0$) at boundaries

## 🏁 Project Roadmap

| Phase | Description | Status |
|-------|-------------|--------|
| **Phase 1** | MATLAB loop-based SIMPLE solver | ✅ Complete |
|             | Vectorized MATLAB implementation | 🚧 In Progress |
| **Phase 2** | Python/NumPy serial port          | 📋 Planned |
|             | Vectorized NumPy solver           | 📋 Planned |
|             | Parallel Python (Numba/Dask)      | 📋 Planned |
| **Phase 3** | OpenFOAM case setup               | 📋 Planned |
|             | STAR-CCM+ case setup              | 📋 Planned |
| **Phase 4** | Validation & analysis             | 📋 Planned |

## 📊 Benchmark Table

| Solver                                   | Language      | Paradigm                | Elapsed Time (s) | Speedup | Status         |
|-------------------------------------------|--------------|-------------------------|------------------|---------|---------------|
| SIMPLE2D_LidDrivenCavity                  | MATLAB       | Serial (Loops)          | 2478.76          | 1x      | ✅ Complete    |
| _SimpleLidCavityVector_                   | MATLAB       | Serial (Vectorized)     | TBD              | TBD     | 🚧 In Progress |
| _lid_cavity_serial.py_                    | Python/NumPy | Serial                  | TBD              | TBD     | 📋 Planned     |
| _lid_cavity_parallel.py_                  | Python       | Parallel                | TBD              | TBD     | 📋 Planned     |
| _OpenFOAM Case_                           | OpenFOAM     | Industrial CFD          | TBD              | TBD     | 📋 Planned     |
| _STAR-CCM+ Case_                          | STAR-CCM+    | Commercial CFD          | TBD              | TBD     | 📋 Planned     |
_Hardware: Intel i7-12700K, 32GB RAM_

## 🚀 Getting Started

### Prerequisites
- **MATLAB:** R2020a or later (for initial implementations)
- **Python:** 3.8+ (NumPy, planned)
- **OpenFOAM:** (planned)
- **STAR-CCM+:** (planned)

## 🤝 Contributing

We welcome all contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines on:
- Adding new solver implementations
- Improving code or documentation
- Adding validation cases
- Reporting issues or suggesting enhancements
Feel free to open an issue or submit a pull request!

## 📜 License

This project is licensed under the MIT License.  
See [LICENSE](LICENSE) for details.

## 📚 References

1. Ghia, U., Ghia, K. N., & Shin, C. T. (1982). _High-Re solutions for incompressible flow using the Navier-Stokes equations and a multigrid method_. J. Comput. Phys., 48(3), 387-411.
2. Patankar, S. V. (1980). _Numerical Heat Transfer and Fluid Flow_. Hemisphere Publishing.
3. Ferziger, J. H., Perić, M., & Street, R. L. (2002). _Computational Methods for Fluid Dynamics_. Springer.

> **Note:** This project is under active development. Check back for updates and new solver releases!
