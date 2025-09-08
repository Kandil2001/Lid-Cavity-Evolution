# Lid Cavity Evolution

*A benchmark suite for unsteady incompressible CFD: From MATLAB fundamentals to industrial applications*

---

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Table of Contents

- [Introduction](#introduction)
- [Why the Lid-Driven Cavity?](#why-the-lid-driven-cavity)
- [Why the SIMPLE Algorithm?](#why-the-simple-algorithm)
- [Governing Equations](#governing-equations)
- [SIMPLE Algorithm Steps](#simple-algorithm-steps)
- [Numerical Methods & Boundary Conditions](#numerical-methods--boundary-conditions)
- [Project Roadmap](#project-roadmap)
- [Benchmark Table](#benchmark-table)
- [Getting Started](#getting-started)
- [Contributing](#contributing)
- [License](#license)
- [References](#references)

---

## 🌟 Introduction

**Lid Cavity Evolution** is an open-source benchmark suite that chronicles the development of the classic lid-driven cavity CFD problem—from foundational MATLAB scripts to industrial-grade solvers. The project emphasizes accuracy, performance, and reproducibility for unsteady incompressible flow simulation.

---

## Why the Lid-Driven Cavity?

- **Standard Test Case:** Simple geometry, well-defined boundaries, and established reference solutions make it ideal for CFD code verification.
- **Rich Physics:** Captures vortex formation, boundary layers, and evolving flow structures.
- **Unsteady Simulation:** Tracks time evolution of flow fields, including transient and nonlinear phenomena.

---

## Why the SIMPLE Algorithm?

Simulating incompressible flows is numerically challenging due to tight coupling of velocity and pressure. The **SIMPLE (Semi-Implicit Method for Pressure-Linked Equations)** algorithm is widely used because it:

- Decouples momentum and continuity equations for robust convergence.
- Uses a staggered grid to prevent checkerboard pressure artifacts.
- Applies under-relaxation for improved stability and efficiency.

---

## 🧮 Governing Equations

**The solver addresses the unsteady, incompressible, 2D Navier-Stokes equations for a square cavity with a moving lid.**

### 1. Continuity (Incompressibility)

$$
\nabla \cdot \mathbf{u} = 0
$$

where $\mathbf{u} = (u, v)$ is the velocity vector.

---

### 2. Momentum (Unsteady, Incompressible, Vector Form)

$$
\frac{\partial \mathbf{u}}{\partial t}
+ (\mathbf{u} \cdot \nabla)\mathbf{u}
= -\nabla p
+ \frac{1}{Re}\nabla^2 \mathbf{u}
$$

- $Re = \frac{UL}{\nu}$ (Reynolds number).  
  $U$: lid velocity, $L$: cavity length, $\nu$: kinematic viscosity.

---

### 3. Expanded 2D Component Form

**x-Momentum:**
$$
\frac{\partial u}{\partial t}
+ u \frac{\partial u}{\partial x}
+ v \frac{\partial u}{\partial y}
= -\frac{\partial p}{\partial x}
+ \frac{1}{Re} \left(
\frac{\partial^2 u}{\partial x^2}
+ \frac{\partial^2 u}{\partial y^2}
\right)
$$

**y-Momentum:**
$$
\frac{\partial v}{\partial t}
+ u \frac{\partial v}{\partial x}
+ v \frac{\partial v}{\partial y}
= -\frac{\partial p}{\partial y}
+ \frac{1}{Re} \left(
\frac{\partial^2 v}{\partial x^2}
+ \frac{\partial^2 v}{\partial y^2}
\right)
$$

---

## 🛠 SIMPLE Algorithm Steps

The SIMPLE algorithm iteratively solves the equations above for incompressible flow.

1. **Predictor Step:**  
   Solve momentum equations for an intermediate velocity $\mathbf{u}^*$, using the current pressure.

2. **Pressure Correction:**  
   Solve the Poisson equation for pressure correction $p'$:
   $$
   \nabla^2 p' = \frac{1}{\Delta t} \nabla \cdot \mathbf{u}^*
   $$

3. **Corrector Step:**  
   Update velocities and pressure:
   $$
   \mathbf{u}^{n+1} = \mathbf{u}^* - \Delta t \nabla p'
   $$
   $$
   p^{n+1} = p^{n} + \alpha p'
   $$
   - $\alpha$: Pressure under-relaxation factor ($0 < \alpha \leq 1$)

---

## 🧑‍💻 Numerical Methods & Boundary Conditions

- **Spatial Discretization:** Second-order central differencing
- **Time Discretization:** First-order implicit Euler
- **Grid:** Staggered (pressure at cell centers, velocities at faces)
- **Boundary Conditions:**
    - **Top lid:** $u = 1$, $v = 0$ (moving wall)
    - **Other walls:** $u = v = 0$ (no-slip)
    - **Pressure:** Neumann ($\frac{\partial p}{\partial n} = 0$) at boundaries

---

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

---

## 📊 Benchmark Table

| Solver                                   | Language      | Paradigm                | Elapsed Time (s) | Speedup | Status         |
|-------------------------------------------|--------------|-------------------------|------------------|---------|---------------|
| SIMPLE2D_LidDrivenCavity                  | MATLAB       | Serial (Loops)          | 2478.76          | 1x      | ✅ Complete    |
| *SimpleLidCavityVector*                   | MATLAB       | Serial (Vectorized)     | TBD              | TBD     | 🚧 In Progress |
| *lid_cavity_serial.py*                    | Python/NumPy | Serial                  | TBD              | TBD     | 📋 Planned     |
| *lid_cavity_parallel.py*                  | Python       | Parallel                | TBD              | TBD     | 📋 Planned     |
| *OpenFOAM Case*                           | OpenFOAM     | Industrial CFD          | TBD              | TBD     | 📋 Planned     |
| *STAR-CCM+ Case*                          | STAR-CCM+    | Commercial CFD          | TBD              | TBD     | 📋 Planned     |

*Hardware: Intel i7-12700K, 32GB RAM*

---

## 🚀 Getting Started

### Prerequisites

- **MATLAB:** R2020a or later (for initial implementations)
- **Python:** 3.8+ (NumPy, planned)
- **OpenFOAM:** (planned)
- **STAR-CCM+:** (planned)

### Running the MATLAB Solver

1. Clone the repository:
    ```bash
    git clone https://github.com/yourusername/lid-cavity-evolution.git
    cd lid-cavity-evolution/matlab
    ```
2. Open `SIMPLE2D_LidDrivenCavity.m` in MATLAB.
3. Run the script and follow prompts (see comments for parameter adjustments).

*Python, OpenFOAM, and STAR-CCM+ instructions will follow with implementation.*

---

## 🤝 Contributing

We welcome all contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines on:

- Adding new solver implementations
- Improving code or documentation
- Adding validation cases
- Reporting issues or suggesting enhancements

Feel free to open an issue or submit a pull request!

---

## 📜 License

This project is licensed under the MIT License.  
See [LICENSE](LICENSE) for details.

---

## 📚 References

1. Ghia, U., Ghia, K. N., & Shin, C. T. (1982). *High-Re solutions for incompressible flow using the Navier-Stokes equations and a multigrid method*. J. Comput. Phys., 48(3), 387-411.
2. Patankar, S. V. (1980). *Numerical Heat Transfer and Fluid Flow*. Hemisphere Publishing.
3. Ferziger, J. H., Perić, M., & Street, R. L. (2002). *Computational Methods for Fluid Dynamics*. Springer.

---

> **Note:** This project is under active development. Check back for updates and new solver releases!
