# Lid Cavity Evolution

**A Benchmark Suite for Incompressible CFD: From MATLAB Fundamentals to Industrial Applications**

---

## 🌟 Introduction

**Lid Cavity Evolution** is an open-source benchmark suite that chronicles the transformation of a classic CFD (Computational Fluid Dynamics) problem—the 2D lid-driven cavity—from basic MATLAB scripts to advanced, industrial-grade implementations. This project is designed not only to compare performance and accuracy across languages, algorithms, and platforms, but also to provide a transparent, reproducible roadmap for anyone interested in CFD solver development, numerical methods, or scientific benchmarking.

### Why the Lid-Driven Cavity?

The lid-driven cavity is a canonical test case in fluid dynamics, valued for its simplicity and rich physical behavior. Its well-defined geometry and boundary conditions make it ideal for validating new solvers, testing algorithms, and benchmarking computational performance.

### Why Benchmark Solver Evolution?

- **Performance Insight:** Understand how coding paradigms (loops, vectorization, parallelization) and language choices (MATLAB, Python, C++, industrial tools) impact speed and scalability.
- **Numerical Reliability:** Compare solution accuracy, convergence, and stability across implementations.
- **Educational Value:** Trace the journey from naïve code to robust, validated engineering tools, demystifying the process for students and practitioners.

---

## 💡 Motivation: The SIMPLE Algorithm

Solving incompressible flow presents unique challenges, especially in coupling pressure and velocity fields. Direct approaches often fail due to numerical instability, slow convergence, or spurious pressure oscillations ("checkerboarding"). The **SIMPLE (Semi-Implicit Method for Pressure-Linked Equations)** algorithm was developed to elegantly overcome these barriers:

- **Pressure-Velocity Coupling:** SIMPLE decouples the solution of momentum and continuity equations, enabling robust convergence for incompressible flows.
- **Staggered Grid Implementation:** Avoids non-physical oscillations by properly arranging pressure and velocity variables.
- **Stability & Scalability:** The iterative predictor-corrector approach, with under-relaxation, ensures stability even on fine grids or high Reynolds numbers.

SIMPLE is the backbone of countless academic and industrial solvers—making it the perfect framework for this benchmark.

---

## 🧮 Mathematical Formulation

- **Continuity:** ∇·u = 0
- **Momentum:** ∂u/∂t + (u·∇)u = -∇p + (1/Re)∇²u

**SIMPLE Steps:**
1. **Predictor:** Solve for intermediate velocities with the current pressure field.
2. **Pressure Correction:** Solve a Poisson equation for pressure correction.
3. **Corrector:** Update velocities and pressure to enforce mass conservation.

**Discretization & Boundaries:**
- **Spatial:** Second-order central differencing  
- **Temporal:** First-order Euler implicit  
- **Grid:** Staggered arrangement  
- **Boundary Conditions:** Moving lid (u=1, v=0), no-slip walls, pressure Neumann BC.

---

## 🏁 Project Roadmap

**Phase 1: MATLAB Foundation**  
- [x] Basic loop-based SIMPLE solver  
- [ ] Vectorized MATLAB implementation  

**Phase 2: Python Evolution**  
- [ ] Serial NumPy port  
- [ ] Vectorized NumPy solver  
- [ ] Parallelized Python (Numba/Dask)  

**Phase 3: Industrial CFD**  
- [ ] OpenFOAM case setup  
- [ ] STAR-CCM+ case setup  

**Phase 4: Validation & Analysis**  
- [ ] Manual hand calculations  
- [ ] Consolidated comparison of performance, accuracy, and resource usage  

---

## 📊 Benchmark Table

| Solver | Language | Paradigm | Elapsed Time (s) | Speedup | Status |
| :--- | :--- | :--- | :--- | :--- | :--- |
| [SIMPLE2D_LidDrivenCavity](matlab/SIMPLE2D_LidDrivenCavity.m) | MATLAB | Serial (Loops) | 2478.76 | 1x | ✅ Complete |
| *SimpleLidCavityVector* | MATLAB | Serial (Vectorized) | TBD | TBD | 🚧 In Progress |
| *lid_cavity_serial.py* | Python/NumPy | Serial | TBD | TBD | 📋 Planned |
| *lid_cavity_parallel.py* | Python | Parallel | TBD | TBD | 📋 Planned |
| *OpenFOAM Case* | OpenFOAM | Industrial CFD | TBD | TBD | 📋 Planned |
| *STAR-CCM+ Case* | STAR-CCM+ | Commercial CFD | TBD | TBD | 📋 Planned |

**Hardware:** Intel i7-12700K, 32GB RAM

---

## 🚀 Getting Started

### Prerequisites
- MATLAB R2020a or later
- Python 3.8+ (for future implementations)
- OpenFOAM (for future implementations)
- STAR-CCM+ (for future implementations)

---

## 🤝 Contributing

We welcome all contributions! See [Contributing Guidelines](CONTRIBUTING.md) for how to:
- Add new solver implementations
- Improve code or documentation
- Add validation cases
- Report issues or suggest enhancements

---

## 📜 License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.

---

## 📚 References

1. Ghia, U., Ghia, K. N., & Shin, C. T. (1982). *High-Re solutions for incompressible flow using the Navier-Stokes equations and a multigrid method*. Journal of Computational Physics, 48(3), 387-411.
2. Patankar, S. V. (1980). *Numerical Heat Transfer and Fluid Flow*. Hemisphere Publishing Corporation.
3. Ferziger, J. H., Perić, M., & Street, R. L. (2002). *Computational Methods for Fluid Dynamics*. Springer.

---

**Note:** This project is under active development. Check back regularly for updates!
