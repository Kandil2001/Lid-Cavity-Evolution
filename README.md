# Lid Cavity Evolution

A comprehensive benchmark suite tracking the evolution of a CFD solver from fundamental MATLAB code to industrial-grade applications, with a focus on performance comparison and validation.

## ⚡ Performance Highlights

- **41-Minute Baseline:** The initial loop-based MATLAB solver provides a crucial performance baseline, completing in **2478.76 seconds** on a 151x151 grid.
- **10x Speedup Goal:** The vectorized implementation targets a **10x reduction** in solving time.
- **Ultimate Goal:** Demonstrate performance scaling across languages (Python) and paradigms (parallel, industrial solvers).

## 📊 Benchmark Results (Preliminary)

*Results compare solving time for an identical problem setup (Re=100, 151x151 grid, total_time=2s).*

| Solver | Language | Paradigm | Elapsed Time (s) | Speedup | Status |
| :--- | :--- | :--- | :--- | :--- | :--- |
| [SIMPLE2D_LidDrivenCavity](matlab/SIMPLE2D_LidDrivenCavity.m) | MATLAB | Serial (Loops) | 2478.76 | 1x | ✅ Complete |
| *SimpleLidCavityVector* | MATLAB | Serial (Vectorized) | TBD | TBD | 🚧 In Progress |
| *lid_cavity_serial.py* | Python/NumPy | Serial | TBD | TBD | 📋 Planned |
| *lid_cavity_parallel.py* | Python | Parallel | TBD | TBD | 📋 Planned |
| *OpenFOAM Case* | OpenFOAM | Industrial CFD | TBD | TBD | 📋 Planned |
| *STAR-CCM+ Case* | STAR-CCM+ | Commercial CFD | TBD | TBD | 📋 Planned |

**Hardware:** Intel i7-12700K, 32GB RAM

## 🧮 Mathematical Formulation

The lid-driven cavity problem solves the incompressible Navier-Stokes equations:

**Continuity equation:**

∇·u = 0


**Momentum equations:**

∂u/∂t + (u·∇)u = -∇p + (1/Re)∇²u


**In expanded form (2D Cartesian coordinates):**

*Continuity:*

∂u/∂x + ∂v/∂y = 0


*x-Momentum:*

∂u/∂t + u·∂u/∂x + v·∂u/∂y = -∂p/∂x + (1/Re)(∂²u/∂x² + ∂²u/∂y²)


*y-Momentum:*

∂v/∂t + u·∂v/∂x + v·∂v/∂y = -∂p/∂y + (1/Re)(∂²v/∂x² + ∂²v/∂y²)


### SIMPLE Algorithm Steps

The SIMPLE (Semi-Implicit Method for Pressure-Linked Equations) algorithm consists of three main steps:

1. **Predictor Step**: Solve momentum equations for intermediate velocities (u*, v*)

    u* = uⁿ + Δt[-(u·∇)u + (1/Re)∇²u - ∇pⁿ]ᵧ

2. **Pressure Correction**: Solve the pressure Poisson equation

    ∇²p' = (∇·u*)/Δt

    where p' is the pressure correction.

3. **Corrector Step**: Update velocities and pressure

    uⁿ⁺¹ = u* - Δt·∇p'
    pⁿ⁺¹ = pⁿ + α·p'

    where α is the pressure under-relaxation factor (0 < α ≤ 1).

### Discretization Schemes

- **Spatial discretization**: Second-order central differencing
- **Temporal discretization**: First-order Euler implicit
- **Pressure-velocity coupling**: SIMPLE algorithm with under-relaxation
- **Grid**: Staggered grid arrangement for pressure and velocity

**Boundary Conditions:**
- Top wall (lid): u = 1, v = 0 (moving lid)
- Other walls: u = 0, v = 0 (no-slip)
- Pressure: Neumann boundary conditions (∂p/∂n = 0)

## 🗺️ Project Roadmap & Progress

### Phase 1: Foundation in MATLAB
- [x] **Basic Solver**: Loop-based SIMPLE algorithm (`SIMPLE2D_LidDrivenCavity.m`)
- [ ] **Vectorized Solver**: Optimized vectorized implementation (`SimpleLidCavityVector.m`)

### Phase 2: Porting & Optimizing in Python
- [ ] **Serial NumPy**: Direct port of the basic solver to Python/NumPy
- [ ] **Vectorized NumPy**: Leveraging NumPy's optimized array operations
- [ ] **Parallel Implementation**: Using parallel processing (Numba, Dask, or Multiprocessing)

### Phase 3: Industrial Solvers
- [ ] **OpenFOAM**: Case setup and simulation using the open-source OpenFOAM toolkit
- [ ] **STAR-CCM+**: Case setup and simulation using the commercial STAR-CCM+ solver

### Phase 4: Validation & Analysis
- [ ] **Hand Calculations**: Manual analysis for a specific setup to validate the physics
- [ ] **Consolidated Results**: Final comparison of solving time, accuracy, and resource usage across all methods

## 📊 Simulation Results

### Final Flow Field (Re=100)
![Velocity Vectors and Streamlines](https://via.placeholder.com/400x300/FFFFFF/000000?text=Velocity+Field+Placeholder)
*Velocity vectors and streamlines showing primary and secondary vortices*

### Velocity and Pressure Contours
![Velocity Magnitude](https://via.placeholder.com/400x300/FFFFFF/000000?text=Velocity+Contour+Placeholder)
*Velocity magnitude contour with clearly defined boundary layers*

![Pressure Field](https://via.placeholder.com/400x300/FFFFFF/000000?text=Pressure+Field+Placeholder)
*Pressure distribution showing low pressure in the core vortex*

### Convergence History
![Residuals Plot](https://via.placeholder.com/400x300/FFFFFF/000000?text=Residuals+Plot+Placeholder)
*Convergence history showing typical SIMPLE algorithm behavior*

## 🚀 Getting Started

### Prerequisites
- MATLAB R2020a or later
- Python 3.8+ (for future implementations)
- OpenFOAM (for future implementations)
- STAR-CCM+ (for future implementations)


## 🤝 Contributing

We welcome contributions to this benchmark suite! Please see our [Contributing Guidelines](CONTRIBUTING.md) for details on how to:
- Add new solver implementations
- Improve existing code
- Add validation cases
- Report issues or suggest enhancements

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📚 References

1. Ghia, U., Ghia, K. N., & Shin, C. T. (1982). High-Re solutions for incompressible flow using the Navier-Stokes equations and a multigrid method. *Journal of Computational Physics*, 48(3), 387-411.

2. Patankar, S. V. (1980). *Numerical Heat Transfer and Fluid Flow*. Hemisphere Publishing Corporation.

3. Ferziger, J. H., Perić, M., & Street, R. L. (2002). *Computational Methods for Fluid Dynamics*. Springer.

---

**Note:** This project is under active development. Check back regularly for updates and new implementations!
