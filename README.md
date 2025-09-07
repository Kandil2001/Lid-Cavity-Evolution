# Lid Cavity Evolution

A comprehensive benchmark suite tracking the evolution of a CFD solver from fundamental MATLAB code to industrial-grade applications, with a focus on performance comparison and validation.

## 🗺️ Project Roadmap & Progress

### Phase 1: Foundation in MATLAB
- [x] **1_BasicSolver**: Loop-based SIMPLE algorithm (`SimpleLidCavity.m`)
- [x] **2_VectorizedSolver**: Optimized vectorized implementation (`SimpleLidCavityVector.m`)
- **Performance Gain (MATLAB)**: [You will add this later, e.g., "10x speedup"]

### Phase 2: Porting & Optimizing in Python
- [ ] **1_SerialNumPy**: Direct port of the basic solver to Python/NumPy
- [ ] **2_VectorizedNumPy**: Leveraging NumPy's optimized array operations
- [ ] **3_Parallel**: Implementation using parallel processing (e.g., Numba, Dask)

### Phase 3: Industrial Solvers
- [ ] **1_OpenFOAM**: Case setup and simulation using the open-source OpenFOAM toolkit.
- [ ] **2_STARCCM+**: Case setup and simulation using the commercial STAR-CCM+ solver.

### Phase 4: Validation & Analysis
- [ ] **HandCalculations**: Manual analysis for a specific setup to validate the physics.
- [ ] **Consolidated Results**: Final comparison of solving time, accuracy, and resource usage across all methods.

## ⏱️ Benchmark Results (Preliminary)

*Results will be updated as each phase is completed. The goal is to compare solving time for an identical problem setup (Grid, Re, Convergence).*

| Solver | Language/Tool | Paradigm | Avg. Time (s) | Hardware | Status |
| :--- | :--- | :--- | :--- | :--- | :--- |
| [SimpleLidCavity](1_Matlab/1_BasicSolver/) | MATLAB | Serial (Loops) | TBD | CPU | ✅ Complete |
| [SimpleLidCavityVector](1_Matlab/2_VectorizedSolver/) | MATLAB | Serial (Vectorized) | TBD | CPU | ✅ Complete |
| *lid_cavity_serial.py* | Python/NumPy | Serial | TBD | CPU | 🚧 In Progress |
| ... | ... | ... | ... | ... | ... |

## 🤝 Contributing
Ideas and feedback are welcome! Please open an issue to discuss.

## 📜 License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
