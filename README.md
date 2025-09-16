# SIMPLE2D_LidDrivenCavity_Optimized.m — MATLAB Iterative SIMPLE Solver

**Branch of the Lid Cavity Evolution Benchmark Suite**  
_Educational, loop-based, and benchmark implementation of the SIMPLE algorithm for unsteady 2D incompressible flow in a square cavity._

---

## Table of Contents

- [About This Solver](#about-this-solver)
- [Key Features](#key-features)
- [Algorithm Summary](#algorithm-summary)
- [Numerical Methods & Boundary Conditions](#numerical-methods--boundary-conditions)
- [Usage Instructions](#usage-instructions)
- [Simulation Outputs](#simulation-outputs)
- [Performance & Convergence](#performance--convergence)
- [Directory Placement](#directory-placement)
- [Contributing](#contributing)
- [License](#license)
- [References](#references)
- [Contact](#contact)

---

## About This Solver

This branch contains a **modular, loop-based MATLAB implementation** of the SIMPLE algorithm for the classic lid-driven cavity problem.  
It serves as the **educational and benchmarking baseline** for more advanced solvers in the [Lid Cavity Evolution](../../README.md) project.

- **File:** `IterativeSolver.m`
- **Location:** `main/matlab/iterative-solver/`
- **Focus:** Algorithm clarity, modularity, and benchmark performance.

---

## Key Features

| Feature                                   | Description                                                    |
|--------------------------------------------|----------------------------------------------------------------|
| Finite Volume, Staggered Grid              | Pressure at cell centers, velocities at faces                  |
| Modular Functions                          | Predictor, corrector, and pressure Poisson steps               |
| Triple-nested loops (no vectorization)     | Transparent, easy to follow algorithm structure                |
| Precomputed constants & preallocated arrays| Efficient memory and calculation management                    |
| Animated visualization & GIF export        | Each variable/scene gets its own GIF; titles include step info |
| Performance reporting                      | Elapsed time and per-step diagnostics                          |
| Strict convergence criteria                | Residual tracking and early stopping                           |
| Final summary plots                        | Velocity, pressure, vorticity, centerline profiles, residuals  |

---

## Algorithm Summary

- **Predictor Step:** Intermediate velocities solved using current pressure.
- **Pressure Correction:** Pressure Poisson equation solved iteratively.
- **Corrector Step:** Velocity and pressure updated.
- **Residual Monitoring:** $u$, $v$, and $p$ residuals checked for convergence.
- **Boundary conditions:** Set after each field update.
- **Each scene is captured in its own figure for GIF creation.**

---

## Numerical Methods & Boundary Conditions

| Aspect                 | Approach / Value                                    |
|------------------------|----------------------------------------------------|
| Spatial Discretization | Second-order central differencing                  |
| Temporal Discretization| First-order implicit Euler                         |
| Grid                   | Staggered, square ($n \times n$)                   |
| Lid (top wall)         | $u=1$, $v=0$ (moving lid)                          |
| Other walls            | $u=v=0$ (no-slip)                                  |
| Pressure BCs           | Homogeneous Neumann ($\partial p/\partial n = 0$)  |

---

## Usage Instructions

1. **Requirements:** MATLAB R2020a or newer.
2. **Setup:**  
   - Place `IterativeSolver.m` in `main/matlab/iterative-solver/`.
   - (Optional) Asset folder for GIF/image outputs.
3. **Configure Parameters:**  
   At the top of the file:
   ```matlab
   Re = 100;           % Reynolds number (try 100, 400, 1000)
   L = 1.0;            % Cavity length
   n = 51;             % Grid size (n x n, e.g., 31/51/101)
   dt = 0.002;         % Time step
   total_time = 1.0;   % Simulated time (seconds)
   alpha_u = 0.7;      % Under-relaxation for velocity (0.5 - 0.8 typical)
   alpha_p = 0.3;      % Under-relaxation for pressure (0.2 - 0.5 typical)
   tol = 1e-6;         % SIMPLE inner iteration tolerance
   max_iter = 300;     % Max SIMPLE inner iterations per time step
   record_gif = true;  % Record GIFs for each scene
    ```
4. **Run the Script:**  
   - In MATLAB, enter `IterativeSolver()` in the Command Window or run the script file.
   - The simulation will begin, displaying progress and performance metrics in the terminal.
   - Animated figures and GIFs for each scene (velocity vectors, magnitude, pressure, streamlines, residuals) are generated automatically.

5. **Check Results:**  
   - After completion, find GIFs (e.g., `iterative_velocity_vectors.gif`, `iterative_pressure_contour.gif`, etc.) and a summary plot (`final_results.png`) in your working directory.
   - The summary plot includes velocity vectors, streamlines, velocity magnitude, pressure field, vorticity, centerline velocity profiles, and convergence history.

---

## Simulation Outputs

| Output Type           | Description                | Example Filename                   |
|-----------------------|---------------------------|------------------------------------|
| Velocity Vectors GIF  | Flow field animation      | `iterative_velocity_vectors.gif`   |
| Velocity Magnitude GIF| Speed contours            | `iterative_velocity_contour.gif`   |
| Pressure GIF          | Pressure distribution     | `iterative_pressure_contour.gif`   |
| Streamlines GIF       | Flow paths                | `iterative_streamlines.gif`        |
| Residuals GIF         | Convergence history       | `iterative_residuals.gif`          |
| Final Results Figure  | Comprehensive summary     | `final_results.png`                |

*GIF/image files are saved in the script directory. Filenames are customizable in the code.*

---

## Performance & Convergence

| Metric                | Typical Value                  |
|-----------------------|-------------------------------|
| Elapsed Time          | ~36,435 seconds (~607 minutes)|
| Avg. Time per Step    | ~72.87 seconds                |
| Final $u$ Residual    | $1.11 \times 10^{-16}$        |
| Final $v$ Residual    | $5.55 \times 10^{-17}$        |
| Final $p$ Residual    | $1.06 \times 10^{-1}$         |

Residuals and convergence are visualized and logged at each time step.

---

## Directory Placement

This solver is part of the MATLAB branch:
```
main/
├── matlab/
│ ├── iterative-solver/
│ │ ├── IterativeSolver.m
│ │ └── README.md
│ ├── vectorized-solver/
│ │ ├── VectorizedSolver.m
│ │ └── README.md
│ └── README.md
├── python/
│ ├── ...
├── assets/
└── ...
```
See [main README](../../README.md) for project-wide context.

---

## Contributing

Contributions, suggestions, and improvements are encouraged!  
Refer to [CONTRIBUTING.md](../../CONTRIBUTING.md).

---

## License

Released under the MIT License.  
See [LICENSE](../../LICENSE).

---

## References

1. Ghia, U., Ghia, K. N., & Shin, C. T. (1982). J. Comput. Phys., 48(3), 387-411.
2. Patankar, S. V. (1980). _Numerical Heat Transfer and Fluid Flow_.
3. Ferziger, J. H., Perić, M., & Street, R. L. (2002). _Computational Methods for Fluid Dynamics_.

---

## Contact

- [GitHub Issues](https://github.com/Kandil2001/Lid-Cavity-Evolution/issues)
- Email: **kandil.ahmed.amr@gmail.com**
- [LinkedIn](https://www.linkedin.com/in/ahmed-kandil01)

---

**This branch provides the baseline for MATLAB implementations in the Lid Cavity Evolution suite.  
For advanced solvers, validations, and industrial CFD, see the main project documentation.**
