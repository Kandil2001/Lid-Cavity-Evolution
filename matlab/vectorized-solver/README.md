<p align="center">
  <a href="https://www.mathworks.com/products/matlab.html">
    <img src="https://github.com/Kandil2001/Lid-Cavity-Evolution/raw/main/logos/matlab.png" width="70"/>
  </a>
</p>

<h1 align="center">🌀 SIMPLE2D Lid-Driven Cavity — MATLAB Vectorized Solver</h1>

<p align="center"><i>High-performance true staggered-grid SIMPLE-style pressure-correction solver for unsteady incompressible CFD</i></p>

<p align="center">
  <a href="https://github.com/Kandil2001/Lid-Cavity-Evolution/blob/main/LICENSE">
    <img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="MIT License"/>
  </a>
  <a href="https://github.com/Kandil2001/Lid-Cavity-Evolution/releases">
    <img src="https://img.shields.io/badge/Version-0.2.0-green.svg" alt="Version"/>
  </a>
  <a href="https://github.com/Kandil2001/Lid-Cavity-Evolution/blob/main/CONTRIBUTING.md">
    <img src="https://img.shields.io/badge/Contributions-Welcome-orange.svg" alt="Contributions Welcome"/>
  </a>
  <a href="https://www.mathworks.com/products/matlab.html">
    <img src="https://img.shields.io/badge/MATLAB-R2020a+-blue.svg" alt="MATLAB"/>
  </a>
</p>

---

## Table of Contents

- [About This Solver](#about-this-solver)
- [Grid Arrangement (Staggered Grid)](#grid-arrangement-staggered-grid)
- [Key Features](#key-features)
- [Why Vectorized?](#why-vectorized)
- [Usage Instructions](#usage-instructions)
- [Numerical Workflow](#numerical-workflow)
- [Simulation Outputs](#simulation-outputs)
- [Performance & Convergence](#performance--convergence)
- [Numerical Notes](#numerical-notes)
- [Directory Placement](#directory-placement)
- [Contributing](#contributing)
- [License](#license)
- [References](#references)
- [Contact](#contact)

---

## About This Solver

This branch contains a **vectorized MATLAB implementation** of a **true staggered-grid SIMPLE-style pressure-correction solver** for the 2D lid-driven cavity problem.

It is designed for:

- **high MATLAB performance**
- **algorithm clarity**
- **direct comparison with the iterative solver**
- **clean modular structure**
- **presentation-ready visualization outputs**

This version preserves the same physical problem and output workflow as the iterative solver, while replacing most internal loop-based field updates with vectorized array operations.

- **File:** `VectorizedSolver.m`
- **Location:** `main/matlab/vectorized-solver/`
- **Focus:** Vectorized speed, staggered-grid consistency, modularity, and benchmarking

---

## Grid Arrangement (Staggered Grid)

This solver uses a **true staggered (MAC) grid**.

### Storage layout

- Pressure is stored at **cell centers**
- `u` velocity is stored at **vertical cell faces**
- `v` velocity is stored at **horizontal cell faces**

```matlab
p = zeros(N, N);      % cell centers
u = zeros(N, N+1);    % vertical faces
v = zeros(N+1, N);    % horizontal faces
````

### Why this matters

A staggered grid improves pressure-velocity coupling and avoids the checkerboard pressure issue commonly associated with collocated storage.

It also matches the classical formulation found in CFD references such as:

* Patankar
* Ferziger & Perić
* MAC / SIMPLE-based educational derivations

---

## Key Features

* **True staggered-grid formulation**
* **Vectorized array operations** for improved MATLAB performance
* Modular solver structure:

  * predictor step
  * pressure Poisson solver
  * corrector step
  * boundary condition application
  * residual evaluation
  * GIF export
* Continuity residual monitoring
* Automatic GIF writing directly to disk
* Final summary plot with:

  * velocity vectors
  * streamlines
  * velocity magnitude
  * pressure
  * vorticity
  * centerline profiles
  * convergence history
* Designed for direct benchmarking against the iterative MATLAB solver

---

## Why Vectorized?

* **Faster execution in MATLAB:** vectorized operations are generally more efficient than explicit nested loops
* **Cleaner bulk updates:** whole field regions are updated at once
* **Benchmarking value:** useful for comparing performance against the iterative implementation
* **Educational value:** shows how CFD operators can be expressed naturally with matrix slices
* **Presentation value:** preserves the same output structure while improving runtime behavior

---

## Usage Instructions

### 1. Requirements

* MATLAB **R2020a or newer**

### 2. Setup

Place the file here:

```text
main/matlab/vectorized-solver/VectorizedSolver.m
```

### 3. Configure Parameters

At the top of the file:

```matlab
Re = 100;                % Reynolds number
L = 1.0;                 % Cavity length
N = 51;                  % Number of pressure cells in each direction
dt = 5e-4;               % Time step
total_time = 1.0;        % Total simulated time

alpha_u = 0.7;           % Under-relaxation for velocity
alpha_p = 0.3;           % Under-relaxation for pressure

tol = 1e-6;              % Outer SIMPLE tolerance
max_iter = 200;          % Max SIMPLE iterations per time step

poisson_tol = 1e-6;      % Pressure Poisson tolerance
poisson_max = 800;       % Pressure Poisson max iterations

record_gif = true;       % Record GIFs? true/false
gif_stride = 1;          % Save every gif_stride step
```

### 4. Run the Solver

```matlab
VectorizedSolver()
```

### 5. Generated Outputs

The solver saves the following automatically in the working directory:

* `vectorized_velocity_vectors.gif`
* `vectorized_velocity_contour.gif`
* `vectorized_pressure_contour.gif`
* `vectorized_streamlines.gif`
* `vectorized_residuals.gif`
* `final_results_vectorized.png`

These files are ready to use in reports, presentations, and GitHub documentation.

---

## Numerical Workflow

Each time step follows the same staggered-grid pressure-correction logic as the iterative version:

1. **Predictor Step**

   * Solve the momentum equations for intermediate face velocities

2. **Pressure Poisson Equation**

   * Solve for pressure correction using the divergence of predicted velocities

3. **Corrector Step**

   * Correct the face velocities using pressure correction
   * Update the pressure field

4. **Under-Relaxation**

   * Stabilize iterative convergence

5. **Boundary Conditions**

   * Enforce moving lid and no-slip wall behavior

6. **Residual Evaluation**

   * `u` residual
   * `v` residual
   * pressure residual
   * continuity residual

---

## Simulation Outputs

### 🌀 Velocity Vectors — Flow Field Animation

<p align="center">
  <img src="https://raw.githubusercontent.com/Kandil2001/Lid-Cavity-Evolution/main/assets/matlab/vectorized_velocity_vectors.gif" alt="Velocity Vectors GIF"/>
</p>

### ⚡ Velocity Magnitude — Speed Contours

<p align="center">
  <img src="https://raw.githubusercontent.com/Kandil2001/Lid-Cavity-Evolution/main/assets/matlab/vectorized_velocity_contour.gif" alt="Velocity Magnitude GIF"/>
</p>

### 🌊 Streamlines — Flow Paths

<p align="center">
  <img src="https://raw.githubusercontent.com/Kandil2001/Lid-Cavity-Evolution/main/assets/matlab/vectorized_streamlines.gif" alt="Streamlines GIF"/>
</p>

### 📊 Pressure Distribution

<p align="center">
  <img src="https://raw.githubusercontent.com/Kandil2001/Lid-Cavity-Evolution/main/assets/matlab/vectorized_pressure_contour.gif" alt="Pressure GIF"/>
</p>

### 📉 Residuals — Convergence History

<p align="center">
  <img src="https://raw.githubusercontent.com/Kandil2001/Lid-Cavity-Evolution/main/assets/matlab/vectorized_residuals.gif" alt="Vectorized Residuals GIF"/>
</p>

*GIF and image files are written directly to disk using fixed filenames for easy reuse in presentations and reports.*

---

## Performance & Convergence

Performance depends strongly on:

* pressure grid size `N`
* Reynolds number
* time step `dt`
* under-relaxation factors
* Poisson iteration count
* GIF recording frequency

### Recommended starting values

```matlab
Re = 100;
N = 51;
dt = 5e-4;
alpha_u = 0.5;
alpha_p = 0.2;
```

These are good starting values for stable early testing. After confirming stable behavior, parameters can be increased gradually.

### Example reporting template

| Metric               | Value         |
| -------------------- | ------------- |
| Pressure grid size   | 51 × 51       |
| Reynolds number      | 100           |
| Time step            | 5e-4          |
| Total simulated time | 1.0 s         |
| Elapsed time         | user measured |
| Avg. time per step   | user measured |

---

## Numerical Notes

This solver is intended primarily for:

* education
* CFD algorithm understanding
* debugging
* validation exercises
* iterative vs vectorized benchmarking
* presentation-quality visualization

It is:

* **not** intended as an industrial CFD solver
* **well suited** as a transparent academic MATLAB implementation of a staggered-grid pressure-correction cavity-flow solver

---

## Directory Placement

```text
main/
├── matlab/
│   ├── iterative-solver/
│   │   ├── IterativeSolver.m
│   │   └── README.md
│   ├── vectorized-solver/
│   │   ├── VectorizedSolver.m
│   │   └── README.md
│   └── README.md
├── python/
├── assets/
└── ...
```

See [main README](../../README.md) for project-wide context.

---

## Contributing

Contributions, suggestions, and improvements are welcome.
See [CONTRIBUTING.md](../../CONTRIBUTING.md).

---

## License

Released under the MIT License.
See [LICENSE](../../LICENSE).

---

## References

1. Ghia, U., Ghia, K. N., & Shin, C. T. (1982). *Journal of Computational Physics*, 48(3), 387–411.
2. Patankar, S. V. (1980). *Numerical Heat Transfer and Fluid Flow*.
3. Ferziger, J. H., Perić, M., & Street, R. L. (2002). *Computational Methods for Fluid Dynamics*.

---

## Contact

* [GitHub Issues](https://github.com/Kandil2001/Lid-Cavity-Evolution/issues)
* Email: **[kandil.ahmed.amr@gmail.com](mailto:kandil.ahmed.amr@gmail.com)**
* [LinkedIn](https://www.linkedin.com/in/ahmed-kandil01)

---

**This solver provides a vectorized, modular, and physically consistent staggered-grid pressure-correction implementation for the lid-driven cavity problem in MATLAB, suitable for learning, benchmarking, debugging, and presentation use.**

A few important improvements from your old README:
- it now correctly matches the **true staggered-grid** solver
- it uses `N` instead of the old `n`
- it reflects **direct GIF writing to disk**
- it removes the old code snippets that no longer matched the solver
- it removes the risky hard-coded timing claims and replaces them with a cleaner template

For the image URLs, I also normalized the first GIF link to the raw GitHub format so all of them match.
```
