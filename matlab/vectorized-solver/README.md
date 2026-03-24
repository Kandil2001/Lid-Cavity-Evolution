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

This version preserves the same physical formulation as the iterative solver while replacing loop-based updates with **vectorized array operations**.

- **File:** `VectorizedSolver.m`
- **Location:** `main/matlab/vectorized-solver/`
- **Focus:** Vectorized performance, staggered-grid consistency, and benchmarking

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

* Prevents **checkerboard pressure modes**
* Improves **pressure–velocity coupling**
* Matches **classical CFD formulations**

---

## Key Features

* **True staggered-grid formulation**
* **Fully vectorized field updates**
* Modular solver structure:

  * predictor step
  * pressure Poisson solver
  * corrector step
  * boundary conditions
  * residual tracking
* Continuity residual monitoring
* Direct **GIF writing to disk (no memory overhead)**
* Final summary plot including:

  * velocity vectors
  * streamlines
  * velocity magnitude
  * pressure field
  * vorticity
  * centerline profiles
  * convergence history
* Designed for direct comparison with the iterative solver

---

## Why Vectorized?

* **Faster in MATLAB:** uses optimized matrix operations
* **Cleaner updates:** operates on full fields instead of loops
* **Scalable:** better suited for larger grids
* **Benchmarking:** allows direct comparison vs iterative version
* **Educational:** shows how CFD operators map to array slicing

---

## Usage Instructions

### 1. Requirements

* MATLAB **R2020a or newer**

### 2. Setup

Place the file:

```text
main/matlab/vectorized-solver/VectorizedSolver.m
```

### 3. Configure Parameters

```matlab
Re = 100;
L = 1.0;
N = 51;

dt = 5e-4;
total_time = 1.0;

alpha_u = 0.7;
alpha_p = 0.3;

tol = 1e-6;
max_iter = 200;

poisson_tol = 1e-6;
poisson_max = 800;

record_gif = true;
gif_stride = 1;
```

### 4. Run

```matlab
VectorizedSolver()
```

### 5. Outputs

Saved automatically:

* `vectorized_velocity_vectors.gif`
* `vectorized_velocity_contour.gif`
* `vectorized_pressure_contour.gif`
* `vectorized_streamlines.gif`
* `vectorized_residuals.gif`
* `final_results_vectorized.png`

---

## Numerical Workflow

Each time step:

1. Predictor step (momentum equations)
2. Pressure Poisson solve
3. Velocity correction
4. Pressure update
5. Under-relaxation
6. Boundary conditions
7. Residual + continuity check

---

## Simulation Outputs

### 🌀 Velocity Vectors

<p align="center">
  <img src="https://raw.githubusercontent.com/Kandil2001/Lid-Cavity-Evolution/main/assets/matlab/vectorized_velocity_vectors.gif"/>
</p>

### ⚡ Velocity Magnitude

<p align="center">
  <img src="https://raw.githubusercontent.com/Kandil2001/Lid-Cavity-Evolution/main/assets/matlab/vectorized_velocity_contour.gif"/>
</p>

### 🌊 Streamlines

<p align="center">
  <img src="https://raw.githubusercontent.com/Kandil2001/Lid-Cavity-Evolution/main/assets/matlab/vectorized_streamlines.gif"/>
</p>

### 📊 Pressure Field

<p align="center">
  <img src="https://raw.githubusercontent.com/Kandil2001/Lid-Cavity-Evolution/main/assets/matlab/vectorized_pressure_contour.gif"/>
</p>

### 📉 Residuals

<p align="center">
  <img src="https://raw.githubusercontent.com/Kandil2001/Lid-Cavity-Evolution/main/assets/matlab/vectorized_residuals.gif"/>
</p>

---

## Performance & Convergence

Performance depends on:

* grid size `N`
* Reynolds number
* time step `dt`
* relaxation factors
* Poisson iterations
* GIF frequency

### Recommended starting values

```matlab
Re = 100;
N = 51;
dt = 5e-4;
alpha_u = 0.5;
alpha_p = 0.2;
```

---

## Numerical Notes

This solver is intended for:

* learning CFD algorithms
* debugging
* validation
* benchmarking (iterative vs vectorized)
* presentation-ready visualization

Not intended for industrial CFD use.

---

## Directory Placement

```text
main/
├── matlab/
│   ├── iterative-solver/
│   ├── vectorized-solver/
│   │   ├── VectorizedSolver.m
│   │   └── README.md
```

---

## Contributing

Contributions are welcome.
See `CONTRIBUTING.md`.

---

## License

MIT License — see `LICENSE`.

---

## References

1. Ghia et al. (1982)
2. Patankar (1980)
3. Ferziger & Perić (2002)

---

## Contact

* GitHub Issues
* [kandil.ahmed.amr@gmail.com](mailto:kandil.ahmed.amr@gmail.com)
* LinkedIn

---

**This solver provides a fast, clean, and physically consistent staggered-grid implementation of the lid-driven cavity problem in MATLAB.**

