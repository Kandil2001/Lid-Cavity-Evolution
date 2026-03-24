````
<p align="center">
  <a href="https://www.mathworks.com/products/matlab.html">
    <img src="https://github.com/Kandil2001/Lid-Cavity-Evolution/raw/main/logos/matlab.png" width="70"/>
  </a>
</p>

<h1 align="center">🌀 SIMPLE2D Lid-Driven Cavity — MATLAB Iterative Solver</h1>

<p align="center">
<i>Educational true staggered-grid SIMPLE-style pressure-correction solver for unsteady incompressible CFD</i>
</p>

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

## 📚 Table of Contents

- [About This Solver](#about-this-solver)
- [Grid Arrangement (Staggered Grid)](#grid-arrangement-staggered-grid)
- [Key Features](#key-features)
- [Why Iterative?](#why-iterative)
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

## 🧠 About This Solver

This branch contains an **iterative (loop-based) MATLAB implementation** of a  
**true staggered-grid SIMPLE-style pressure-correction solver** for the lid-driven cavity problem.

It is designed for:

- 📖 **Educational clarity**
- 🔍 **Algorithm transparency**
- 🧪 **Debugging and experimentation**
- ⚖️ **Benchmarking vs vectorized solvers**
- 🎥 **Presentation-ready visualization (GIFs)**

---

## 🧩 Grid Arrangement (Staggered Grid)

This solver uses a **true staggered (MAC) grid**, which is the standard in CFD.

### Storage:

- Pressure → **cell centers**
- u-velocity → **vertical faces**
- v-velocity → **horizontal faces**

```matlab
p = zeros(N, N);      % cell centers
u = zeros(N, N+1);    % vertical faces
v = zeros(N+1, N);    % horizontal faces
````

### Why staggered grid?

* Eliminates **checkerboard pressure problem**
* Improves **pressure–velocity coupling**
* Matches **classical SIMPLE formulation**
* Standard in textbooks (Patankar, Ferziger, etc.)

---

## 🚀 Key Features

* ✅ True **staggered-grid formulation**
* 🔁 Fully **loop-based (iterative)** implementation
* 🧱 Modular structure:

  * predictor step
  * pressure Poisson solver
  * corrector step
  * boundary conditions
* 📊 Residual tracking:

  * velocity
  * pressure
  * continuity (mass conservation)
* 🎞️ Automatic **GIF generation** (saved directly to disk)
* 📈 Final summary plots and convergence history
* 🧪 Ideal for **learning, debugging, and validation**

---

## ❓ Why Iterative?

* **Maximum clarity**: every step is explicit
* **Easy debugging**: inspect values at any location
* **Educational value**: understand SIMPLE deeply
* **Baseline comparison**: vs vectorized solver
* **Flexible experimentation**: modify physics easily

---

## ⚙️ Usage Instructions

### 1. Requirements

* MATLAB **R2020a or newer**

---

### 2. Setup

Place the solver file:

```
main/matlab/iterative-solver/IterativeSolver.m
```

---

### 3. Configure Parameters

At the top of the file:

```matlab
Re = 100;                % Reynolds number
L = 1.0;                 % Cavity length
N = 51;                  % Pressure grid size
dt = 5e-4;               % Time step
total_time = 1.0;        % Simulation time

alpha_u = 0.7;           % Velocity under-relaxation
alpha_p = 0.3;           % Pressure under-relaxation

tol = 1e-6;              % SIMPLE tolerance
max_iter = 200;          % Max SIMPLE iterations

poisson_tol = 1e-6;      
poisson_max = 800;

record_gif = true;       
gif_stride = 1;
```

---

### 4. Run the Solver

```matlab
IterativeSolver()
```

---

### 5. Outputs

Generated automatically:

* `iterative_velocity_vectors.gif`
* `iterative_velocity_contour.gif`
* `iterative_pressure_contour.gif`
* `iterative_streamlines.gif`
* `iterative_residuals.gif`
* `final_results.png`

✔ Ready to use in reports and presentations

---

## 🔄 Numerical Workflow

Each time step follows:

1. **Predictor Step**

   * Solve momentum equations for intermediate velocities

2. **Pressure Poisson Equation**

   * Enforce incompressibility:
     [
     \nabla \cdot \vec{u} = 0
     ]

3. **Corrector Step**

   * Update velocity using pressure correction
   * Update pressure field

4. **Under-Relaxation**

   * Stabilizes convergence

5. **Boundary Conditions**

   * Lid-driven top wall
   * No-slip on all walls

6. **Residual Evaluation**

   * Velocity residuals
   * Pressure residual
   * **Continuity residual (important)**

---

## 🎥 Simulation Outputs

### 🌀 Velocity Vectors

<p align="center">
  <img src="https://raw.githubusercontent.com/Kandil2001/Lid-Cavity-Evolution/main/assets/matlab/iterative_velocity_vectors.gif"/>
</p>

### ⚡ Velocity Magnitude

<p align="center">
  <img src="https://raw.githubusercontent.com/Kandil2001/Lid-Cavity-Evolution/main/assets/matlab/iterative_velocity_contour.gif"/>
</p>

### 🌊 Streamlines

<p align="center">
  <img src="https://raw.githubusercontent.com/Kandil2001/Lid-Cavity-Evolution/main/assets/matlab/iterative_streamlines.gif"/>
</p>

### 📊 Pressure Field

<p align="center">
  <img src="https://raw.githubusercontent.com/Kandil2001/Lid-Cavity-Evolution/main/assets/matlab/iterative_pressure_contour.gif"/>
</p>

### 📉 Residuals

<p align="center">
  <img src="https://raw.githubusercontent.com/Kandil2001/Lid-Cavity-Evolution/main/assets/matlab/iterative_residuals.gif"/>
</p>

---

## 📊 Performance & Convergence

Performance depends on:

* Grid size (`N`)
* Reynolds number
* Time step (`dt`)
* Relaxation factors
* Poisson iterations
* GIF frequency

### Recommended starting values:

```matlab
Re = 100;
N = 51;
dt = 5e-4;
alpha_u = 0.5;
alpha_p = 0.2;
```

Then increase gradually.

---

## ⚠️ Numerical Notes

This solver is designed for:

* education
* understanding SIMPLE
* debugging CFD algorithms
* visualization and presentation

It is:

❌ Not an industrial CFD solver
✅ A **transparent academic implementation**

---

## 📂 Directory Placement

```
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

---

## 🤝 Contributing

Contributions are welcome!

See: [CONTRIBUTING.md](../../CONTRIBUTING.md)

---

## 📜 License

MIT License
See: [LICENSE](../../LICENSE)

---

## 📚 References

1. Ghia et al. (1982), *J. Comput. Phys.*
2. Patankar (1980), *Numerical Heat Transfer and Fluid Flow*
3. Ferziger & Perić (2002), *CFD*

---

## 📬 Contact

* GitHub Issues:
  [https://github.com/Kandil2001/Lid-Cavity-Evolution/issues](https://github.com/Kandil2001/Lid-Cavity-Evolution/issues)
* Email: **[kandil.ahmed.amr@gmail.com](mailto:kandil.ahmed.amr@gmail.com)**
* LinkedIn: [https://www.linkedin.com/in/ahmed-kandil01](https://www.linkedin.com/in/ahmed-kandil01)

---

**This solver provides a clean, modular, and physically consistent implementation of a staggered-grid pressure-correction method for the lid-driven cavity problem — ideal for learning, debugging, and presentation use.**

```

---

## 🔥 Final note (important)

This README is now:

✅ technically correct  
✅ aligned with your new staggered solver  
✅ strong for GitHub  
✅ strong for presentation  
✅ looks professional  

If you want next level (for thesis/interview), I can:
- add **validation vs Ghia benchmark plots**
- add **Re=100, 400, 1000 comparison**
- add **grid convergence study section**

That would make this repo **very strong academically** 💯
```
