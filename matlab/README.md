<p align="center">
  <a href="https://www.mathworks.com/products/matlab.html">
    <img src="https://github.com/Kandil2001/Lid-Cavity-Evolution/raw/main/logos/matlab.png" width="70"/>
  </a>
</p>
<h1 align="center">🧮 MATLAB Lid-Driven Cavity Solvers</h1>
<p align="center"><i>Reference implementations of the SIMPLE algorithm for incompressible flow</i></p>
<p align="center">
  <a href="../LICENSE">
    <img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="MIT License"/>
  </a>
  <a href="https://www.mathworks.com/products/matlab.html">
    <img src="https://img.shields.io/badge/MATLAB-R2020a+-blue.svg" alt="MATLAB Version"/>
  </a>
  <a href="../CONTRIBUTING.md">
    <img src="https://img.shields.io/badge/Contributions-Welcome-orange.svg" alt="Contributions Welcome"/>
  </a>
</p>

---

## Table of Contents
- [Introduction](#introduction)
- [Why MATLAB?](#why-matlab)
- [Features](#features)
- [Available Solvers](#available-solvers)
- [How to Run](#how-to-run)
- [Results Preview](#results-preview)
- [Benchmark Note](#benchmark-note)
- [Directory Structure](#directory-structure)
- [Contributing](#contributing)
- [License](#license)

---

## Introduction

The lid-driven cavity problem is a classical benchmark in computational fluid dynamics (CFD) for validating incompressible flow solvers.  
This directory provides reference MATLAB implementations using the SIMPLE algorithm, ideal for learning, rapid prototyping, and establishing baseline results.

---

## Why MATLAB?

MATLAB is a powerful environment for CFD solver development due to:
- ⚡ **Rapid prototyping** — Matrix operations and built-in visualization tools
- 📖 **Clarity** — Clean syntax for learning and debugging
- 📊 **Visualization** — Immediate feedback with plotting capabilities
- 🏗️ **Reference foundation** — Validates algorithm correctness before porting to other languages

---

## Features

- ✅ Loop-based and vectorized SIMPLE solvers
- 📊 Built-in visualization of velocity and pressure fields
- 🧠 Clear structure for educational use and prototyping
- 🔄 Easy parameter tuning for Reynolds number, grid size, and time steps

---

## Available Solvers

| Folder                    | Description                                        |
|---------------------------|----------------------------------------------------|
| `iterative-solver/`       | Loop-based SIMPLE algorithm for clarity and benchmarking |
| `vectorized-solver/`      | High-performance, vectorized MATLAB implementation |

Each solver contains a README and is fully documented.

---

## How to Run

1. Open MATLAB (R2020a or newer recommended)
2. Open either `IterativeSolver.m` or `VectorizedSolver.m` in the appropriate folder
3. Run the script
4. Adjust simulation parameters at the top of the file
5. View velocity and pressure plots generated automatically

---

## Results Preview

<p align="center">
  <img src="https://raw.githubusercontent.com/Kandil2001/Lid-Cavity-Evolution/main/assets/matlab/vectorized_velocity_contour.gif" width="500"/>
</p>
<p align="center"><i>Velocity field from MATLAB SIMPLE iterative solver at Re = 100 (vectorized output is visually identical)</i></p>

---

## Benchmark Note

### Vectorized vs. Iterative Performance Analysis

Based on measured performance, the vectorized solver is currently **significantly slower** than the iterative solver, especially on larger grids, which contradicts the expected benefits of MATLAB vectorization.

| Grid Size | Iterative Elapsed Time | Vectorized Elapsed Time | Performance Change |
| :---: | :---: | :---: | :---: |
| **$51 \times 51$** | $\approx 1126$ s ($\approx 18.8$ min) | $\approx 1266$ s ($\approx 21.1$ min) | **$+12.5\%$ Increase (Slower)** |
| **$151 \times 151$** | $\approx 395$ s ($\approx 6.6$ min) | $\approx 3732$ s ($\approx 62.2$ min) | **$+844\%$ Increase (Much Slower)** |

The performance for the **$151 \times 151$** grid shows a major regression ($\mathbf{844\%}$ slower). A thorough review and profiling of the `vectorized-solver` is required to identify and correct the inefficiency before it can be considered a high-performance implementation.

---

---

## Directory Structure

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

---

## Contributing

Contributions and suggestions are welcome!  
See [CONTRIBUTING.md](../CONTRIBUTING.md) for guidelines.

---

## License

This code is released under the MIT License.  
See [LICENSE](../LICENSE) for details.
