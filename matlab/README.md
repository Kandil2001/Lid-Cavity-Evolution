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
  <img src="https://github.com/Kandil2001/Lid-Cavity-Evolution/raw/main/assets/velocityiter.gif" width="500"/>
</p>
<p align="center"><i>Velocity field from MATLAB SIMPLE iterative solver at Re = 100 (vectorized output is visually identical)</i></p>

---

## Benchmark Note

- ⚡ **Performance:** The vectorized solver is typically **10x faster** than the loop-based version on larger grids, thanks to MATLAB’s optimized matrix operations.
- Both solvers produce identical results, allowing direct speed and accuracy comparisons.

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
