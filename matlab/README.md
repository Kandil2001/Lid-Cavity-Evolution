<p align="center">
  <a href="https://www.python.org/">
    <img src="https://github.com/Kandil2001/Lid-Cavity-Evolution/raw/main/logos/python.png" width="70"/>
  </a>
</p>
<h1 align="center">🐍 Python Lid-Driven Cavity Solvers</h1>
<p align="center"><i>Reference implementations of the SIMPLE algorithm for incompressible flow</i></p>
<p align="center">
  <a href="../LICENSE">
    <img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="MIT License"/>
  </a>
  <a href="https://www.python.org/">
    <img src="https://img.shields.io/badge/Python-3.8+-blue.svg" alt="Python Version"/>
  </a>
  <a href="../CONTRIBUTING.md">
    <img src="https://img.shields.io/badge/Contributions-Welcome-orange.svg" alt="Contributions Welcome"/>
  </a>
</p>

---

## Table of Contents
- [Introduction](#introduction)
- [Why Python?](#why-python)
- [Features](#features)
- [Available Solvers](#available-solvers)
- [How to Run](#how-to-run)
- [Results Preview](#results-preview)
- [Benchmark Note](#benchmark-note)
- [Contributing](#contributing)
- [License](#license)

---

## Introduction

This directory contains reference Python implementations of the SIMPLE algorithm for the classic lid-driven cavity CFD benchmark. Python’s readability and extensive scientific ecosystem make it ideal for learning, research, and rapid prototyping of incompressible flow solvers.

---

## Why Python?

Python is a leading platform for scientific computing and CFD development because:
- 🛠️ **Open source** — Free, widely supported, and cross-platform
- 🧮 **Powerful libraries** — NumPy/SciPy for fast arrays and numerics
- 📈 **Visualization** — Matplotlib, seaborn, and interactive plotting
- 🤝 **Integration** — Easy to interface with C/C++, Fortran, and parallel frameworks
- 💡 **Education** — Clear, concise syntax and a massive support community

---

## Features

- ✅ Loop-based and vectorized SIMPLE solvers (NumPy-based)
- 🚦 Future support for parallel (MPI, OpenMP, Numba/Dask) computation
- 📊 Built-in visualization of velocity and pressure fields
- 🧠 Clear structure for educational use and prototyping
- 🔄 Easy parameter tuning for Reynolds number, grid size, and time steps

---

## Available Solvers

| File/Folder                                                                 | Description                                      |
|------------------------------------------------------------------------------|--------------------------------------------------|
| [IterativeSolver.py](./serial/iterative/IterativeSolver.py)                  | Classic SIMPLE algorithm, explicit Python loops  |
| [VectorizedSolver.py](../serial/vectorized/VectorizedSolver.py)               | Vectorized Python (NumPy) implementation        |
| [MPISolver.py](./parallel/mpi/MPISolver.py)                                  | (Planned) MPI-based parallel solver             |
| [OpenMPSolver.py](./parallel/openmp/OpenMPSolver.py)                         | (Planned) OpenMP-based parallel solver          |

Each solver includes inline documentation and usage instructions.

---

## How to Run

1. **Install requirements** (recommended: use a virtual environment):
    ```bash
    python3 -m venv .venv
    source .venv/bin/activate  # or .venv\Scripts\activate on Windows
    pip install numpy matplotlib
    ```
2. **Run a solver:**
    ```bash
    cd serial/iterative
    python IterativeSolver.py
    ```
    or
    ```bash
    cd serial/vectorized
    python VectorizedSolver.py
    ```
3. **Adjust simulation parameters** at the top of the script as desired.
4. **View outputs:** The script will display velocity and pressure plots automatically.

---

## Results Preview

<p align="center">
  <img src="https://github.com/Kandil2001/Lid-Cavity-Evolution/raw/main/assets/velocitypy.gif" width="500"/>
</p>
<p align="center"><i>Velocity field from Python SIMPLE solver at Re = 100 (output from the vectorized solver; iterative output is visually identical)</i></p>

---

## Benchmark Note

- **Performance:** The vectorized NumPy solver is typically **5–20x faster** than the pure Python loop-based solver for large grids, thanks to optimized array operations.
- Both solvers produce identical results, allowing direct speed and accuracy comparisons.

---

## Contributing

Contributions and suggestions are welcome!  
See [CONTRIBUTING.md](./CONTRIBUTING.md) for guidelines.

---

## License

This code is released under the MIT License.  
See [LICENSE](./LICENSE) for details.
