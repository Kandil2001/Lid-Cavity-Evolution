
<p align="center">
  <a href="https://www.python.org/">
    <img src="https://github.com/Kandil2001/Lid-Cavity-Evolution/raw/main/logos/python.png" width="70"/>
  </a>
</p>
<h1 align="center">🐍 Python Serial SIMPLE Solvers</h1>
<p align="center"><i>Reference implementations of the SIMPLE algorithm for the lid-driven cavity problem — single-core (serial) versions</i></p>
<p align="center">
  <a href="../../../LICENSE">
    <img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="MIT License"/>
  </a>
  <a href="https://www.python.org/">
    <img src="https://img.shields.io/badge/Python-3.8+-blue.svg" alt="Python Version"/>
  </a>
  <a href="../../../CONTRIBUTING.md">
    <img src="https://img.shields.io/badge/Contributions-Welcome-orange.svg" alt="Contributions Welcome"/>
  </a>
</p>

---

## Table of Contents
- [Introduction](#introduction)
- [Features](#features)
- [Available Solvers](#available-solvers)
- [How to Run](#how-to-run)
- [Results Preview](#results-preview)
- [Benchmark Note](#benchmark-note)
- [Contributing](#contributing)
- [License](#license)

---

## Introduction

This folder contains **serial (single-core) Python implementations** of the SIMPLE algorithm for the lid-driven cavity CFD benchmark.  
The solvers here were originally developed in their own branches (`iterative` and `vectorized`) and are now organized under a single `serial` branch for easier access and comparison.

The **iterative** solver uses explicit Python loops for clarity, while the **vectorized** solver uses NumPy array operations for higher performance.

---

## Features

- ✅ **Iterative (loop-based)** solver for educational clarity
- ✅ **Vectorized (NumPy)** solver for speed
- 📊 Built-in visualization of velocity, pressure, and residual convergence
- 🧠 Modular code structure with clear functions for each SIMPLE step
- 🔄 Easy parameter tuning (Reynolds number, grid size, time step)
- 📝 Inline documentation

---

## Available Solvers

| File                                                                                  | Description                                      |
|---------------------------------------------------------------------------------------|--------------------------------------------------|
| [`IterativeSolver.py`](./iterative/IterativeSolver.py)                      | Loop-based SIMPLE algorithm (explicit iteration) |
| [`VectorizedSolver.py`](./vectorized/VectorizedSolver.py)                  | Fully vectorized NumPy implementation for speed  |

---

## How to Run

1. **Install requirements** (in a virtual environment is recommended):
    ```bash
    python3 -m venv .venv
    source .venv/bin/activate  # Windows: .venv\Scripts\activate
    pip install numpy matplotlib
    ```

2. **Run the iterative solver**:
    ```bash
    cd iterative
    python IterativeSolver.py
    ```

3. **Run the vectorized solver**:
    ```bash
    cd ../vectorized
    python VectorizedSolver.py
    ```

4. **Adjust parameters**:  
   At the top of each script, you can change:
   ```python
   Re = 100
   n = 51
   dt = 0.002
   total_time = 1.0
   record_gif = True
   ```

## View Results

After running a solver, the program will:
- Display velocity vector plots
- Display pressure contours
- Show residual convergence plots in real-time

If `record_gif = True`:
- GIF animations of the simulation will be saved in the current working directory.
- GIFs will include velocity vectors, velocity magnitude contours, pressure fields, streamlines, and residual plots.

## Results Preview

<p align="center">
  <img src="https://github.com/Kandil2001/Lid-Cavity-Evolution/raw/main/assets/velocitypy.gif" width="500"/>
</p>
<p align="center"><i>Velocity field from the Python vectorized solver at Re = 100 (iterative output is visually identical)</i></p>

## Benchmark Note

- ⚡ The vectorized solver is typically **5–20× faster** than the loop-based version on large grids.
- Both solvers produce identical results — the vectorized version simply runs faster by leveraging optimized NumPy array operations.
- Use the iterative solver if you want to understand the algorithm step-by-step, or the vectorized solver for larger simulations and faster runtimes.

## Contributing

Contributions and suggestions are welcome!  
See [CONTRIBUTING.md](../../../CONTRIBUTING.md) for guidelines.

Ways to contribute:
- Add new solver variants
- Improve performance or readability
- Enhance visualizations
- Report bugs or suggest features

## License

This code is released under the MIT License.  
See [LICENSE](../../../LICENSE) for details.

You are free to:
- Use the code in personal or commercial projects
- Modify and redistribute it
- Provided you include the original license and attribution
