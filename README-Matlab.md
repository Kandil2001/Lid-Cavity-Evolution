# MATLAB Lid-Driven Cavity Solvers

This directory contains MATLAB implementations of the 2D lid-driven cavity CFD benchmark using the SIMPLE algorithm.

[![License: MIT](../LICENSE)]

## Table of Contents

- [Introduction](#introduction)
- [Why MATLAB?](#why-matlab)
- [Available Solvers](#available-solvers)
- [How to Run](#how-to-run)
- [Contributing](#contributing)
- [License](#license)

## Introduction

The lid-driven cavity problem is a classical benchmark in computational fluid dynamics (CFD) for validating incompressible flow solvers. This directory provides reference MATLAB implementations using the SIMPLE algorithm, suitable for learning, rapid prototyping, and establishing baseline results.

## Why MATLAB?

MATLAB is an excellent environment for CFD solver development because:

- **Rapid prototyping:** Easy matrix operations and built-in visualization.
- **Clarity:** Clean syntax helps learning and debugging.
- **Visualization:** Immediate feedback with built-in plotting tools.
- **Reference foundation:** Validates algorithm correctness before porting to other languages or optimizing for performance.

Starting with MATLAB allows contributors to focus on the core algorithm and compare results easily.

## Available Solvers

- **Iterative Solver:**  
  Main script: [`matlab/IterativeSolver.m`](IterativeSolver.m)  
  README: [`matlab/IterativeSolver_README.md`](IterativeSolver_README.md) (if you create this, link here)

- **Vectorized Solver:**  
  Main script: [`matlab/VectorizedSolver.m`](VectorizedSolver.m)  
  README: [`matlab/VectorizedSolver_README.md`](VectorizedSolver_README.md) (if you create this, link here)

> Each solver has its own script and (optionally) its own README describing usage, features, and simulation parameters.

## How to Run

1. Open MATLAB (R2020a or newer recommended).
2. Navigate to the `matlab/` directory.
3. Open either `IterativeSolver.m` or `VectorizedSolver.m`.
4. Run the script (`F5` or type the script name in the command window).
5. Adjust simulation parameters at the top of the script if desired.
6. For details, check the comments inside each script or the corresponding README.

## Contributing

Contributions and suggestions are welcome!  
See the main [CONTRIBUTING.md](../CONTRIBUTING.md) for guidelines.

## License

This code is released under the MIT License.  
See the [LICENSE](../LICENSE) for details.

---

**For solver-specific guidance, see each solver's README or script comments.**
