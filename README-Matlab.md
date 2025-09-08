# MATLAB Lid-Driven Cavity Solvers

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](../LICENSE)

This directory contains MATLAB implementations of the 2D lid-driven cavity CFD benchmark using the SIMPLE algorithm.

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

## Available Solvers

- [IterativeSolver.m](IterativeSolver.m): Classic SIMPLE algorithm using explicit iteration.
- [VectorizedSolver.m](VectorizedSolver.m): Vectorized MATLAB implementation for improved performance.

Each solver script contains its own documentation and usage instructions in the comments.

## How to Run

1. Open MATLAB (R2020a or newer recommended).
2. Open either `IterativeSolver.m` or `VectorizedSolver.m`.
3. Run the script.
4. Adjust simulation parameters at the top of the script as desired.

## Contributing

Contributions and suggestions are welcome!  
See [CONTRIBUTING.md](../CONTRIBUTING.md) for guidelines.

## License

This code is released under the MIT License.  
See [LICENSE](../LICENSE) for details.

---
