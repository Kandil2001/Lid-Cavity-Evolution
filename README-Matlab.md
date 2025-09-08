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

- [Iterative-Solver-Matlab](Iterative-Solver-Matlab/README.md): Classic SIMPLE algorithm using explicit loops for educational clarity and benchmarking.
- [Vectorized-Solver-Matlab](Vectorized-Solver-Matlab/README.md): Vectorized MATLAB implementation for improved performance and scalability.

> **Each solver directory contains its own README with usage details, features, and simulation parameters.**  
> *Please ensure you use the exact folder names as shown above (case-sensitive).*

## How to Run

1. Open MATLAB (R2020a or newer recommended).
2. Navigate to the solver directory (`Iterative-Solver-Matlab` or `Vectorized-Solver-Matlab`).
3. Open the main `.m` file and run it (press `F5` or type the script name in the command window).
4. Adjust simulation parameters at the top of the script if desired.
5. Refer to the solver's README for details on outputs and visualization.

## Contributing

Contributions and suggestions are welcome!  
See the main [CONTRIBUTING.md](../CONTRIBUTING.md) for guidelines.

## License

This code is released under the MIT License.  
See the [LICENSE](../LICENSE) for details.

---
