# MATLAB Lid-Driven Cavity Solvers

This directory contains MATLAB implementations of the 2D lid-driven cavity CFD benchmark using the SIMPLE algorithm.

## Table of Contents

- [Introduction](#introduction)
- [Why MATLAB?](#why-matlab)
- [Available Solvers](#available-solvers)
- [How to Run](#how-to-run)
- [Contributing](#contributing)
- [License](#license)

## Introduction

The lid-driven cavity problem is a classical CFD benchmark for validating incompressible flow solvers. This directory provides MATLAB reference implementations for rapid prototyping and educational purposes.

## Why MATLAB?

MATLAB is an ideal environment for CFD solver development because:

- **Rapid prototyping:** Easy matrix operations and built-in visualization.
- **Clarity:** Clean syntax aids learning and debugging.
- **Visualization:** Built-in plotting for immediate feedback.
- **Reference foundation:** Validates methods before porting to other languages.

Starting in MATLAB helps establish correct algorithms and reference results for future optimized or parallel versions.

## Available Solvers

- [Iterative-Solver-Matlab](Iterative-Solver-Matlab/README.md): Classic SIMPLE algorithm with explicit loops for educational clarity.
- [Vectorized-Solver-Matlab](Vectorized-Solver-Matlab/README.md): Vectorized MATLAB implementation for improved performance and scalability.

> **Each solver directory contains its own README with usage details, features, and parameters.**

## How to Run

1. Open MATLAB (R2020a or newer recommended).
2. Navigate to the solver directory of your choice.
3. Open the main `.m` file and run it (press `F5` or type the script name in the command window).
4. Adjust simulation parameters at the top of the script if desired.

## Contributing

Contributions and suggestions are welcome!  
See the main [CONTRIBUTING.md](../CONTRIBUTING.md) for guidelines.

## License

This code is released under the MIT License.  
See the [LICENSE](../LICENSE) for details.

---

**For solver-specific guidance, see each solver's README in its folder.**
