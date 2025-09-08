# IterativeSolver.m — SIMPLE Algorithm for 2D Lid-Driven Cavity (MATLAB)

This solver implements the SIMPLE algorithm with explicit triple-nested loops to solve the classic 2D lid-driven cavity problem in MATLAB.

## Table of Contents

- [Overview](#overview)
- [Key Features](#key-features)
- [Setup and Usage](#setup-and-usage)
- [Main Code Structure](#main-code-structure)
- [Example Output GIFs](#example-output-gifs)
- [Solving Time](#solving-time)
- [Residuals and Convergence](#residuals-and-convergence)
- [License](#license)

---

## Overview

The lid-driven cavity problem is a standard CFD benchmark for incompressible, unsteady flow. This script uses the SIMPLE algorithm (Semi-Implicit Method for Pressure-Linked Equations) on a staggered grid, with explicit iteration for educational clarity.

---

## Key Features

- Finite Volume discretization on staggered grid
- Explicit (triple-nested) loops for clarity
- Real-time animated visualization of:
  - Velocity vectors
  - Velocity magnitude contours
  - Pressure contours
  - Streamlines
  - Residuals
- GIF output for each visualization
- Benchmarking: reports elapsed simulation time

---

## Setup and Usage

1. Open MATLAB (R2020a or newer recommended).
2. Open `IterativeSolver.m` in the `Matlab` folder.
3. Adjust simulation parameters at the top of the script as needed:
   ```matlab
   Re = 100;         % Reynolds number
   n = 151;          % Grid size (n x n)
   dt = 0.0005;      % Time step
   total_time = 2;   % Total simulation time (seconds)
