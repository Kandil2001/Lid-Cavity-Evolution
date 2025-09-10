# SIMPLE2D_LidDrivenCavity_Optimized.m — Optimized Iterative SIMPLE Solver

This script implements the SIMPLE algorithm for unsteady 2D lid-driven cavity flow on a staggered grid, using an optimized, pure loop-based approach for clarity and performance. Real-time visualization and GIF output are included.

## Table of Contents
- [Overview](#overview)
- [Key Features](#key-features)
- [Why Is This Implementation Optimized?](#why-is-this-implementation-optimized)
- [How to Use](#how-to-use)
- [Key Code Snippets](#key-code-snippets)
- [Simulation Output](#simulation-output)
- [Performance and Running Time](#performance-and-running-time)
- [Residuals and Convergence](#residuals-and-convergence)
- [License](#license)

## Overview

This solver is a performance-optimized, non-vectorized implementation of the SIMPLE algorithm for the canonical lid-driven cavity CFD benchmark. It is designed for benchmarking and for transparent understanding of each algorithmic step.

## Key Features

- Finite Volume discretization on a staggered grid
- Optimized triple-nested loops (no vectorization)
- Minimal memory overhead and precomputed constants
- Real-time animated visualization (velocity, pressure, streamlines, residuals)
- Optional GIF export for all scenes
- Elapsed time and per-step performance reporting

## Why Is This Implementation Optimized?

This solver is crafted for maximum efficiency within the constraints of MATLAB’s loop-based coding. Optimizations include:

- **Precomputing constants** (e.g., $1/dx$, $dt \cdot \alpha$) outside inner loops to avoid redundant calculations.
- **Minimal memory overhead:** All arrays are preallocated and expanded only if necessary.
- **Tight triple-nested loops:** Loops are structured for clarity and cache efficiency, with no unnecessary function calls or dynamic allocation inside the core.
- **No vectorization:** While MATLAB is fastest with vectorized code, this is the fastest possible pure-loop version—a direct performance baseline for educational use and comparison with more advanced approaches.

This approach is ideal for benchmarking, teaching, and as a reference point before moving to vectorized or compiled solvers.

## How to Use

1. Open MATLAB (R2020a or newer recommended).
2. Open `SIMPLE2D_LidDrivenCavity_Optimized.m` in the `Matlab` folder.
3. Adjust simulation parameters at the top as needed:
    ```matlab
    Re = 100;           % Reynolds number
    n = 151;            % Grid size
    dt = 0.001;         % Time step
    total_time = 2;     % Total simulation time (seconds)
    record_gif = true;  % Set to true to record GIFs
    ```
4. Run the script.
5. Check figures and terminal output for convergence and performance.

## Key Code Snippets

**Initialization:**
```matlab
Re = 100; L = 1.0; n = 151;
dx = L/(n-1); dt = 0.001;
u = zeros(n); v = zeros(n); p = zeros(n);
u(end,:) = 1; % Lid boundary condition
```
**Main Time-Stepping and SIMPLE Loop:**
```matlab
while time < total_time
    for iter = 1:max_iter
        [u_star, v_star] = predictor_step_fast(u, v, p, dx, dy, dt, nu, alpha_u);
        p_prime = solve_pressure_poisson_fast(u_star, v_star, dx, dy, dt, tol, max_iter);
        [u, v, p] = corrector_step_fast(u_star, v_star, p, p_prime, dx, dy, dt, alpha_p);
        % Check convergence
    end
    % Visualization & GIF saving
    time = time + dt;
end
```
**GIF Creation:**
```matlab
if record_gif && mod(step, gif_frame_interval) == 0
    subplot(2,3,1); gif_scenes.velocity_vectors.frames = [gif_scenes.velocity_vectors.frames, getframe(gcf)];
    subplot(2,3,2); gif_scenes.velocity_contour.frames = [gif_scenes.velocity_contour.frames, getframe(gcf)];
    % ...repeat for other scenes...
end
```

## Simulation Output

**Velocity Vectors Animation:**  
![Velocity Vectors](velocity_vectors.gif)

**Velocity Magnitude Contour:**  
![Velocity Magnitude](velocity_contour.gif)

**Pressure Contour:**  
![Pressure Field](pressure_contour.gif)

**Residuals Convergence:**  
![Residuals](residuals.gif)

*Replace the filenames above with your actual GIF/image filenames if different.  
To display your GIFs, upload them to the same folder as this README.*

## Performance and Running Time

- **Total elapsed real time:**  

Elapsed real time: ___ seconds (___ minutes)

- **Average time per time step:**  

Average time per step: ___ seconds

## Residuals and Convergence

- **Residual history plot:**  
![Residuals Convergence](residuals.gif)

- **Final residual values:**  

Final residuals: u = ___, v = ___, p = ___


## License

This code is released under the MIT License.  
See [LICENSE](../LICENSE) for details.


*For questions, feedback, or contributions, see the main project repository.*
