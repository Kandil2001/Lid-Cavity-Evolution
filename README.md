# SIMPLE2D_LidDrivenCavity_Vectorized.m — Vectorized SIMPLE Solver

[View the MATLAB source code: SIMPLE2D_LidDrivenCavity_Vectorized.m](Vectorized-Solver.m)

This script implements the SIMPLE algorithm for unsteady 2D lid-driven cavity flow on a staggered grid, using a fully vectorized approach for maximum MATLAB performance. Structure, boundary conditions, and outputs match the iterative version for direct comparison. Real-time visualization and GIF output are included.

## Table of Contents
- [Overview](#overview)
- [Key Features](#key-features)
- [Why Vectorized?](#why-vectorized)
- [How to Use](#how-to-use)
- [Key Code Snippets](#key-code-snippets)
- [Simulation Output](#simulation-output)
- [Performance and Running Time](#performance-and-running-time)
- [Residuals and Convergence](#residuals-and-convergence)
- [Contributing](#contributing)
- [License](#license)

## Overview

This solver is a high-performance, vectorized MATLAB implementation of the SIMPLE algorithm for the classic lid-driven cavity CFD benchmark. It is designed to match the structure and outputs of the iterative (loop-based) solver for fair benchmarking, while demonstrating MATLAB's array operation advantages.

## Key Features

- Finite Volume discretization on a staggered grid
- Fully vectorized matrix operations for maximum speed in MATLAB
- Modular function structure for clarity and extension
- Real-time animated visualization (velocity, pressure, streamlines, residuals)
- Optional GIF export for all scenes
- Identical boundary conditions, plotting, and output metrics as the iterative version
- Elapsed time and per-step performance reporting

## Why Vectorized?

This implementation leverages MATLAB’s strength in vector and matrix operations to dramatically improve computation speed over explicit for-loops. Major operations are expressed as array calculations, minimizing loop overhead and enabling larger grid sizes and faster testing.

- **Faster than loops:** MATLAB is optimized for vector/matrix operations, making this approach highly efficient.
- **Direct benchmarking:** Structure and outputs are kept identical to the iterative version, so users can compare speed and accuracy directly.
- **Educational:** Shows how classic CFD algorithms can be efficiently mapped to MATLAB’s strengths.

## How to Use

1. Open MATLAB (R2020a or newer recommended).
2. Open `SIMPLE2D_LidDrivenCavity_Vectorized.m` in this folder.
3. Adjust simulation parameters at the top as needed:
    ```matlab
    Re = 100;           % Reynolds number
    n = 151;            % Grid size
    dt = 0.0005;        % Time step
    total_time = 2;     % Total simulation time (seconds)
    record_gif = true;  % Set to true to record GIFs
    ```
4. Run the script.
5. Check figures and terminal output for convergence and performance.

## Key Code Snippets

**Initialization:**
```matlab
% Simulation parameters
Re = 100;                  % Reynolds number
L = 1.0;                   % Cavity size
n = 151;                   % Grid size (n x n)
dx = L/(n-1); dy = dx;
dt = 0.0005;               % Time step
nu = 1/Re;                 % Kinematic viscosity
total_time = 2;            % Total simulation time

% Grid and variables
[X, Y] = meshgrid(0:dx:L, 0:dy:L);
u = zeros(n); v = zeros(n); p = zeros(n);
u(end,:) = 1;              % Moving lid BC
```
**Main Time-Stepping and SIMPLE Loop:**
```matlab
while time < total_time
    step = step + 1;

    % SIMPLE inner iteration (vectorized)
    for iter = 1:max_iter
        [u_star, v_star] = predictor_step_vec(u, v, p, dx, dy, dt, nu, alpha_u);
        p_prime = solve_pressure_poisson_vec(u_star, v_star, dx, dy, dt, tol, max_iter);
        [u, v, p] = corrector_step_vec(u_star, v_star, p, p_prime, dx, dy, dt, alpha_p);

        % Check convergence
        res_u = max(max(abs(u - u_prev)));
        res_v = max(max(abs(v - v_prev)));
        res_p = max(max(abs(p - p_prev)));
        if max([res_u, res_v, res_p]) < tol
            break;
        end
    end

    % Visualization, GIF capture, etc.
    time = time + dt;
end
```
**GIF Creation:**
```matlab
if record_gif && mod(step, gif_frame_interval) == 0
    subplot(2,3,1); gif_scenes.velocity_vectors.frames = [gif_scenes.velocity_vectors.frames, getframe(gcf)];
    subplot(2,3,2); gif_scenes.velocity_contour.frames = [gif_scenes.velocity_contour.frames, getframe(gcf)];
    subplot(2,3,3); gif_scenes.pressure_contour.frames = [gif_scenes.pressure_contour.frames, getframe(gcf)];
    subplot(2,3,5); gif_scenes.residuals.frames = [gif_scenes.residuals.frames, getframe(gcf)];
    subplot(2,3,4); gif_scenes.streamlines.frames = [gif_scenes.streamlines.frames, getframe(gcf)];
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

## Contributing

Contributions, suggestions, and improvements are welcome!  
If you’d like to propose changes or edits, please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines on how to participate.

## License

This code is released under the MIT License.  
See [LICENSE](LICENSE) for details.


*For questions, feedback, or contributions, see the main project repository.*
