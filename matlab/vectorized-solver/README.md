<p align="center">
  <a href="https://www.mathworks.com/products/matlab.html">
    <img src="https://github.com/Kandil2001/Lid-Cavity-Evolution/raw/main/logos/matlab.png" width="70"/>
  </a>
</p>
<h1 align="center">🌀 SIMPLE2D Lid-Driven Cavity — MATLAB Vectorized Solver Branch</h1>
<p align="center"><i>High-performance vectorized SIMPLE algorithm for unsteady incompressible CFD</i></p>
<p align="center">
  <a href="https://github.com/Kandil2001/Lid-Cavity-Evolution/blob/main/LICENSE">
    <img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="MIT License"/>
  </a>
  <a href="https://github.com/Kandil2001/Lid-Cavity-Evolution/releases">
    <img src="https://img.shields.io/badge/Version-0.1.0-green.svg" alt="Version"/>
  </a>
  <a href="https://github.com/Kandil2001/Lid-Cavity-Evolution/blob/main/CONTRIBUTING.md">
    <img src="https://img.shields.io/badge/Contributions-Welcome-orange.svg" alt="Contributions Welcome"/>
  </a>
  <a href="https://www.mathworks.com/products/matlab.html">
    <img src="https://img.shields.io/badge/MATLAB-R2020a+-blue.svg" alt="MATLAB"/>
  </a>
</p>

---

## Table of Contents

- [About This Solver](#about-this-solver)
- [Key Features](#key-features)
- [Why Vectorized?](#why-vectorized)
- [Usage Instructions](#usage-instructions)
- [Key Code Snippets](#key-code-snippets)
- [Simulation Outputs](#simulation-outputs)
- [Performance & Convergence](#performance--convergence)
- [Directory Placement](#directory-placement)
- [Contributing](#contributing)
- [License](#license)
- [References](#references)
- [Contact](#contact)

---

## About This Solver

This branch contains a **fully vectorized MATLAB implementation** of the SIMPLE algorithm for the lid-driven cavity problem.  
It is designed for direct benchmarking against the iterative (loop-based) version, using **identical structure, boundary conditions, and outputs** for fair comparison and maximum MATLAB performance.

- **File:** `VectorizedSolver.m`
- **Location:** `main/matlab/vectorized-solver/`
- **Focus:** Vectorized speed, educational clarity, and modular structure.

---

## Key Features

- Finite volume discretization on a staggered grid
- **Fully vectorized** operations for superior MATLAB performance
- Modular function structure for maintainability and clarity
- Real-time visualization and optional GIF export
- Residual tracking, convergence monitoring, and performance reporting
- Identical interface and outputs to the iterative solver for direct benchmarking

---

## Why Vectorized?

- **Maximum speed in MATLAB:** Exploits MATLAB’s optimization for array/matrix operations, drastically reducing runtime for large grids.
- **Direct benchmarking:** Uses the same interface and output as the loop version, so you can compare speed and accuracy directly.
- **Educational value:** Demonstrates how classic CFD algorithms can be efficiently mapped to MATLAB’s strengths.

---

## Usage Instructions

1. **Requirements:** MATLAB R2020a or newer.
2. **Setup:**  
   Place `VectorizedSolver.m` in `main/matlab/vectorized-solver/`.
3. **Configure Parameters:**  
   At the top of the file:
   ```matlab
   Re = 100;                % Reynolds number
   n = 51;                  % Grid size
   dt = 0.002;              % Time step
   total_time = 1.0;        % Simulation time (s)
   alpha_u = 0.7;           % Under-relaxation for velocity
   alpha_p = 0.3;           % Under-relaxation for pressure
   tol = 1e-6;              % SIMPLE inner iteration tolerance
   max_iter = 300;          % Max SIMPLE iterations per time step
   record_gif = true;       % Enable GIF recording
    ```
4. **Run the Script:**  
   - In MATLAB, enter `VectorizedSolver()` in the Command Window or run the script file.
   - The simulation will begin, displaying progress and performance metrics in the terminal.
   - Animated figures and GIFs for each scene (velocity vectors, magnitude, pressure, streamlines, residuals) are generated automatically.

5. **Check Results:**  
   - After completion, find GIFs (e.g., `vectorized_velocity_vectors.gif`, `vectorized_pressure_contour.gif`, etc.) and a summary plot (`final_results.png`) in your working directory.
   - The summary plot includes velocity vectors, streamlines, velocity magnitude, pressure field, vorticity, centerline velocity profiles, and convergence history.

---

## Key Code Snippets

Below are essential excerpts from `VectorizedSolver.m`, demonstrating the modular and transparent structure of this solver.

### 1. User-Adjustable Simulation Parameters

```matlab
% Set physical and numerical parameters
Re = 100;                % Reynolds number
n = 51;                  % Grid size
dt = 0.002;              % Time step
total_time = 1.0;        % Simulation time
alpha_u = 0.7;           % Under-relaxation (velocity)
alpha_p = 0.3;           % Under-relaxation (pressure)
tol = 1e-6;              % Convergence tolerance
max_iter = 300;          % Max SIMPLE iterations per time step
record_gif = true;       % Enable GIF creation
```

### 2. Main Time-Stepping and SIMPLE Iteration Loop

```matlab
while time < total_time
    step = step + 1;
    for iter = 1:max_iter
        u_prev = u; v_prev = v; p_prev = p;
        [u_star, v_star] = predictor_step_vectorized(u, v, p, dx, dy, dt, nu, alpha_u);
        p_prime = solve_pressure_poisson_vectorized(u_star, v_star, dx, dy, dt, tol, max_iter);
        [u, v, p] = corrector_step_vectorized(u_star, v_star, p, p_prime, dx, dy, dt, alpha_p);
        % Boundary conditions
        u(:,1) = 0; u(:,end) = 0; u(1,:) = 0; u(end,:) = 1;
        v(:,1) = 0; v(:,end) = 0; v(1,:) = 0; v(end,:) = 0;
        p(1,:) = p(2,:); p(end,:) = p(end-1,:);
        p(:,1) = p(:,2); p(:,end) = p(:,end-1);
        % Residuals/convergence
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

### 3. Vectorized Predictor Step

```matlab
function [u_star, v_star] = predictor_step_vectorized(u, v, p, dx, dy, dt, nu, alpha)
n = size(u,1);
ii = 2:n-1; jj = 2:n-1;
% u-momentum: vectorized update
du2dx = ((u(jj,ii) + u(jj,ii+1)).^2 - (u(jj,ii-1) + u(jj,ii)).^2) / (4*dx);
duvdy = ((v(jj,ii) + v(jj,ii+1)) .* (u(jj,ii) + u(jj+1,ii)) ...
        - (v(jj-1,ii) + v(jj-1,ii+1)) .* (u(jj-1,ii) + u(jj,ii))) / (4*dy);
d2udx2 = (u(jj,ii+1) - 2*u(jj,ii) + u(jj,ii-1)) / dx^2;
d2udy2 = (u(jj+1,ii) - 2*u(jj,ii) + u(jj-1,ii)) / dy^2;
dpdx = (p(jj,ii+1) - p(jj,ii)) / dx;
u_star = u;
u_star(jj,ii) = u(jj,ii) + alpha * dt * (-du2dx - duvdy - dpdx + nu*(d2udx2 + d2udy2));
% v-momentum: similar vectorized update
% ...
end
```

### 4. Vectorized Pressure Poisson Solver

```matlab
function p_prime = solve_pressure_poisson_vectorized(u_star, v_star, dx, dy, dt, tol, max_iter)
n = size(u_star,1);
ii = 2:n-1; jj = 2:n-1;
p_prime = zeros(n);
for iter = 1:max_iter
    p_old = p_prime;
    rhs = ((u_star(jj,ii) - u_star(jj,ii-1))/dx + (v_star(jj,ii) - v_star(jj-1,ii))/dy) / dt;
    p_prime(jj,ii) = 0.25 * (p_prime(jj,ii+1) + p_prime(jj,ii-1) ...
                            + p_prime(jj+1,ii) + p_prime(jj-1,ii) - dx^2 * rhs);
    % Neumann BCs
    p_prime(1,:) = p_prime(2,:);
    p_prime(end,:) = p_prime(end-1,:);
    p_prime(:,1) = p_prime(:,2);
    p_prime(:,end) = p_prime(:,end-1);
    if max(max(abs(p_prime - p_old))) < tol
        break;
    end
end
end
```

### 5. Vectorized Corrector Step

```matlab
function [u, v, p] = corrector_step_vectorized(u_star, v_star, p, p_prime, dx, dy, dt, alpha)
n = size(p,1);
ii = 2:n-1; jj = 2:n-1;
u = u_star; v = v_star;
u(jj,ii) = u_star(jj,ii) - alpha * dt / dx * (p_prime(jj,ii+1) - p_prime(jj,ii));
v(jj,ii) = v_star(jj,ii) - alpha * dt / dy * (p_prime(jj+1,ii) - p_prime(jj,ii));
p = p + alpha * p_prime;
end
```

### 6. GIF Creation and Post-Processing

```matlab
function create_gifs(gif_scenes)
    scenes = fieldnames(gif_scenes);
    for i = 1:length(scenes)
        scene_name = scenes{i};
        scene_data = gif_scenes.(scene_name);
        if ~isempty(scene_data.frames)
            fprintf('Creating %s...\n', scene_data.filename);
            for j = 1:length(scene_data.frames)
                [A, map] = rgb2ind(scene_data.frames(j).cdata, 256);
                if j == 1
                    imwrite(A, map, scene_data.filename, 'gif', ...
                           'LoopCount', Inf, 'DelayTime', 0.1);
                else
                    imwrite(A, map, scene_data.filename, 'gif', ...
                           'WriteMode', 'append', 'DelayTime', 0.1);
                end
            end
            fprintf('Saved: %s\n', scene_data.filename);
        end
    end
end
```

## Simulation Outputs

### 🌀 Velocity Vectors — Flow Field Animation
![Velocity Vectors GIF](https://github.com/Kandil2001/Lid-Cavity-Evolution/blob/efb371647516be5995487e8319755c7db985d059/assets/matlab/vectorized_velocity_vectors.gif)

---

### ⚡ Velocity Magnitude — Speed Contours
![Velocity Magnitude GIF](https://raw.githubusercontent.com/Kandil2001/Lid-Cavity-Evolution/main/assets/matlab/vectorized_velocity_contour.gif)

---

### 🌊 Streamlines — Flow Paths
![Streamlines GIF](https://raw.githubusercontent.com/Kandil2001/Lid-Cavity-Evolution/main/assets/matlab/vectorized_streamlines.gif)

---

### 📊 Pressure Distribution
![Pressure GIF](https://raw.githubusercontent.com/Kandil2001/Lid-Cavity-Evolution/main/assets/matlab/vectorized_pressure_contour.gif)


*GIF/image files are saved in the script directory. Filenames are customizable in the code.*

---

## Performance and Convergence

| Metric                | Typical Value                  |
|-----------------------|-------------------------------|
| Elapsed Time          | ~1266,38 seconds (~21,11 minutes)|
| Avg. Time per Step    | ~2,5328 seconds                |


### 📉 Residuals — Convergence History
![Residuals GIF](https://raw.githubusercontent.com/Kandil2001/Lid-Cavity-Evolution/main/assets/matlab/vectorized_residuals.gif)

---

## Directory Placement

This solver is part of the MATLAB branch:
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
See [main README](../../README.md) for project-wide context.

---

## Contributing

Contributions, suggestions, and improvements are encouraged!  
Refer to [CONTRIBUTING.md](../../CONTRIBUTING.md).

---

## License

Released under the MIT License.  
See [LICENSE](../../LICENSE).

---

## References

1. Ghia, U., Ghia, K. N., & Shin, C. T. (1982). J. Comput. Phys., 48(3), 387-411.
2. Patankar, S. V. (1980). _Numerical Heat Transfer and Fluid Flow_.
3. Ferziger, J. H., Perić, M., & Street, R. L. (2002). _Computational Methods for Fluid Dynamics_.

---

## Contact

- [GitHub Issues](https://github.com/Kandil2001/Lid-Cavity-Evolution/issues)
- Email: **kandil.ahmed.amr@gmail.com**
- [LinkedIn](https://www.linkedin.com/in/ahmed-kandil01)

---

**This branch provides the baseline for MATLAB implementations in the Lid Cavity Evolution suite.  
For advanced solvers, validations, and industrial CFD, see the main project documentation.**
