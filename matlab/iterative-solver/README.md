<p align="center">
  <a href="https://www.mathworks.com/products/matlab.html">
    <img src="https://github.com/Kandil2001/Lid-Cavity-Evolution/raw/main/logos/matlab.png" width="70"/>
  </a>
</p>
<h1 align="center">🌀 SIMPLE2D Lid-Driven Cavity — MATLAB Iterative Solver Branch</h1>
<p align="center"><i>Educational loop-based SIMPLE algorithm for unsteady incompressible CFD</i></p>
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
- [Why Iterative?](#why-iterative)
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

This branch contains an **iterative (loop-based) MATLAB implementation** of the SIMPLE algorithm for the lid-driven cavity problem.  
It is designed for educational clarity and direct benchmarking against the vectorized version, using **identical structure, boundary conditions, and outputs** for fair comparison and algorithm understanding.

- **File:** `IterativeSolver.m`
- **Location:** `main/matlab/iterative-solver/`
- **Focus:** Educational clarity, algorithm transparency, and performance benchmarking

---

## Key Features

- Finite volume discretization on a staggered grid
- **Triple-nested loop implementation** for maximum educational clarity
- Modular function structure for maintainability and understanding
- Real-time visualization and optional GIF export
- Residual tracking, convergence monitoring, and performance reporting

---

## Why Iterative?

- **Educational clarity:** Explicit loops make the algorithm steps transparent and easy to follow
- **Debugging friendly:** Step-by-step execution allows for easy inspection of intermediate values
- **Baseline performance:** Provides reference timing for vectorized optimization comparisons
- **Algorithm understanding:** Helps students and researchers understand the fundamental SIMPLE steps
- **Foundation for learning:** Serves as a starting point for understanding CFD algorithm implementation

---

## Usage Instructions

1. **Requirements:** MATLAB R2020a or newer.
2. **Setup:**  
   Place `IterativeSolver.m` in `main/matlab/iterative-solver/`.
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
   - In MATLAB, enter `IterativeSolver()` in the Command Window or run the script file.
   - The simulation will begin, displaying progress and performance metrics in the terminal.
   - Animated figures and GIFs for each scene (velocity vectors, magnitude, pressure, streamlines, residuals) are generated automatically.

5. **Check Results:**  
   - After completion, find GIFs (e.g., `iterative_velocity_vectors.gif`, `iterative_pressure_contour.gif`, etc.) and a summary plot (`final_results.png`) in your working directory.
   - The summary plot includes velocity vectors, streamlines, velocity magnitude, pressure field, vorticity, centerline velocity profiles, and convergence history.

---

## Key Code Snippets

Below are essential excerpts from `IterativeSolver.m`, demonstrating the modular and transparent structure of this solver.

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
        [u_star, v_star] = predictor_step(u, v, p, dx, dy, dt, nu, alpha_u);
        p_prime = solve_pressure_poisson(u_star, v_star, dx, dy, dt, tol, max_iter);
        [u, v, p] = corrector_step(u_star, v_star, p, p_prime, dx, dy, dt, alpha_p);
        % Apply boundary conditions and check convergence
        u(:,1) = 0; u(:,end) = 0; u(1,:) = 0; u(end,:) = 1; % Lid
        v(:,1) = 0; v(:,end) = 0; v(1,:) = 0; v(end,:) = 0; % No-slip
        p(1,:) = p(2,:); p(end,:) = p(end-1,:);
        p(:,1) = p(:,2); p(:,end) = p(:,end-1);
        res_u = max(max(abs(u - u_prev)));
        res_v = max(max(abs(v - v_prev)));
        res_p = max(max(abs(p - p_prev)));
        if max([res_u, res_v, res_p]) < tol
            break;
        end
    end
    % Advance simulation time, record results/GIFs...
    time = time + dt;
end
```

### 3. Predictor Step (Momentum Equations)

```matlab
function [u_star, v_star] = predictor_step(u, v, p, dx, dy, dt, nu, alpha)
n = size(u,1);
u_star = u; v_star = v;
for j = 2:n-1
    for i = 2:n-1
        du2dx = ((u(j,i)+u(j,i+1))^2 - (u(j,i-1)+u(j,i))^2) / (4*dx);
        duvdy = ((v(j,i)+v(j,i+1))*(u(j,i)+u(j+1,i)) ...
                -(v(j-1,i)+v(j-1,i+1))*(u(j-1,i)+u(j,i))) / (4*dy);
        d2udx2 = (u(j,i+1)-2*u(j,i)+u(j,i-1)) / dx^2;
        d2udy2 = (u(j+1,i)-2*u(j,i)+u(j-1,i)) / dy^2;
        dpdx = (p(j,i+1) - p(j,i)) / dx;
        u_star(j,i) = u(j,i) + alpha * dt * (-du2dx - duvdy - dpdx + nu*(d2udx2 + d2udy2));
    end
end
% (Similar loop for v_star)
end
```

### 4. Pressure Poisson Equation Solver

```matlab
function p_prime = solve_pressure_poisson(u_star, v_star, dx, dy, dt, tol, max_iter)
n = size(u_star,1);
p_prime = zeros(n);
for iter = 1:max_iter
    p_old = p_prime;
    for j = 2:n-1
        for i = 2:n-1
            rhs = ((u_star(j,i) - u_star(j,i-1))/dx + (v_star(j,i) - v_star(j-1,i))/dy) / dt;
            p_prime(j,i) = 0.25 * (p_prime(j,i+1) + p_prime(j,i-1) ...
                                   + p_prime(j+1,i) + p_prime(j-1,i) - dx^2 * rhs);
        end
    end
    % Neumann boundary conditions
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

### 5. Corrector Step (Velocity and Pressure Update)

```matlab
function [u, v, p] = corrector_step(u_star, v_star, p, p_prime, dx, dy, dt, alpha)
n = size(p,1);
u = u_star; v = v_star;
for j = 2:n-1
    for i = 2:n-1
        u(j,i) = u_star(j,i) - alpha * dt / dx * (p_prime(j,i+1) - p_prime(j,i));
        v(j,i) = v_star(j,i) - alpha * dt / dy * (p_prime(j+1,i) - p_prime(j,i));
    end
end
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
        for j = 1:length(scene_data.frames)
            [A, map] = rgb2ind(scene_data.frames(j).cdata, 256);
            if j == 1
                imwrite(A, map, scene_data.filename, 'gif', 'LoopCount', Inf, 'DelayTime', 0.1);
            else
                imwrite(A, map, scene_data.filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
            end
        end
    end
end
end
```
---

## Simulation Outputs

### 🌀 Velocity Vectors — Flow Field Animation
<p align="center">
  <img src="https://raw.githubusercontent.com/Kandil2001/Lid-Cavity-Evolution/main/assets/matlab/iterative_velocity_vectors.gif" alt="Velocity Vectors GIF"/>
</p>

### ⚡ Velocity Magnitude — Speed Contours
<p align="center">
  <img src="https://raw.githubusercontent.com/Kandil2001/Lid-Cavity-Evolution/main/assets/matlab/iterative_velocity_contour.gif" alt="Velocity Magnitude GIF"/>
</p>

### 🌊 Streamlines — Flow Paths
<p align="center">
  <img src="https://raw.githubusercontent.com/Kandil2001/Lid-Cavity-Evolution/main/assets/matlab/iterative_streamlines.gif" alt="Streamlines GIF"/>
</p>

### 📊 Pressure Distribution
<p align="center">
  <img src="https://raw.githubusercontent.com/Kandil2001/Lid-Cavity-Evolution/main/assets/matlab/iterative_pressure_contour.gif" alt="Pressure GIF"/>
</p>

_GIF/image files are saved in the script directory. Filenames are customizable in the code._

---

## Performance and Convergence

**This is for 51x51 grid size**
| Metric                | Typical Value                  |
|-----------------------|--------------------------------|
| Elapsed Time          | ~1126,03 seconds (~18,77 minutes) |
| Avg. Time per Step    | ~2,2521 seconds                |

**This is for 151x151 grid size**
| Metric                | Typical Value                  |
|-----------------------|--------------------------------|
| Elapsed Time          | ~395,43 seconds (~6,59 minutes) |
| Avg. Time per Step    | ~0,7909 seconds                |

### 📉 Residuals — Convergence History
<p align="center">
  <img src="https://raw.githubusercontent.com/Kandil2001/Lid-Cavity-Evolution/main/assets/matlab/iterative_residuals.gif" alt="Residuals GIF"/>
</p>

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
