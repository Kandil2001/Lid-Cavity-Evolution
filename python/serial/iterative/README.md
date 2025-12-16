<p align="center">
  <a href="https://www.python.org/">
    <img src="https://github.com/Kandil2001/Lid-Cavity-Evolution/raw/main/logos/python.png" width="70"/>
  </a>
</p>
<h1 align="center">🌀 IterativeSolver.py — Python Iterative SIMPLE Solver</h1>
<p align="center"><i>Educational loop-based SIMPLE algorithm for unsteady incompressible CFD</i></p>
<p align="center">
  <a href="../../../../LICENSE">
    <img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="MIT License"/>
  </a>
  <a href="https://github.com/Kandil2001/Lid-Cavity-Evolution/releases">
    <img src="https://img.shields.io/badge/Version-0.1.0-green.svg" alt="Version"/>
  </a>
  <a href="https://github.com/Kandil2001/Lid-Cavity-Evolution/blob/main/CONTRIBUTING.md">
    <img src="https://img.shields.io/badge/Contributions-Welcome-orange.svg" alt="Contributions Welcome"/>
  </a>
  <a href="https://www.python.org/">
    <img src="https://img.shields.io/badge/Python-3.8+-blue.svg" alt="Python"/>
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

This implementation contains an **iterative (loop-based) Python version** of the SIMPLE algorithm for the lid-driven cavity problem.  
It is designed for educational clarity and direct benchmarking against the vectorized version, using **identical structure, boundary conditions, and outputs** for fair comparison and algorithm understanding.

- **File:** `IterativeSolver.py`
- **Location:** `main/python/serial/iterative/`
- **Focus:** Educational clarity, algorithm transparency, and performance benchmarking

---

## Key Features

- Finite volume discretization on a staggered grid
- **Triple-nested loop implementation** for maximum educational clarity
- Modular function structure for maintainability and understanding
- Real-time visualization and optional GIF export
- Residual tracking, convergence monitoring, and performance reporting
- Identical interface and outputs to vectorized solver for direct benchmarking

---

## Why Iterative?

- **Educational clarity:** Explicit loops make the algorithm steps transparent and easy to follow
- **Debugging friendly:** Step-by-step execution allows for easy inspection of intermediate values
- **Baseline performance:** Provides reference timing for vectorized optimization comparisons
- **Algorithm understanding:** Helps students and researchers understand the fundamental SIMPLE steps
- **Foundation for learning:** Serves as a starting point for understanding CFD algorithm implementation

---

## Usage Instructions

1. **Requirements:** Python 3.8+ with NumPy, Matplotlib, and imageio
2. **Setup:**  
   Place `IterativeSolver.py` in `main/python/serial/iterative/`
3. **Configure Parameters:**  
   At the top of the file:
   ```python
   Re = 100          # Reynolds number
   n = 51            # Grid size
   dt = 0.002        # Time step
   total_time = 1.0  # Simulation time (s)
   alpha_u = 0.7     # Under-relaxation for velocity
   alpha_p = 0.3     # Under-relaxation for pressure
   tol = 1e-6        # SIMPLE inner iteration tolerance
   max_iter = 300    # Max SIMPLE iterations per time step
   record_gif = True # Enable GIF recording
    ```
4. **Run the Script:**  
   - In terminal: `python IterativeSolver.py`
   - The simulation will begin, displaying progress and performance metrics
   - Animated figures and GIFs for each scene are generated automatically
5. **Check Results:**  
   - After completion, find GIFs and a summary plot in your working directory
   - The summary plot includes velocity vectors, streamlines, velocity magnitude, pressure field, vorticity, centerline profiles, and convergence history

---

## Key Code Snippets

Below are essential excerpts from `IterativeSolver.py`, demonstrating the transparent, loop-based structure of this solver.

### 1. User-Adjustable Simulation Parameters

```python
# Set physical and numerical parameters
Re = 100          # Reynolds number
n = 51            # Grid size
dt = 0.002        # Time step
total_time = 1.0  # Simulation time
alpha_u = 0.7     # Under-relaxation (velocity)
alpha_p = 0.3     # Under-relaxation (pressure)
tol = 1e-6        # Convergence tolerance
max_iter = 300    # Max SIMPLE iterations per time step
record_gif = True # Enable GIF creation
```

### 2. Main Time-Stepping and SIMPLE Iteration Loop

```python
while sim_time < total_time:
    step += 1
    for iter_ in range(1, max_iter + 1):
        u_prev = u.copy()
        v_prev = v.copy()
        p_prev = p.copy()
        
        # SIMPLE algorithm steps (iterative loops)
        u_star, v_star = predictor_step_iter(u, v, p, dx, dy, dt, nu, alpha_u)
        p_prime = solve_pressure_poisson_iter(u_star, v_star, dx, dy, dt, tol, max_iter)
        u, v, p = corrector_step_iter(u_star, v_star, p, p_prime, dx, dy, dt, alpha_p)
        
        # Boundary conditions
        u[:, 0] = 0; u[:, -1] = 0; u[0, :] = 0; u[-1, :] = 1
        v[:, 0] = 0; v[:, -1] = 0; v[0, :] = 0; v[-1, :] = 0
        p[0, :] = p[1, :]; p[-1, :] = p[-2, :]
        p[:, 0] = p[:, 1]; p[:, -1] = p[:, -2]
        
        # Residuals/convergence
        res_u = np.max(np.abs(u - u_prev))
        res_v = np.max(np.abs(v - v_prev))
        res_p = np.max(np.abs(p - p_prev))
        if max(res_u, res_v, res_p) < tol:
            break

    # Visualization, GIF capture, etc.
    sim_time += dt
```

### 3. Iterative Predictor Step

```python
def predictor_step_iter(u, v, p, dx, dy, dt, nu, alpha):
    n = u.shape[0]
    u_star = u.copy()
    v_star = v.copy()
    inv_4dx = 1 / (4 * dx)
    inv_4dy = 1 / (4 * dy)
    inv_dx_sq = 1 / dx ** 2
    inv_dy_sq = 1 / dy ** 2
    alpha_dt = alpha * dt
    
    # u-momentum (explicit loops)
    for j in range(1, n - 1):
        for i in range(1, n - 1):
            du2dx = ((u[j, i] + u[j, i + 1]) ** 2 - (u[j, i - 1] + u[j, i]) ** 2) * inv_4dx
            duvdy = ((v[j, i] + v[j, i + 1]) * (u[j, i] + u[j + 1, i])
                   - (v[j - 1, i] + v[j - 1, i + 1]) * (u[j - 1, i] + u[j, i])) * inv_4dy
            d2udx2 = (u[j, i + 1] - 2 * u[j, i] + u[j, i - 1]) * inv_dx_sq
            d2udy2 = (u[j + 1, i] - 2 * u[j, i] + u[j - 1, i]) * inv_dy_sq
            dpdx = (p[j, i + 1] - p[j, i]) / dx
            u_star[j, i] = u[j, i] + alpha_dt * (-du2dx - duvdy - dpdx + nu * (d2udx2 + d2udy2))
    
    # v-momentum (explicit loops)
    for j in range(1, n - 1):
        for i in range(1, n - 1):
            dv2dy = ((v[j, i] + v[j + 1, i]) ** 2 - (v[j - 1, i] + v[j, i]) ** 2) * inv_4dy
            duvdx = ((u[j + 1, i] + u[j, i]) * (v[j, i + 1] + v[j, i])
                   - (u[j + 1, i - 1] + u[j, i - 1]) * (v[j, i] + v[j, i - 1])) * inv_4dx
            d2vdx2 = (v[j, i + 1] - 2 * v[j, i] + v[j, i - 1]) * inv_dx_sq
            d2vdy2 = (v[j + 1, i] - 2 * v[j, i] + v[j - 1, i]) * inv_dy_sq
            dpdy = (p[j + 1, i] - p[j, i]) / dy
            v_star[j, i] = v[j, i] + alpha_dt * (-duvdx - dv2dy - dpdy + nu * (d2vdx2 + d2vdy2))
    
    return u_star, v_star
```
### 4. Iterative Pressure Poisson Solver

```python
def solve_pressure_poisson_iter(u_star, v_star, dx, dy, dt, tol, max_iter):
    n = u_star.shape[0]
    p_prime = np.zeros((n, n))
    inv_dx = 1 / dx
    inv_dy = 1 / dy
    dt_rhs_factor = 1 / dt
    laplacian_factor = 0.25
    
    for it in range(max_iter):
        p_old = p_prime.copy()
        for j in range(1, n - 1):
            for i in range(1, n - 1):
                rhs = ((u_star[j, i] - u_star[j, i - 1]) * inv_dx +
                       (v_star[j, i] - v_star[j - 1, i]) * inv_dy) * dt_rhs_factor
                p_prime[j, i] = laplacian_factor * (
                    p_prime[j, i + 1] + p_prime[j, i - 1] +
                    p_prime[j + 1, i] + p_prime[j - 1, i] - dx ** 2 * rhs)
        
        # Homogeneous Neumann BCs
        p_prime[0, :] = p_prime[1, :]
        p_prime[-1, :] = p_prime[-2, :]
        p_prime[:, 0] = p_prime[:, 1]
        p_prime[:, -1] = p_prime[:, -2]
        
        if np.max(np.abs(p_prime - p_old)) < tol:
            break
            
    return p_prime
```

### 5. Iterative Corrector Step

```python
def corrector_step_iter(u_star, v_star, p, p_prime, dx, dy, dt, alpha):
    n = p.shape[0]
    u = u_star.copy()
    v = v_star.copy()
    alpha_dt_dx = alpha * dt / dx
    alpha_dt_dy = alpha * dt / dy
    
    for j in range(1, n - 1):
        for i in range(1, n - 1):
            u[j, i] = u_star[j, i] - alpha_dt_dx * (p_prime[j, i + 1] - p_prime[j, i])
            v[j, i] = v_star[j, i] - alpha_dt_dy * (p_prime[j + 1, i] - p_prime[j, i])
    
    p = p + alpha * p_prime
    return u, v, p
```

### 6. GIF Creation and Post-Processing
```python
def fig_to_image(fig):
    """Convert a Matplotlib figure to a numpy RGB image."""
    fig.canvas.draw()
    w, h = fig.canvas.get_width_height()
    img = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8).reshape(h, w, 3)
    return img

# GIF recording structure
if record_gif:
    gif_scenes = {
        "velocity_vectors": [],
        "velocity_contour": [],
        "pressure_contour": [],
        "residuals": [],
        "streamlines": []
    }
    gif_filenames = {
        "velocity_vectors": "iterative_velocity_vectors.gif",
        "velocity_contour": "iterative_velocity_contour.gif",
        "pressure_contour": "iterative_pressure_contour.gif",
        "residuals": "iterative_residuals.gif",
        "streamlines": "iterative_streamlines.gif"
    }

# Create GIFs from frames
if record_gif:
    print('Creating GIF files...')
    for key, frames in gif_scenes.items():
        filename = gif_filenames[key]
        imageio.mimsave(filename, frames, duration=0.1)
        print(f'Saved: {filename}')
```

---


# Simulation Outputs Section

## Simulation Outputs

### 🌀 Velocity Vectors — Flow Field Animation
<p align="center">
  <img src="https://github.com/Kandil2001/Lid-Cavity-Evolution/raw/main/assets/python/iterative_velocity_vectors.gif" alt="Velocity Vectors GIF"/>
</p>

### ⚡ Velocity Magnitude — Speed Contours
<p align="center">
  <img src="https://github.com/Kandil2001/Lid-Cavity-Evolution/raw/main/assets/python/iterative_velocity_contour.gif" alt="Velocity Magnitude GIF"/>
</p>

### 🌊 Streamlines — Flow Paths
<p align="center">
  <img src="https://github.com/Kandil2001/Lid-Cavity-Evolution/raw/main/assets/python/iterative_streamlines.gif" alt="Streamlines GIF"/>
</p>

### 📊 Pressure Distribution
<p align="center">
  <img src="https://github.com/Kandil2001/Lid-Cavity-Evolution/raw/main/assets/python/iterative_pressure_contour.gif" alt="Pressure GIF"/>
</p>

_GIF files are saved in the working directory with descriptive names_

---

## Performance and Convergence

| Metric                | Typical Value                  |
|-----------------------|-------------------------------|
| Elapsed Time          | ~86036.84 seconds (~1433.95 minutes)|
| Avg. Time per Step    | ~172.0737 seconds              |

### 📉 Residuals — Convergence History
<p align="center">
  <img src="https://github.com/Kandil2001/Lid-Cavity-Evolution/raw/main/assets/python/iterative_residuals.gif" alt="Iterative Residuals GIF"/>
</p>

---

## Directory Placement

This solver is part of the Python serial implementations:
```
main/
├── python/
│ ├── serial/
│ │ ├── iterative-solver/
│ │ │ ├── IterativeSolver.py
│ │ │ └── README.md
│ │ ├── vectorized-solver/
│ │ │ ├── VectorizedSolver.py
│ │ │ └── README.md
│ │ └── README.md
│ ├── parallel/
│ └── README.md
├── matlab/
│ ├── ...
└── assets/
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

**This implementation provides the Python iterative baseline for the Lid Cavity Evolution suite.  
For vectorized Python versions, MATLAB implementations, and parallel solvers, see the main project documentation.**
