<p align="center">
  <a href="https://www.python.org/">
    <img src="https://github.com/Kandil2001/Lid-Cavity-Evolution/raw/main/logos/python.png" width="70"/>
  </a>
</p>
<h1 align="center">🌀 VectorizedSolver.py — Python Vectorized SIMPLE Solver</h1>
<p align="center"><i>High-performance vectorized SIMPLE algorithm for unsteady incompressible CFD</i></p>
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

This implementation contains a **fully vectorized Python/NumPy version** of the SIMPLE algorithm for the lid-driven cavity problem.  
It is designed for direct benchmarking against both the Python iterative version and MATLAB implementations, using **identical structure, boundary conditions, and outputs** for fair comparison and maximum performance.

- **File:** `VectorizedSolver.py`
- **Location:** `main/python/serial/vectorized/`
- **Focus:** Vectorized speed, educational clarity, and cross-language benchmarking

---

## Key Features

- Finite volume discretization on a staggered grid
- **Fully vectorized NumPy operations** for superior Python performance
- Modular function structure for maintainability and clarity
- Real-time visualization and optional GIF export
- Residual tracking, convergence monitoring, and performance reporting
- Identical interface and outputs to iterative solver for direct benchmarking

---

## Why Vectorized?

- **Maximum speed in Python:** Exploits NumPy's optimization for array operations, drastically reducing runtime for large grids
- **Direct benchmarking:** Uses the same interface and output as the loop version for accurate speed comparison
- **Educational value:** Demonstrates how to efficiently map CFD algorithms to Python's scientific computing ecosystem
- **Foundation for optimization:** Serves as baseline for further optimizations (Numba, Cython, etc.)

---

## Usage Instructions

1. **Requirements:** Python 3.8+ with NumPy, Matplotlib, and imageio
2. **Setup:**  
   Place `VectorizedSolver.py` in `main/python/serial/vectorized/`
3. **Configure Parameters:**  
   At the top of the file:
   ```python
   Re = 100         # Reynolds number
   n = 51           # Grid size
   dt = 0.002       # Time step
   total_time = 1.0 # Simulation time (s)
   alpha_u = 0.7    # Under-relaxation for velocity
   alpha_p = 0.3    # Under-relaxation for pressure
   tol = 1e-6       # SIMPLE inner iteration tolerance
   max_iter = 300   # Max SIMPLE iterations per time step
   record_gif = True # Enable GIF recording
   ```
4. **Run the Script:**  
   - In terminal: `python VectorizedSolver.py`
   - The simulation will begin, displaying progress and performance metrics
   - Animated figures and GIFs for each scene are generated automatically
5. **Check Results:**  
   - After completion, find GIFs and a summary plot in your working directory
   - The summary plot includes velocity vectors, streamlines, velocity magnitude, pressure field, vorticity, centerline profiles, and convergence history

---

## Key Code Snippets

Below are essential excerpts from `VectorizedSolver.py`, demonstrating the modular and transparent structure of this solver.

### 1. User-Adjustable Simulation Parameters

```python
# Set physical and numerical parameters
Re = 100         # Reynolds number
n = 51           # Grid size
dt = 0.002       # Time step
total_time = 1.0 # Simulation time
alpha_u = 0.7    # Under-relaxation (velocity)
alpha_p = 0.3    # Under-relaxation (pressure)
tol = 1e-6       # Convergence tolerance
max_iter = 300   # Max SIMPLE iterations per time step
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
        
        # SIMPLE algorithm steps
        u_star, v_star = predictor_step_vectorized(u, v, p, dx, dy, dt, nu, alpha_u)
        p_prime = solve_pressure_poisson_vectorized(u_star, v_star, dx, dy, dt, tol, max_iter)
        u, v, p = corrector_step_vectorized(u_star, v_star, p, p_prime, dx, dy, dt, alpha_p)
        
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

### 3. Vectorized Predictor Step

```python
def predictor_step_vectorized(u, v, p, dx, dy, dt, nu, alpha):
    n = u.shape[0]
    u_star = u.copy()
    v_star = v.copy()
    inv_4dx = 1 / (4 * dx)
    inv_4dy = 1 / (4 * dy)
    inv_dx_sq = 1 / dx ** 2
    inv_dy_sq = 1 / dy ** 2
    alpha_dt = alpha * dt

    jj, ii = np.meshgrid(np.arange(1, n-1), np.arange(1, n-1), indexing='ij')

    # u-momentum vectorized
    du2dx = ((u[jj, ii] + u[jj, ii+1]) ** 2 - (u[jj, ii-1] + u[jj, ii]) ** 2) * inv_4dx
    duvdy = ((v[jj, ii] + v[jj, ii+1]) * (u[jj, ii] + u[jj+1, ii])
           - (v[jj-1, ii] + v[jj-1, ii+1]) * (u[jj-1, ii] + u[jj, ii])) * inv_4dy
    d2udx2 = (u[jj, ii+1] - 2 * u[jj, ii] + u[jj, ii-1]) * inv_dx_sq
    d2udy2 = (u[jj+1, ii] - 2 * u[jj, ii] + u[jj-1, ii]) * inv_dy_sq
    dpdx = (p[jj, ii+1] - p[jj, ii]) / dx
    u_star[jj, ii] = u[jj, ii] + alpha_dt * (-du2dx - duvdy - dpdx + nu * (d2udx2 + d2udy2))

    # v-momentum vectorized
    dv2dy = ((v[jj, ii] + v[jj+1, ii]) ** 2 - (v[jj-1, ii] + v[jj, ii]) ** 2) * inv_4dy
    duvdx = ((u[jj+1, ii] + u[jj, ii]) * (v[jj, ii+1] + v[jj, ii])
           - (u[jj+1, ii-1] + u[jj, ii-1]) * (v[jj, ii] + v[jj, ii-1])) * inv_4dx
    d2vdx2 = (v[jj, ii+1] - 2 * v[jj, ii] + v[jj, ii-1]) * inv_dx_sq
    d2vdy2 = (v[jj+1, ii] - 2 * v[jj, ii] + v[jj-1, ii]) * inv_dy_sq
    dpdy = (p[jj+1, ii] - p[jj, ii]) / dy
    v_star[jj, ii] = v[jj, ii] + alpha_dt * (-duvdx - dv2dy - dpdy + nu * (d2vdx2 + d2vdy2))

    return u_star, v_star
```

### 4. Vectorized Pressure Poisson Solver

```python
def solve_pressure_poisson_vectorized(u_star, v_star, dx, dy, dt, tol, max_iter):
    n = u_star.shape[0]
    p_prime = np.zeros((n, n))
    
    jj, ii = np.meshgrid(np.arange(1, n-1), np.arange(1, n-1), indexing='ij')
    
    for _ in range(max_iter):
        p_old = p_prime.copy()
        
        rhs = ((u_star[jj, ii] - u_star[jj, ii-1]) / dx + 
               (v_star[jj, ii] - v_star[jj-1, ii]) / dy) / dt
        
        p_prime[jj, ii] = 0.25 * (p_prime[jj, ii+1] + p_prime[jj, ii-1] +
                                 p_prime[jj+1, ii] + p_prime[jj-1, ii] - dx**2 * rhs)
        
        # Neumann BCs
        p_prime[0, :] = p_prime[1, :]
        p_prime[-1, :] = p_prime[-2, :]
        p_prime[:, 0] = p_prime[:, 1]
        p_prime[:, -1] = p_prime[:, -2]
        
        if np.max(np.abs(p_prime - p_old)) < tol:
            break
            
    return p_prime
```

### 5. GIF Creation and Post-Processing

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
        "velocity_vectors": "vectorized_velocity_vectors.gif",
        "velocity_contour": "vectorized_velocity_contour.gif",
        "pressure_contour": "vectorized_pressure_contour.gif",
        "residuals": "vectorized_residuals.gif",
        "streamlines": "vectorized_streamlines.gif"
    }

# Create GIFs from frames
if record_gif:
    print('Creating GIF files...')
    for key, frames in gif_scenes.items():
        filename = gif_filenames[key]
        imageio.mimsave(filename, frames, duration=0.1)
        print(f'Saved: {filename}')
```

### 6. Vectorized Corrector Step

```python
def corrector_step_vectorized(u_star, v_star, p, p_prime, dx, dy, dt, alpha):
    n = p.shape[0]
    alpha_dt_dx = alpha * dt / dx
    alpha_dt_dy = alpha * dt / dy
    u = u_star.copy()
    v = v_star.copy()
    jj, ii = np.meshgrid(np.arange(1, n-1), np.arange(1, n-1), indexing='ij')
    u[jj, ii] = u_star[jj, ii] - alpha_dt_dx * (p_prime[jj, ii+1] - p_prime[jj, ii])
    v[jj, ii] = v_star[jj, ii] - alpha_dt_dy * (p_prime[jj+1, ii] - p_prime[jj, ii])
    p = p + alpha * p_prime
    return u, v, p
```

---

## Simulation Outputs

### 🌀 Velocity Vectors — Flow Field Animation
<p align="center">
  <img src="https://github.com/Kandil2001/Lid-Cavity-Evolution/raw/main/assets/python/vectorized_velocity_vectors.gif" alt="Velocity Vectors GIF"/>
</p>

### ⚡ Velocity Magnitude — Speed Contours
<p align="center">
  <img src="https://github.com/Kandil2001/Lid-Cavity-Evolution/raw/main/assets/python/vectorized_velocity_contour.gif" alt="Velocity Magnitude GIF"/>
</p>

### 🌊 Streamlines — Flow Paths
<p align="center">
  <img src="https://github.com/Kandil2001/Lid-Cavity-Evolution/raw/main/assets/python/vectorized_streamlines.gif" alt="Streamlines GIF"/>
</p>

### 📊 Pressure Distribution
<p align="center">
  <img src="https://github.com/Kandil2001/Lid-Cavity-Evolution/raw/main/assets/python/vectorized_pressure_contour.gif" alt="Pressure GIF"/>
</p>

_GIF files are saved in the working directory with descriptive names_

---

## Performance and Convergence

| Metric                | Typical Value                  |
|-----------------------|-------------------------------|
| Elapsed Time          | ~250-500 seconds (~4-8 minutes)|
| Avg. Time per Step    | ~0.5-1.0 seconds              |
| Speedup vs Iterative  | 5-20x faster                  |

### 📉 Residuals — Convergence History
<p align="center">
  <img src="https://github.com/Kandil2001/Lid-Cavity-Evolution/raw/main/assets/python/vectorized_residuals.gif" alt="Vectorized Residuals GIF"/>
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
See [main README](../../../../README.md) for project-wide context

---

## Contributing

Contributions, suggestions, and improvements are encouraged!  
Refer to [CONTRIBUTING.md](../../../../CONTRIBUTING.md)

---

## License

Released under the MIT License.  
See [LICENSE](../../../../LICENSE)

---

## References

1. Ghia, U., Ghia, K. N., & Shin, C. T. (1982). J. Comput. Phys., 48(3), 387-411
2. Patankar, S. V. (1980). _Numerical Heat Transfer and Fluid Flow_
3. Ferziger, J. H., Perić, M., & Street, R. L. (2002). _Computational Methods for Fluid Dynamics_

---

## Contact

- [GitHub Issues](https://github.com/Kandil2001/Lid-Cavity-Evolution/issues)
- Email: **kandil.ahmed.amr@gmail.com**
- [LinkedIn](https://www.linkedin.com/in/ahmed-kandil01)

---

**This implementation provides the Python vectorized baseline for the Lid Cavity Evolution suite.  
For iterative Python versions, MATLAB implementations, and parallel solvers, see the main project documentation.**
