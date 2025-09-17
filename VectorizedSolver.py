"""
VectorizedSolver.py

SIMPLE 2D Lid-Driven Cavity Solver (Finite Volume, Staggered Grid, Vectorized, Python Version)

Author: Your Name (your.email@domain.com)
License: MIT

Usage:
    - Requires: numpy, matplotlib, imageio
    - Run: python VectorizedSolver.py
    - Adjust parameters at the top of the file as needed
    - GIFs and summary plots are generated in the working directory

Features:
    - Fully vectorized, modular implementation for Python/NumPy
    - GIF recording: each variable/scene gets its own GIF, titles include time step and SIMPLE iteration
    - Final summary plot with velocity, pressure, vorticity, centerline profiles, residuals
    - Residual tracking and performance reporting
    - Structure and interface match iterative version for direct benchmarking
"""

import numpy as np
import matplotlib.pyplot as plt
import imageio
import os
from time import time

# USER-ADJUSTABLE PARAMETERS
Re = 100         # Reynolds number
L = 1.0          # Cavity length
n = 51           # Grid size
dt = 0.002       # Time step
total_time = 1.0 # Total simulated time
alpha_u = 0.7    # Under-relaxation for velocity
alpha_p = 0.3    # Under-relaxation for pressure
tol = 1e-6       # SIMPLE inner iteration tolerance
max_iter = 300   # Max SIMPLE inner iterations per step
record_gif = True # Record GIFs

def main():
    # INITIALIZATION
    nu = 1 / Re
    dx = L / (n - 1)
    dy = dx
    max_steps = int(np.ceil(total_time / dt))
    res_u_arr = np.zeros(max_steps)
    res_v_arr = np.zeros(max_steps)
    res_p_arr = np.zeros(max_steps)

    x = np.linspace(0, L, n)
    y = np.linspace(0, L, n)
    X, Y = np.meshgrid(x, y)

    u = np.zeros((n, n))
    v = np.zeros((n, n))
    p = np.zeros((n, n))
    u[-1, :] = 1.0  # Lid moves right

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

    print(f"Starting Vectorized SIMPLE Lid Driven Cavity Simulation...")
    print(f"Grid size: {n}x{n}, Re: {Re}")
    print(f"dt: {dt}, tol: {tol}, max_iter: {max_iter}")
    start_time = time()
    sim_time = 0.0
    step = 0

    # MAIN TIME STEPPING LOOP
    while sim_time < total_time:
        step += 1
        u_old = u.copy()
        v_old = v.copy()
        p_old = p.copy()

        # SIMPLE Inner Iterations (vectorized)
        for iter_ in range(1, max_iter + 1):
            u_prev = u.copy()
            v_prev = v.copy()
            p_prev = p.copy()
            u_star, v_star = predictor_step_vectorized(u, v, p, dx, dy, dt, nu, alpha_u)
            p_prime = solve_pressure_poisson_vectorized(u_star, v_star, dx, dy, dt, tol, max_iter)
            u, v, p = corrector_step_vectorized(u_star, v_star, p, p_prime, dx, dy, dt, alpha_p)

            # Boundary conditions (set after field update)
            u[:, 0] = 0; u[:, -1] = 0; u[0, :] = 0; u[-1, :] = 1
            v[:, 0] = 0; v[:, -1] = 0; v[0, :] = 0; v[-1, :] = 0
            p[0, :] = p[1, :]; p[-1, :] = p[-2, :]
            p[:, 0] = p[:, 1]; p[:, -1] = p[:, -2]

            # SIMPLE residuals
            res_u = np.max(np.abs(u - u_prev))
            res_v = np.max(np.abs(v - v_prev))
            res_p = np.max(np.abs(p - p_prev))
            if max(res_u, res_v, res_p) < tol:
                break

        res_u_arr[step - 1] = res_u
        res_v_arr[step - 1] = res_v
        res_p_arr[step - 1] = res_p

        # GIF CAPTURE: Each scene is a separate figure with iteration info
        if record_gif:
            # Velocity vectors
            fig1, ax1 = plt.subplots(figsize=(5, 5))
            ax1.quiver(X[::4, ::4], Y[::4, ::4], u[::4, ::4], v[::4, ::4], scale=2, color='k')
            ax1.set_title(f"Velocity Vectors\nStep {step}, SIMPLE iter {iter_}, Time = {sim_time:.3f} s")
            ax1.set_xlabel('X'); ax1.set_ylabel('Y'); ax1.axis('equal')
            ax1.set_xlim([0, L]); ax1.set_ylim([0, L])
            ax1.grid(True)
            fig1.tight_layout()
            gif_scenes["velocity_vectors"].append(fig_to_image(fig1))
            plt.close(fig1)

            # Velocity magnitude contour
            fig2, ax2 = plt.subplots(figsize=(5, 5))
            velMag = np.sqrt(u ** 2 + v ** 2)
            c = ax2.contourf(X, Y, velMag, 20)
            plt.colorbar(c, ax=ax2)
            ax2.set_title(f"Velocity Magnitude\nStep {step}, SIMPLE iter {iter_}, Time = {sim_time:.3f} s")
            ax2.set_xlabel('X'); ax2.set_ylabel('Y'); ax2.axis('equal')
            ax2.set_xlim([0, L]); ax2.set_ylim([0, L])
            fig2.tight_layout()
            gif_scenes["velocity_contour"].append(fig_to_image(fig2))
            plt.close(fig2)

            # Pressure contour
            fig3, ax3 = plt.subplots(figsize=(5, 5))
            c = ax3.contourf(X, Y, p, 20)
            plt.colorbar(c, ax=ax3)
            ax3.set_title(f"Pressure Field\nStep {step}, SIMPLE iter {iter_}, Time = {sim_time:.3f} s")
            ax3.set_xlabel('X'); ax3.set_ylabel('Y'); ax3.axis('equal')
            ax3.set_xlim([0, L]); ax3.set_ylim([0, L])
            fig3.tight_layout()
            gif_scenes["pressure_contour"].append(fig_to_image(fig3))
            plt.close(fig3)

            # Streamlines
            fig4, ax4 = plt.subplots(figsize=(5, 5))
            startx = np.linspace(0, L, 15)
            starty = np.linspace(0, L, 15)
            sx, sy = np.meshgrid(startx, starty)
            ax4.streamplot(X, Y, u, v, start_points=np.array([sx.flatten(), sy.flatten()]).T, color='b')
            ax4.set_title(f"Streamlines\nStep {step}, SIMPLE iter {iter_}, Time = {sim_time:.3f} s")
            ax4.set_xlabel('X'); ax4.set_ylabel('Y'); ax4.axis('equal')
            ax4.set_xlim([0, L]); ax4.set_ylim([0, L])
            fig4.tight_layout()
            gif_scenes["streamlines"].append(fig_to_image(fig4))
            plt.close(fig4)

            # Residuals
            fig5, ax5 = plt.subplots(figsize=(5, 5))
            ax5.semilogy(range(1, step + 1), res_u_arr[:step], 'r-', label='u-res')
            ax5.semilogy(range(1, step + 1), res_v_arr[:step], 'g-', label='v-res')
            ax5.semilogy(range(1, step + 1), res_p_arr[:step], 'b-', label='p-res')
            ax5.set_xlabel('Time Step')
            ax5.set_ylabel('Residual (log scale)')
            ax5.set_title(f"Convergence History\nStep {step}, SIMPLE iter {iter_}, Time = {sim_time:.3f} s")
            ax5.legend(loc='upper right')
            ax5.grid(True)
            fig5.tight_layout()
            gif_scenes["residuals"].append(fig_to_image(fig5))
            plt.close(fig5)

        # Advance time
        sim_time += dt

    # CREATE GIFs FROM FRAMES
    if record_gif:
        print('Creating GIF files...')
        for key, frames in gif_scenes.items():
            filename = gif_filenames[key]
            imageio.mimsave(filename, frames, duration=0.1)
            print(f'Saved: {filename}')

    # FINAL SUMMARY AND PLOTS
    elapsed_time = time() - start_time
    print('\nSimulation complete.')
    print(f'Elapsed real time: {elapsed_time:.2f} seconds ({elapsed_time / 60:.2f} minutes).')
    print(f'Total time steps: {step}, Final time: {sim_time:.4f} s')
    print(f'Average time per step: {elapsed_time / step:.4f} seconds')

    plot_final_results(X, Y, u, v, p, res_u_arr[:step], res_v_arr[:step], res_p_arr[:step])

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

def solve_pressure_poisson_vectorized(u_star, v_star, dx, dy, dt, tol, max_iter):
    n = u_star.shape[0]
    p_prime = np.zeros((n, n))
    inv_dx = 1 / dx
    inv_dy = 1 / dy
    dt_rhs_factor = 1 / dt
    laplacian_factor = 0.25
    jj, ii = np.meshgrid(np.arange(1, n-1), np.arange(1, n-1), indexing='ij')
    for _ in range(max_iter):
        p_old = p_prime.copy()
        rhs = ((u_star[jj, ii] - u_star[jj, ii-1]) * inv_dx +
               (v_star[jj, ii] - v_star[jj-1, ii]) * inv_dy) * dt_rhs_factor
        p_prime[jj, ii] = laplacian_factor * (
            p_prime[jj, ii+1] + p_prime[jj, ii-1] +
            p_prime[jj+1, ii] + p_prime[jj-1, ii] - dx ** 2 * rhs)
        # Homogeneous Neumann BCs
        p_prime[0, :] = p_prime[1, :]
        p_prime[-1, :] = p_prime[-2, :]
        p_prime[:, 0] = p_prime[:, 1]
        p_prime[:, -1] = p_prime[:, -2]
        if np.max(np.abs(p_prime - p_old)) < tol:
            break
    return p_prime

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

def plot_final_results(X, Y, u, v, p, res_u, res_v, res_p):
    fig = plt.figure(figsize=(16, 10))
    # Velocity vectors and streamlines
    ax1 = fig.add_subplot(2, 3, 1)
    ax1.quiver(X[::3, ::3], Y[::3, ::3], u[::3, ::3], v[::3, ::3], scale=1.5, color='k')
    startx = np.linspace(0, 1, 20)
    starty = np.linspace(0, 1, 20)
    sx, sy = np.meshgrid(startx, starty)
    ax1.streamplot(X, Y, u, v, start_points=np.array([sx.flatten(), sy.flatten()]).T, color='b')
    ax1.set_title('Velocity Vectors and Streamlines')
    ax1.set_xlabel('X'); ax1.set_ylabel('Y')
    ax1.axis('equal'); ax1.grid(True)
    # Velocity magnitude contour
    ax2 = fig.add_subplot(2, 3, 2)
    velMag = np.sqrt(u ** 2 + v ** 2)
    c2 = ax2.contourf(X, Y, velMag, 20)
    fig.colorbar(c2, ax=ax2)
    ax2.set_title('Velocity Magnitude')
    ax2.set_xlabel('X'); ax2.set_ylabel('Y'); ax2.axis('equal')
    # Pressure contour
    ax3 = fig.add_subplot(2, 3, 3)
    c3 = ax3.contourf(X, Y, p, 20)
    fig.colorbar(c3, ax=ax3)
    ax3.set_title('Pressure Field')
    ax3.set_xlabel('X'); ax3.set_ylabel('Y'); ax3.axis('equal')
    # Vorticity contour
    ax4 = fig.add_subplot(2, 3, 4)
    vorticity = compute_vorticity(X, Y, u, v)
    c4 = ax4.contourf(X, Y, vorticity, 20)
    fig.colorbar(c4, ax=ax4)
    ax4.set_title('Vorticity Field')
    ax4.set_xlabel('X'); ax4.set_ylabel('Y'); ax4.axis('equal')
    # Centerline velocity profiles
    ax5 = fig.add_subplot(2, 3, 5)
    ax5.plot(u[n // 2, :], Y[n // 2, :], 'b-', linewidth=2, label='Vertical Centerline')
    ax5.plot(u[:, n // 2], X[:, n // 2], 'r-', linewidth=2, label='Horizontal Centerline')
    ax5.set_title('Centerline Velocity Profiles')
    ax5.set_xlabel('Velocity'); ax5.set_ylabel('Position')
    ax5.legend(); ax5.grid(True)
    # Residual history
    ax6 = fig.add_subplot(2, 3, 6)
    ax6.semilogy(range(1, len(res_u) + 1), res_u, 'r-', label='u-residual')
    ax6.semilogy(range(1, len(res_v) + 1), res_v, 'g-', label='v-residual')
    ax6.semilogy(range(1, len(res_p) + 1), res_p, 'b-', label='p-residual')
    ax6.set_title('Convergence History')
    ax6.set_xlabel('Time Step'); ax6.set_ylabel('Residual (log scale)')
    ax6.legend(); ax6.grid(True)
    plt.tight_layout()
    plt.savefig('final_results.png')
    plt.show()

def compute_vorticity(X, Y, u, v):
    dudy, dudx = np.gradient(u, Y[0, 1] - Y[0, 0], X[1, 0] - X[0, 0])
    dvdy, dvdx = np.gradient(v, Y[0, 1] - Y[0, 0], X[1, 0] - X[0, 0])
    return dvdx - dudy

def fig_to_image(fig):
    """Convert a Matplotlib figure to a numpy RGB image."""
    fig.canvas.draw()
    w, h = fig.canvas.get_width_height()
    img = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8).reshape(h, w, 3)
    return img

if __name__ == "__main__":
    main()
