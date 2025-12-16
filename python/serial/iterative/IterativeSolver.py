import numpy as np
import matplotlib.pyplot as plt
import imageio
import time


def IterativeSolver():
    Re = 100               # Reynolds number
    L = 1.0                # Cavity length
    n = 51                 # Grid size (n x n)
    dt = 0.002             # Time step
    total_time = 1.0       # Total simulated time
    alpha_u = 0.7          # Under-relaxation for velocity
    alpha_p = 0.3          # Under-relaxation for pressure
    tol = 1e-6             # SIMPLE inner iteration tolerance
    max_iter = 300         # Max SIMPLE inner iterations per step
    record_gif = True      # Record GIFs

    nu = 1.0 / Re
    dx = L / (n - 1)
    dy = dx  

    X, Y = np.meshgrid(np.linspace(0.0, L, n), np.linspace(0.0, L, n))

    u = np.zeros((n, n))
    v = np.zeros((n, n))
    p = np.zeros((n, n))

    u[-1, :] = 1.0

    max_steps = int(np.ceil(total_time / dt))
    res_u_arr = np.zeros(max_steps)
    res_v_arr = np.zeros(max_steps)
    res_p_arr = np.zeros(max_steps)

    gif_scenes = {}
    if record_gif:
        gif_scenes = {
            'velocity_vectors': {'frames': [], 'filename': 'iterative_velocity_vectors.gif'},
            'velocity_contour': {'frames': [], 'filename': 'iterative_velocity_contour.gif'},
            'pressure_contour': {'frames': [], 'filename': 'iterative_pressure_contour.gif'},
            'residuals':        {'frames': [], 'filename': 'iterative_residuals.gif'},
            'streamlines':      {'frames': [], 'filename': 'iterative_streamlines.gif'},
        }

    print("Starting SIMPLE Lid Driven Cavity Simulation...")
    print(f"Grid size: {n}x{n}, Re: {Re}")
    print(f"dt: {dt:.4g}, tol: {tol:g}, max_iter: {max_iter}")

    start_time = time.time()
    step = 0
    t = 0.0

    while t < total_time:
        step += 1

        for iter in range(1, max_iter + 1):
            u_prev = u.copy()
            v_prev = v.copy()
            p_prev = p.copy()

            u_star, v_star = predictor_step(u, v, p, dx, dy, dt, nu, alpha_u)

            p_prime = solve_pressure_poisson(u_star, v_star, dx, dy, dt, tol, max_iter)

            u, v, p = corrector_step(u_star, v_star, p, p_prime, dx, dy, dt, alpha_p)

            u[:, 0] = 0.0
            u[:, -1] = 0.0
            u[0, :] = 0.0
            u[-1, :] = 1.0

            v[:, 0] = 0.0
            v[:, -1] = 0.0
            v[0, :] = 0.0
            v[-1, :] = 0.0

            p[0, :]  = p[1, :]
            p[-1, :] = p[-2, :]
            p[:, 0]  = p[:, 1]
            p[:, -1] = p[:, -2]

            res_u = np.max(np.abs(u - u_prev))
            res_v = np.max(np.abs(v - v_prev))
            res_p = np.max(np.abs(p - p_prev))

            if max(res_u, res_v, res_p) < tol:
                break

        res_u_arr[step - 1] = res_u
        res_v_arr[step - 1] = res_v
        res_p_arr[step - 1] = res_p

        if record_gif:
            capture_frames(
                gif_scenes, X, Y, u, v, p,
                res_u_arr[:step], res_v_arr[:step], res_p_arr[:step],
                step, iter, t
            )

        t += dt

    if record_gif:
        print("Creating GIF files...")
        create_gifs(gif_scenes)

    elapsed = time.time() - start_time
    print("\nSimulation complete.")
    print(f"Elapsed real time: {elapsed:.2f} seconds ({elapsed/60.0:.2f} minutes).")
    print(f"Total time steps: {step}, Final time: {t:.4f} s")
    print(f"Average time per step: {elapsed/step:.4f} seconds")

    plot_final_results(X, Y, u, v, p,
                       res_u_arr[:step], res_v_arr[:step], res_p_arr[:step])


def predictor_step(u, v, p, dx, dy, dt, nu, alpha):
    n = u.shape[0]
    u_star = u.copy()
    v_star = v.copy()

    inv_4dx = 1.0 / (4.0 * dx)
    inv_4dy = 1.0 / (4.0 * dy)
    inv_dx_sq = 1.0 / (dx * dx)
    inv_dy_sq = 1.0 / (dy * dy)
    alpha_dt = alpha * dt

    for j in range(1, n - 1):
        for i in range(1, n - 1):
            du2dx = ((u[j, i] + u[j, i + 1])**2 -
                     (u[j, i - 1] + u[j, i])**2) * inv_4dx
            duvdy = ((v[j, i] + v[j, i + 1]) * (u[j, i] + u[j + 1, i]) -
                     (v[j - 1, i] + v[j - 1, i + 1]) * (u[j - 1, i] + u[j, i])) * inv_4dy
            d2udx2 = (u[j, i + 1] - 2.0 * u[j, i] + u[j, i - 1]) * inv_dx_sq
            d2udy2 = (u[j + 1, i] - 2.0 * u[j, i] + u[j - 1, i]) * inv_dy_sq
            dpdx = (p[j, i + 1] - p[j, i]) / dx

            u_star[j, i] = u[j, i] + alpha_dt * (
                -du2dx - duvdy - dpdx + nu * (d2udx2 + d2udy2)
            )

    for j in range(1, n - 1):
        for i in range(1, n - 1):
            dv2dy = ((v[j, i] + v[j + 1, i])**2 -
                     (v[j - 1, i] + v[j, i])**2) * inv_4dy
            duvdx = ((u[j + 1, i] + u[j, i]) * (v[j, i + 1] + v[j, i]) -
                     (u[j + 1, i - 1] + u[j, i - 1]) * (v[j, i] + v[j, i - 1])) * inv_4dx
            d2vdx2 = (v[j, i + 1] - 2.0 * v[j, i] + v[j, i - 1]) * inv_dx_sq
            d2vdy2 = (v[j + 1, i] - 2.0 * v[j, i] + v[j - 1, i]) * inv_dy_sq
            dpdy = (p[j + 1, i] - p[j, i]) / dy

            v_star[j, i] = v[j, i] + alpha_dt * (
                -duvdx - dv2dy - dpdy + nu * (d2vdx2 + d2vdy2)
            )

    return u_star, v_star

def solve_pressure_poisson(u_star, v_star, dx, dy, dt, tol, max_iter):
    n = u_star.shape[0]
    p_prime = np.zeros((n, n))

    inv_dx = 1.0 / dx
    inv_dy = 1.0 / dy
    dt_rhs_factor = 1.0 / dt
    laplacian_factor = 0.25  # for dx = dy

    for _ in range(max_iter):
        p_old = p_prime.copy()

        for j in range(1, n - 1):
            for i in range(1, n - 1):
                rhs = ((u_star[j, i] - u_star[j, i - 1]) * inv_dx +
                       (v_star[j, i] - v_star[j - 1, i]) * inv_dy) * dt_rhs_factor

                p_prime[j, i] = laplacian_factor * (
                    p_prime[j, i + 1] + p_prime[j, i - 1] +
                    p_prime[j + 1, i] + p_prime[j - 1, i] -
                    dx * dx * rhs
                )

        p_prime[0, :] = p_prime[1, :]
        p_prime[-1, :] = p_prime[-2, :]
        p_prime[:, 0] = p_prime[:, 1]
        p_prime[:, -1] = p_prime[:, -2]

        if np.max(np.abs(p_prime - p_old)) < tol:
            break

    return p_prime


def corrector_step(u_star, v_star, p, p_prime, dx, dy, dt, alpha):
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

def capture_frames(gif_scenes, X, Y, u, v, p,
                   res_u, res_v, res_p,
                   step, iter, time_s):

    fig, ax = plt.subplots()
    ax.quiver(X[::4, ::4], Y[::4, ::4],
              u[::4, ::4], v[::4, ::4], scale=None)
    ax.set_title(f"Velocity Vectors\nStep {step}, SIMPLE iter {iter}, Time = {time_s:.3f} s")
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_aspect('equal')
    ax.set_xlim(X.min(), X.max())
    ax.set_ylim(Y.min(), Y.max())
    ax.grid(True)
    gif_scenes['velocity_vectors']['frames'].append(get_frame(fig))
    plt.close(fig)

    velMag = np.sqrt(u**2 + v**2)
    fig, ax = plt.subplots()
    cs = ax.contourf(X, Y, velMag, 20)
    fig.colorbar(cs)
    ax.set_title(f"Velocity Magnitude\nStep {step}, SIMPLE iter {iter}, Time = {time_s:.3f} s")
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_aspect('equal')
    ax.set_xlim(X.min(), X.max())
    ax.set_ylim(Y.min(), Y.max())
    gif_scenes['velocity_contour']['frames'].append(get_frame(fig))
    plt.close(fig)

    fig, ax = plt.subplots()
    cs = ax.contourf(X, Y, p, 20)
    fig.colorbar(cs)
    ax.set_title(f"Pressure Field\nStep {step}, SIMPLE iter {iter}, Time = {time_s:.3f} s")
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_aspect('equal')
    ax.set_xlim(X.min(), X.max())
    ax.set_ylim(Y.min(), Y.max())
    gif_scenes['pressure_contour']['frames'].append(get_frame(fig))
    plt.close(fig)

    fig, ax = plt.subplots()
    startx = np.linspace(0.0, 1.0, 15)
    starty = np.linspace(0.0, 1.0, 15)
    sx, sy = np.meshgrid(startx, starty)
    start_points = np.vstack([sx.ravel(), sy.ravel()]).T
    ax.streamplot(X, Y, u, v, start_points=start_points, density=1.0)
    ax.set_title(f"Streamlines\nStep {step}, SIMPLE iter {iter}, Time = {time_s:.3f} s")
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_aspect('equal')
    ax.set_xlim(X.min(), X.max())
    ax.set_ylim(Y.min(), Y.max())
    ax.grid(True)
    gif_scenes['streamlines']['frames'].append(get_frame(fig))
    plt.close(fig)

    fig, ax = plt.subplots()
    steps = np.arange(1, len(res_u) + 1)
    ax.semilogy(steps, res_u, '-r', label='u-res')
    ax.semilogy(steps, res_v, '-g', label='v-res')
    ax.semilogy(steps, res_p, '-b', label='p-res')
    ax.set_xlabel('Time Step')
    ax.set_ylabel('Residual (log scale)')
    ax.set_title(f"Convergence History\nStep {step}, SIMPLE iter {iter}, Time = {time_s:.3f} s")
    ax.grid(True, which='both')
    ax.legend(loc='upper right')
    gif_scenes['residuals']['frames'].append(get_frame(fig))
    plt.close(fig)


def get_frame(fig):
    fig.canvas.draw()
    w, h = fig.canvas.get_width_height()
    buf = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8)
    return buf.reshape((h, w, 3))


def create_gifs(gif_scenes):
    for scene_name, scene_data in gif_scenes.items():
        frames = scene_data['frames']
        filename = scene_data['filename']
        if frames:
            imageio.mimsave(filename, frames, duration=0.1)
            print(f"Saved: {filename}")


def curl(X, Y, u, v):
    dy = Y[1, 0] - Y[0, 0]
    dx = X[0, 1] - X[0, 0]
    dudy, dudx = np.gradient(u, dy, dx)
    dvdy, dvdx = np.gradient(v, dy, dx)
    return dvdx - dudy


def plot_final_results(X, Y, u, v, p, res_u, res_v, res_p):
    n = u.shape[0]
    fig, axs = plt.subplots(2, 3, figsize=(15, 9))
    fig.subplots_adjust(wspace=0.3, hspace=0.3)

    ax = axs[0, 0]
    ax.quiver(X[::3, ::3], Y[::3, ::3],
              u[::3, ::3], v[::3, ::3], scale=None)
    startx = np.linspace(0.0, 1.0, 20)
    starty = np.linspace(0.0, 1.0, 20)
    sx, sy = np.meshgrid(startx, starty)
    start_points = np.vstack([sx.ravel(), sy.ravel()]).T
    ax.streamplot(X, Y, u, v, start_points=start_points, density=1.0)
    ax.set_title('Velocity Vectors and Streamlines')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_aspect('equal')
    ax.set_xlim(X.min(), X.max())
    ax.set_ylim(Y.min(), Y.max())
    ax.grid(True)

    ax = axs[0, 1]
    velMag = np.sqrt(u**2 + v**2)
    cs = ax.contourf(X, Y, velMag, 20)
    fig.colorbar(cs, ax=ax)
    ax.set_title('Velocity Magnitude')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_aspect('equal')
    ax.set_xlim(X.min(), X.max())
    ax.set_ylim(Y.min(), Y.max())

    ax = axs[0, 2]
    cs = ax.contourf(X, Y, p, 20)
    fig.colorbar(cs, ax=ax)
    ax.set_title('Pressure Field')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_aspect('equal')
    ax.set_xlim(X.min(), X.max())
    ax.set_ylim(Y.min(), Y.max())

    ax = axs[1, 0]
    vort = curl(X, Y, u, v)
    cs = ax.contourf(X, Y, vort, 20)
    fig.colorbar(cs, ax=ax)
    ax.set_title('Vorticity Field')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_aspect('equal')
    ax.set_xlim(X.min(), X.max())
    ax.set_ylim(Y.min(), Y.max())

    ax = axs[1, 1]
    mid = n // 2
    ax.plot(u[mid, :], Y[mid, :], 'b-', linewidth=2, label='Vertical Centerline')
    ax.plot(u[:, mid], X[:, mid], 'r-', linewidth=2, label='Horizontal Centerline')
    ax.set_title('Centerline Velocity Profiles')
    ax.set_xlabel('Velocity')
    ax.set_ylabel('Position')
    ax.legend()
    ax.grid(True)

    ax = axs[1, 2]
    steps = np.arange(1, len(res_u) + 1)
    ax.semilogy(steps, res_u, 'r-', label='u-residual')
    ax.semilogy(steps, res_v, 'g-', label='v-residual')
    ax.semilogy(steps, res_p, 'b-', label='p-residual')
    ax.set_title('Convergence History')
    ax.set_xlabel('Time Step')
    ax.set_ylabel('Residual (log scale)')
    ax.legend(loc='upper right')
    ax.grid(True, which='both')

    fig.text(0.02, 0.02,
             f"Re = 100, Grid = {n}x{n}",
             fontsize=9,
             bbox=dict(facecolor='white', edgecolor='black'))

    plt.tight_layout()
    plt.savefig('final_results.png', dpi=150)
    plt.show()


if __name__ == "__main__":
    IterativeSolver()
