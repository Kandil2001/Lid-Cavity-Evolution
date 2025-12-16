import numpy as np
import matplotlib.pyplot as plt
import imageio
import time

def VectorizedSolver():
    # USER PARAMETERS
    Re = 100
    L = 1.0
    n = 51
    dt = 0.002
    total_time = 1.0
    alpha_u = 0.7
    alpha_p = 0.3
    tol = 1e-6
    max_iter = 300
    record_gif = True

    # INITIALIZATION
    nu = 1.0 / Re
    dx = L / (n - 1)
    dy = dx
    max_steps = int(np.ceil(total_time / dt))

    res_u_arr = np.zeros(max_steps)
    res_v_arr = np.zeros(max_steps)
    res_p_arr = np.zeros(max_steps)

    # Grid generation
    X, Y = np.meshgrid(np.linspace(0, L, n), np.linspace(0, L, n))
    
    # Fields
    u = np.zeros((n, n))
    v = np.zeros((n, n))
    p = np.zeros((n, n))

    # Lid boundary condition
    u[-1, :] = 1.0

    gif_scenes = {}
    if record_gif:
        gif_scenes = {
            'velocity_vectors': [],
            'velocity_contour': [],
            'pressure_contour': [],
            'streamlines': [],
            'residuals': []
        }

    print(f"Starting Vectorized SIMPLE Lid Driven Cavity Simulation...\nGrid size: {n}x{n}, Re: {Re}")
    start_time = time.time()
    step = 0
    sim_time = 0.0

    while sim_time < total_time:
        step += 1
        u_old, v_old, p_old = u.copy(), v.copy(), p.copy()

        # SIMPLE inner iterations
        for iter in range(1, max_iter + 1):
            u_prev, v_prev, p_prev = u.copy(), v.copy(), p.copy()

            # 1. Predictor step (Momentum)
            u_star, v_star = predictor_step_vectorized(u, v, p, dx, dy, dt, nu, alpha_u)
            
            # 2. Pressure Poisson Equation
            p_prime = solve_pressure_poisson_vectorized(u_star, v_star, dx, dy, dt, tol, max_iter)
            
            # 3. Corrector step
            u, v, p = corrector_step_vectorized(u_star, v_star, p, p_prime, dx, dy, dt, alpha_p)

            # Apply Boundary conditions
            # No-slip walls
            u[:, 0] = 0; u[:, -1] = 0; u[0, :] = 0
            u[-1, :] = 1.0 # Lid

            v[:, 0] = 0; v[:, -1] = 0
            v[0, :] = 0; v[-1, :] = 0

            # Pressure Boundary Conditions (Neumann: dp/dn = 0)
            p[0, :] = p[1, :]
            p[-1, :] = p[-2, :]
            p[:, 0] = p[:, 1]
            p[:, -1] = p[:, -2]

            # Calculate residuals
            res_u = np.max(np.abs(u - u_prev))
            res_v = np.max(np.abs(v - v_prev))
            res_p = np.max(np.abs(p - p_prev))

            if max(res_u, res_v, res_p) < tol:
                break

        res_u_arr[step-1] = res_u
        res_v_arr[step-1] = res_v
        res_p_arr[step-1] = res_p

        # Capture frames for GIF
        if record_gif:
            capture_frames(gif_scenes, X, Y, u, v, p, res_u_arr[:step], res_v_arr[:step], res_p_arr[:step], step, iter, sim_time)

        sim_time += dt

    if record_gif:
        print("Creating GIF files...")
        create_gifs(gif_scenes)

    elapsed = time.time() - start_time
    print(f"\nSimulation complete.\nElapsed: {elapsed:.2f} s ({elapsed/60:.2f} min)")
    print(f"Total steps: {step}, Final time: {sim_time:.4f} s, Avg. step time: {elapsed/step:.4f} s")

    plot_final_results(X, Y, u, v, p, res_u_arr[:step], res_v_arr[:step], res_p_arr[:step])


def predictor_step_vectorized(u, v, p, dx, dy, dt, nu, alpha):
    inv_4dx = 1.0 / (4.0 * dx)
    inv_4dy = 1.0 / (4.0 * dy)
    inv_dx_sq = 1.0 / (dx * dx)
    inv_dy_sq = 1.0 / (dy * dy)
    alpha_dt = alpha * dt

    u_star = u.copy()
    v_star = v.copy()

    ui = u[1:-1, 1:-1]
    vi = v[1:-1, 1:-1]
    pi = p[1:-1, 1:-1]

    # convection in u
    du2dx = ((ui + u[1:-1, 2:])**2 - (u[1:-1, 0:-2] + ui)**2) * inv_4dx
    duvdy = ((v[1:-1, 1:-1] + v[1:-1, 2:]) * (ui + u[2:, 1:-1]) -
             (v[0:-2, 1:-1] + v[0:-2, 2:]) * (u[0:-2, 1:-1] + ui)) * inv_4dy

    d2udx2 = (u[1:-1, 2:] - 2.0 * ui + u[1:-1, 0:-2]) * inv_dx_sq
    d2udy2 = (u[2:, 1:-1] - 2.0 * ui + u[0:-2, 1:-1]) * inv_dy_sq

    dpdx = (p[1:-1, 2:] - pi) / dx

    u_star[1:-1, 1:-1] = ui + alpha_dt * (
        -du2dx - duvdy - dpdx + nu * (d2udx2 + d2udy2)
    )

    # convection in v
    dv2dy = ((vi + v[2:, 1:-1])**2 - (v[0:-2, 1:-1] + vi)**2) * inv_4dy
    duvdx = ((u[2:, 1:-1] + ui) * (v[1:-1, 2:] + vi) -
             (u[2:, 0:-2] + u[1:-1, 0:-2]) * (vi + v[1:-1, 0:-2])) * inv_4dx

    d2vdx2 = (v[1:-1, 2:] - 2.0 * vi + v[1:-1, 0:-2]) * inv_dx_sq
    d2vdy2 = (v[2:, 1:-1] - 2.0 * vi + v[0:-2, 1:-1]) * inv_dy_sq

    dpdy = (p[2:, 1:-1] - pi) / dy

    v_star[1:-1, 1:-1] = vi + alpha_dt * (
        -duvdx - dv2dy - dpdy + nu * (d2vdx2 + d2vdy2)
    )

    return u_star, v_star


def solve_pressure_poisson_vectorized(u_star, v_star, dx, dy, dt, tol, max_iter):
    n = u_star.shape[0]
    p_prime = np.zeros((n, n))
    inv_dx = 1.0 / dx
    inv_dy = 1.0 / dy
    dt_rhs_factor = 1.0 / dt
    laplacian_factor = 0.25

    # Pre-calculate RHS using explicit slices instead of slice objects
    # RHS = (du/dx + dv/dy) / dt
    # u_star diff: center - left (backward difference approximation for divergence here)
    div_u = (u_star[1:-1, 1:-1] - u_star[1:-1, 0:-2]) * inv_dx
    div_v = (v_star[1:-1, 1:-1] - v_star[0:-2, 1:-1]) * inv_dy
    rhs = (div_u + div_v) * dt_rhs_factor
    
    rhs_term = dx**2 * rhs

    for _ in range(max_iter):
        p_old = p_prime.copy()
        
        # Poisson update: p[i,j] = 0.25 * (p[i,j+1] + p[i,j-1] + p[i+1,j] + p[i-1,j] - dx^2 * RHS)
        p_prime[1:-1, 1:-1] = laplacian_factor * (
            p_prime[1:-1, 2:] + p_prime[1:-1, 0:-2] +
            p_prime[2:, 1:-1] + p_prime[0:-2, 1:-1] - rhs_term
        )

        # Neumann BC
        p_prime[0, :] = p_prime[1, :]
        p_prime[-1, :] = p_prime[-2, :]
        p_prime[:, 0] = p_prime[:, 1]
        p_prime[:, -1] = p_prime[:, -2]

        if np.max(np.abs(p_prime - p_old)) < tol:
            break

    return p_prime


def corrector_step_vectorized(u_star, v_star, p, p_prime, dx, dy, dt, alpha):
    alpha_dt_dx = alpha * dt / dx
    alpha_dt_dy = alpha * dt / dy

    u = u_star.copy()
    v = v_star.copy()

    # Velocity correction: u = u* - alpha*dt*dp'/dx
    # Forward difference for pressure grad
    u[1:-1, 1:-1] = u_star[1:-1, 1:-1] - alpha_dt_dx * (p_prime[1:-1, 2:] - p_prime[1:-1, 1:-1])
    
    # v = v* - alpha*dt*dp'/dy
    v[1:-1, 1:-1] = v_star[1:-1, 1:-1] - alpha_dt_dy * (p_prime[2:, 1:-1] - p_prime[1:-1, 1:-1])
    
    # Pressure update
    p = p + alpha * p_prime

    return u, v, p


def capture_frames(gif_scenes, X, Y, u, v, p, res_u, res_v, res_p, step, iter, time_s):
    # Only capture every 5 steps to save time/memory, unless it's very early
    if step % 5 != 0 and step > 10:
        return

    fig, ax = plt.subplots()
    # Subsample vectors for cleaner plot
    stride = 4
    ax.quiver(X[::stride, ::stride], Y[::stride, ::stride], u[::stride, ::stride], v[::stride, ::stride])
    ax.set_title(f"Velocity Vectors\nStep {step}, Iter {iter}, Time={time_s:.3f}s")
    gif_scenes['velocity_vectors'].append(get_frame(fig))
    plt.close(fig)

    fig, ax = plt.subplots()
    velMag = np.sqrt(u**2 + v**2)
    cs = ax.contourf(X, Y, velMag, 20, cmap='viridis')
    fig.colorbar(cs)
    ax.set_title(f"Velocity Magnitude\nStep {step}, Iter {iter}, Time={time_s:.3f}s")
    gif_scenes['velocity_contour'].append(get_frame(fig))
    plt.close(fig)

    fig, ax = plt.subplots()
    cs = ax.contourf(X, Y, p, 20, cmap='coolwarm')
    fig.colorbar(cs)
    ax.set_title(f"Pressure Field\nStep {step}, Iter {iter}, Time={time_s:.3f}s")
    gif_scenes['pressure_contour'].append(get_frame(fig))
    plt.close(fig)

    fig, ax = plt.subplots()
    ax.semilogy(range(1, step+1), res_u, '-r', label='u')
    ax.semilogy(range(1, step+1), res_v, '-g', label='v')
    ax.semilogy(range(1, step+1), res_p, '-b', label='p')
    ax.legend()
    ax.grid(True)
    ax.set_title(f"Convergence History\nStep {step}, Iter {iter}, Time={time_s:.3f}s")
    gif_scenes['residuals'].append(get_frame(fig))
    plt.close(fig)


def get_frame(fig):
    fig.canvas.draw()
    # Use tobytes() which is standard in newer Matplotlib versions
    # tostring_rgb() is deprecated/removed in recent versions
    try:
        data = fig.canvas.tostring_rgb()
    except AttributeError:
        data = fig.canvas.tobytes()
        
    image = np.frombuffer(data, dtype='uint8')
    image = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
    return image


def create_gifs(gif_scenes):
    print("Saving GIFs...")
    for scene_name, frames in gif_scenes.items():
        if frames:
            filename = f"{scene_name}.gif"
            try:
                imageio.mimsave(filename, frames, fps=10, loop=0)
                print(f"Saved: {filename}")
            except Exception as e:
                print(f"Could not save {filename}: {e}")


def curl(X, Y, u, v):
    # Calculate vorticity
    dudy, dudx = np.gradient(u, Y[1,0]-Y[0,0], X[0,1]-X[0,0], axis=(0,1))
    dvdy, dvdx = np.gradient(v, Y[1,0]-Y[0,0], X[0,1]-X[0,0], axis=(0,1))
    return dvdx - dudy


def plot_final_results(X, Y, u, v, p, res_u, res_v, res_p):
    n = u.shape[0]
    fig, axs = plt.subplots(2, 3, figsize=(16, 10))

    # 1. Vectors
    axs[0,0].quiver(X[::3,::3], Y[::3,::3], u[::3,::3], v[::3,::3])
    axs[0,0].set_title("Velocity Vectors")
    axs[0,0].set_aspect('equal')

    # 2. Velocity Magnitude
    velMag = np.sqrt(u**2 + v**2)
    cs = axs[0,1].contourf(X, Y, velMag, 20, cmap='viridis')
    fig.colorbar(cs, ax=axs[0,1])
    axs[0,1].set_title("Velocity Magnitude")
    axs[0,1].set_aspect('equal')

    # 3. Pressure
    cs = axs[0,2].contourf(X, Y, p, 20, cmap='coolwarm')
    fig.colorbar(cs, ax=axs[0,2])
    axs[0,2].set_title("Pressure Field")
    axs[0,2].set_aspect('equal')

    # 4. Vorticity
    vort = curl(X, Y, u, v)
    cs = axs[1,0].contourf(X, Y, vort, 20, cmap='RdBu')
    fig.colorbar(cs, ax=axs[1,0])
    axs[1,0].set_title("Vorticity")
    axs[1,0].set_aspect('equal')

    # 5. Centerline Velocity Profiles
    # Vertical centerline (u velocity)
    axs[1,1].plot(u[:, n//2], Y[:, n//2], 'b-o', markersize=3, label='u (vertical cut)')
    # Horizontal centerline (v velocity)
    axs[1,1].plot(X[n//2, :], v[n//2, :], 'r-o', markersize=3, label='v (horizontal cut)')
    axs[1,1].legend()
    axs[1,1].grid(True)
    axs[1,1].set_title("Centerline Velocity Profiles")

    # 6. Residuals
    axs[1,2].semilogy(range(1, len(res_u)+1), res_u, 'r-', label='u')
    axs[1,2].semilogy(range(1, len(res_v)+1), res_v, 'g-', label='v')
    axs[1,2].semilogy(range(1, len(res_p)+1), res_p, 'b-', label='p')
    axs[1,2].legend()
    axs[1,2].grid(True)
    axs[1,2].set_title("Residuals")

    plt.tight_layout()
    plt.savefig('final_results.png')
    print("Final results saved to 'final_results.png'")
    # plt.show() # Commented out to prevent blocking in non-interactive environments


if __name__ == "__main__":
    VectorizedSolver()
