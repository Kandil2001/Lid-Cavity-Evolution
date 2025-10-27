import numpy as np
import matplotlib.pyplot as plt
import imageio
import time

# ---------------------------------------
# SIMPLE 2D Lid Driven Cavity Solver (Finite Volume, Staggered Grid)
# ---------------------------------------

def IterativeSolver():
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
    X, Y = np.meshgrid(np.linspace(0, L, n), np.linspace(0, L, n))
    u = np.zeros((n, n))
    v = np.zeros((n, n))
    p = np.zeros((n, n))

    # Lid boundary
    u[-1, :] = 1

    max_steps = int(np.ceil(total_time / dt))
    res_u_arr = np.zeros(max_steps)
    res_v_arr = np.zeros(max_steps)
    res_p_arr = np.zeros(max_steps)

    gif_scenes = {}
    if record_gif:
        gif_scenes = {
            'velocity_vectors': [],
            'velocity_contour': [],
            'pressure_contour': [],
            'streamlines': [],
            'residuals': []
        }

    print(f"Starting SIMPLE Lid Driven Cavity Simulation...\nGrid size: {n}x{n}, Re: {Re}")
    start_time = time.time()

    step = 0
    t = 0.0

    # MAIN LOOP
    while t < total_time:
        step += 1
        u_old, v_old, p_old = u.copy(), v.copy(), p.copy()

        # SIMPLE inner iterations
        for iter in range(1, max_iter + 1):
            u_prev, v_prev, p_prev = u.copy(), v.copy(), p.copy()

            u_star, v_star = predictor_step(u, v, p, dx, dy, dt, nu, alpha_u)
            p_prime = solve_pressure_poisson(u_star, v_star, dx, dy, dt, tol, max_iter)
            u, v, p = corrector_step(u_star, v_star, p, p_prime, dx, dy, dt, alpha_p)

            # Boundary conditions
            u[:, 0] = u[:, -1] = 0
            u[0, :] = 0
            u[-1, :] = 1
            v[:, 0] = v[:, -1] = 0
            v[0, :] = v[-1, :] = 0
            p[0, :] = p[1, :]
            p[-1, :] = p[-2, :]
            p[:, 0] = p[:, 1]
            p[:, -1] = p[:, -2]

            res_u = np.max(np.abs(u - u_prev))
            res_v = np.max(np.abs(v - v_prev))
            res_p = np.max(np.abs(p - p_prev))

            if max(res_u, res_v, res_p) < tol:
                break

        res_u_arr[step - 1] = res_u
        res_v_arr[step - 1] = res_v
        res_p_arr[step - 1] = res_p

        # GIF capture
        if record_gif:
            capture_frames(gif_scenes, X, Y, u, v, p, res_u_arr[:step], res_v_arr[:step], res_p_arr[:step], step, iter, t)

        t += dt

    # Save GIFs
    if record_gif:
        print("Creating GIF files...")
        create_gifs(gif_scenes)

    elapsed = time.time() - start_time
    print(f"\nSimulation complete.\nElapsed: {elapsed:.2f} s ({elapsed/60:.2f} min)")
    print(f"Total steps: {step}, Final time: {t:.4f} s, Avg. step time: {elapsed/step:.4f} s")

    plot_final_results(X, Y, u, v, p, res_u_arr[:step], res_v_arr[:step], res_p_arr[:step])


def predictor_step(u, v, p, dx, dy, dt, nu, alpha):
    n = u.shape[0]
    u_star = u.copy()
    v_star = v.copy()
    inv_4dx = 1.0 / (4 * dx)
    inv_4dy = 1.0 / (4 * dy)
    inv_dx_sq = 1.0 / dx**2
    inv_dy_sq = 1.0 / dy**2
    alpha_dt = alpha * dt

    # u-momentum
    for j in range(1, n - 1):
        for i in range(1, n - 1):
            du2dx = ((u[j,i] + u[j,i+1])**2 - (u[j,i-1] + u[j,i])**2) * inv_4dx
            duvdy = ((v[j,i] + v[j,i+1]) * (u[j,i] + u[j+1,i]) -
                     (v[j-1,i] + v[j-1,i+1]) * (u[j-1,i] + u[j,i])) * inv_4dy
            d2udx2 = (u[j,i+1] - 2*u[j,i] + u[j,i-1]) * inv_dx_sq
            d2udy2 = (u[j+1,i] - 2*u[j,i] + u[j-1,i]) * inv_dy_sq
            dpdx = (p[j,i+1] - p[j,i]) / dx
            u_star[j,i] = u[j,i] + alpha_dt * (-du2dx - duvdy - dpdx + nu * (d2udx2 + d2udy2))

    # v-momentum
    for j in range(1, n - 1):
        for i in range(1, n - 1):
            dv2dy = ((v[j,i] + v[j+1,i])**2 - (v[j-1,i] + v[j,i])**2) * inv_4dy
            duvdx = ((u[j+1,i] + u[j,i]) * (v[j,i+1] + v[j,i]) -
                     (u[j+1,i-1] + u[j,i-1]) * (v[j,i] + v[j,i-1])) * inv_4dx
            d2vdx2 = (v[j,i+1] - 2*v[j,i] + v[j,i-1]) * inv_dx_sq
            d2vdy2 = (v[j+1,i] - 2*v[j,i] + v[j-1,i]) * inv_dy_sq
            dpdy = (p[j+1,i] - p[j,i]) / dy
            v_star[j,i] = v[j,i] + alpha_dt * (-duvdx - dv2dy - dpdy + nu * (d2vdx2 + d2vdy2))

    return u_star, v_star


def solve_pressure_poisson(u_star, v_star, dx, dy, dt, tol, max_iter):
    n = u_star.shape[0]
    p_prime = np.zeros((n, n))
    inv_dx = 1.0 / dx
    inv_dy = 1.0 / dy
    dt_rhs_factor = 1.0 / dt
    laplacian_factor = 0.25

    for _ in range(max_iter):
        p_old = p_prime.copy()
        for j in range(1, n - 1):
            for i in range(1, n - 1):
                rhs = ((u_star[j,i] - u_star[j,i-1]) * inv_dx +
                       (v_star[j,i] - v_star[j-1,i]) * inv_dy) * dt_rhs_factor
                p_prime[j,i] = laplacian_factor * (p_prime[j,i+1] + p_prime[j,i-1] +
                                                   p_prime[j+1,i] + p_prime[j-1,i] - dx**2 * rhs)

        # Neumann BC
        p_prime[0,:] = p_prime[1,:]
        p_prime[-1,:] = p_prime[-2,:]
        p_prime[:,0] = p_prime[:,1]
        p_prime[:,-1] = p_prime[:, -2]

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
            u[j,i] = u_star[j,i] - alpha_dt_dx * (p_prime[j,i+1] - p_prime[j,i])
            v[j,i] = v_star[j,i] - alpha_dt_dy * (p_prime[j+1,i] - p_prime[j,i])

    p = p + alpha * p_prime
    return u, v, p


def capture_frames(gif_scenes, X, Y, u, v, p, res_u, res_v, res_p, step, iter, time_s):
    # Velocity vectors
    fig, ax = plt.subplots()
    ax.quiver(X[::4,::4], Y[::4,::4], u[::4,::4], v[::4,::4])
    ax.set_title(f"Velocity Vectors\nStep {step}, Iter {iter}, Time={time_s:.3f}s")
    gif_scenes['velocity_vectors'].append(get_frame(fig))
    plt.close(fig)

    # Velocity magnitude contour
    fig, ax = plt.subplots()
    velMag = np.sqrt(u**2 + v**2)
    cs = ax.contourf(X, Y, velMag, 20)
    fig.colorbar(cs)
    ax.set_title(f"Velocity Magnitude\nStep {step}, Iter {iter}, Time={time_s:.3f}s")
    gif_scenes['velocity_contour'].append(get_frame(fig))
    plt.close(fig)

    # Pressure contour
    fig, ax = plt.subplots()
    cs = ax.contourf(X, Y, p, 20)
    fig.colorbar(cs)
    ax.set_title(f"Pressure Field\nStep {step}, Iter {iter}, Time={time_s:.3f}s")
    gif_scenes['pressure_contour'].append(get_frame(fig))
    plt.close(fig)

    # Streamlines
    fig, ax = plt.subplots()
    startx = np.linspace(0, 1, 15)
    starty = np.linspace(0, 1, 15)
    sx, sy = np.meshgrid(startx, starty)
    ax.streamplot(X, Y, u, v, start_points=np.array([sx.flatten(), sy.flatten()]).T)
    ax.set_title(f"Streamlines\nStep {step}, Iter {iter}, Time={time_s:.3f}s")
    gif_scenes['streamlines'].append(get_frame(fig))
    plt.close(fig)

    # Residuals
    fig, ax = plt.subplots()
    ax.semilogy(range(1, step+1), res_u, '-r', label='u')
    ax.semilogy(range(1, step+1), res_v, '-g', label='v')
    ax.semilogy(range(1, step+1), res_p, '-b', label='p')
    ax.legend()
    ax.set_title(f"Convergence History\nStep {step}, Iter {iter}, Time={time_s:.3f}s")
    gif_scenes['residuals'].append(get_frame(fig))
    plt.close(fig)


def get_frame(fig):
    fig.canvas.draw()
    image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
    image = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
    return image


def create_gifs(gif_scenes):
    for scene_name, frames in gif_scenes.items():
        if frames:
            filename = f"{scene_name}.gif"
            imageio.mimsave(filename, frames, duration=0.1)
            print(f"Saved: {filename}")


def plot_final_results(X, Y, u, v, p, res_u, res_v, res_p):
    n = u.shape[0]
    fig, axs = plt.subplots(2, 3, figsize=(15, 10))

    axs[0,0].quiver(X[::3,::3], Y[::3,::3], u[::3,::3], v[::3,::3])
    axs[0,0].set_title("Velocity Vectors")

    velMag = np.sqrt(u**2 + v**2)
    cs = axs[0,1].contourf(X, Y, velMag, 20)
    fig.colorbar(cs, ax=axs[0,1])
    axs[0,1].set_title("Velocity Magnitude")

    cs = axs[0,2].contourf(X, Y, p, 20)
    fig.colorbar(cs, ax=axs[0,2])
    axs[0,2].set_title("Pressure Field")

    vort = curl(X, Y, u, v)
    cs = axs[1,0].contourf(X, Y, vort, 20)
    fig.colorbar(cs, ax=axs[1,0])
    axs[1,0].set_title("Vorticity")

    axs[1,1].plot(u[n//2,:], Y[n//2,:], 'b-', label='Vertical')
    axs[1,1].plot(u[:,n//2], X[:,n//2], 'r-', label='Horizontal')
    axs[1,1].legend()
    axs[1,1].set_title("Centerline Velocity")

    axs[1,2].semilogy(range(1, len(res_u)+1), res_u, 'r-', label='u')
    axs[1,2].semilogy(range(1, len(res_v)+1), res_v, 'g-', label='v')
    axs[1,2].semilogy(range(1, len(res_p)+1), res_p, 'b-', label='p')
    axs[1,2].legend()
    axs[1,2].set_title("Residuals")

    plt.tight_layout()
    plt.savefig('final_results.png')
    plt.show()


def curl(X, Y, u, v):
    dudy, dudx = np.gradient(u, Y[1,0]-Y[0,0], X[0,1]-X[0,0])
    dvdy, dvdx = np.gradient(v, Y[1,0]-Y[0,0], X[0,1]-X[0,0])
    return dvdx - dudy


if __name__ == "__main__":
    IterativeSolver()
