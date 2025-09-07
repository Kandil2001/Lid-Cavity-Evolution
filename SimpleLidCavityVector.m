function SIMPLE2D_LidDrivenCavity_Vectorized()
%% Simulation parameters
clear; clc; close all;
Re = 100;               % Reynolds number
L = 1.0;                % Cavity length
n = 101;                % Grid size (n x n)
dx = L/(n-1); dy = dx;  % Square cells
dt = 0.0005;            % Time step
nu = 1/Re;              % Kinematic viscosity
alpha_u = 0.5;          % Under-relaxation for velocity
alpha_p = 0.2;          % Under-relaxation for pressure
tol = 1e-5;             % SIMPLE convergence tolerance
max_iter = 500;         % Max SIMPLE iterations per time step
total_time = 2;        % Total simulation time

[X, Y] = meshgrid(0:dx:L, 0:dy:L);
u = zeros(n); v = zeros(n); p = zeros(n);
u(end,:) = 1;           % Lid moves to the right

time = 0; step = 0;
max_steps = ceil(total_time/dt);
res_u_arr = nan(max_steps,1); res_v_arr = nan(max_steps,1); res_p_arr = nan(max_steps,1);

figure('Name','SIMPLE Lid Driven Cavity (Vectorized)','Units','normalized','Position',[0.05 0.1 0.9 0.8]);
subplot(2,2,1); title('Velocity Vectors');
subplot(2,2,2); title('Velocity Magnitude Contour');
subplot(2,2,3); title('Pressure Contour');
subplot(2,2,4); title('Residuals');

fprintf('Starting SIMPLE Lid Driven Cavity Simulation (vectorized)...\n');
fprintf('SIMPLE tolerance: %g\n', tol);

tic;
while time < total_time
    step = step + 1;
    u_old = u; v_old = v; p_old = p;

    % --- SIMPLE Inner Iteration (with tolerance) ---
    for iter = 1:max_iter
        u_prev = u; v_prev = v; p_prev = p;
        [u_star, v_star] = predictor_step_vec(u, v, p, dx, dy, dt, nu, alpha_u);
        p_prime = solve_pressure_poisson_vec(u_star, v_star, dx, dy, dt, tol, max_iter);
        [u, v, p] = corrector_step_vec(u_star, v_star, p, p_prime, dx, dy, dt, alpha_p);

        res_u = max(max(abs(u - u_prev)));
        res_v = max(max(abs(v - v_prev)));
        res_p = max(max(abs(p - p_prev)));

        if max([res_u, res_v, res_p]) < tol
            break;
        end
    end

    res_u_arr(step) = res_u; res_v_arr(step) = res_v; res_p_arr(step) = res_p;

    % --- Real-Time Animation: Update every time step ---
    subplot(2,2,1); cla;
    Xs = X(1:4:end,1:4:end); Ys = Y(1:4:end,1:4:end);
    Us = u(1:4:end,1:4:end); Vs = v(1:4:end,1:4:end);
    quiver(Xs, Ys, Us, Vs, 2, 'k');
    title('Velocity Vectors'); xlabel('X'); ylabel('Y'); axis equal tight;

    subplot(2,2,2); cla;
    velMag = sqrt(u.^2 + v.^2);
    contourf(X, Y, velMag, 20, 'LineColor','none'); colorbar;
    title('Velocity Magnitude Contour'); xlabel('X'); ylabel('Y'); axis equal tight;

    subplot(2,2,3); cla;
    contourf(X, Y, p, 20, 'LineColor','none'); colorbar;
    title('Pressure Contour'); xlabel('X'); ylabel('Y'); axis equal tight;

    subplot(2,2,4); cla;
    semilogy(1:step, res_u_arr(1:step), '-r', 'DisplayName', 'u-res'); hold on;
    semilogy(1:step, res_v_arr(1:step), '-g', 'DisplayName', 'v-res');
    semilogy(1:step, res_p_arr(1:step), '-b', 'DisplayName', 'p-res');
    hold off;
    xlabel('Step'); ylabel('Residual');
    title('Residuals');
    legend('Location','northeast');
    xlim([1 max_steps]); ylim([1e-8 1]); grid on;

    drawnow;
    fprintf('Simulation time: %.4f / %.2f\n', time, total_time);

    time = time + dt;
end

elapsedTime = toc;
fprintf('Simulation complete.\n');
fprintf('Elapsed real time: %.2f seconds.\n', elapsedTime);
end

% Vectorized Predictor Step
function [u_star, v_star] = predictor_step_vec(u, v, p, dx, dy, dt, nu, alpha)
n = size(u,1);
u_star = u; v_star = v;
% Interior indices
idx = 2:n-1; idy = 2:n-1;

% u-momentum equation (excluding boundaries)
du2dx  = ((u(idy,idx) + u(idy,idx+1)).^2 - (u(idy,idx-1) + u(idy,idx)).^2) / (4*dx);
duvdy  = ((v(idy,idx) + v(idy,idx+1)) .* (u(idy,idx) + u(idy+1,idx)) ...
        - (v(idy-1,idx) + v(idy-1,idx+1)) .* (u(idy-1,idx) + u(idy,idx))) / (4*dy);
d2udx2 = (u(idy,idx+1) - 2*u(idy,idx) + u(idy,idx-1)) / dx^2;
d2udy2 = (u(idy+1,idx) - 2*u(idy,idx) + u(idy-1,idx)) / dy^2;
dpdx   = (p(idy,idx+1) - p(idy,idx)) / dx;
u_star(idy,idx) = u(idy,idx) + alpha * dt * (-du2dx - duvdy - dpdx + nu*(d2udx2 + d2udy2));

% v-momentum equation (excluding boundaries)
dv2dy  = ((v(idy,idx) + v(idy+1,idx)).^2 - (v(idy-1,idx) + v(idy,idx)).^2) / (4*dy);
duvdx  = ((u(idy+1,idx) + u(idy,idx)) .* (v(idy,idx+1) + v(idy,idx)) ...
        - (u(idy+1,idx-1) + u(idy,idx-1)) .* (v(idy,idx) + v(idy,idx-1))) / (4*dx);
d2vdx2 = (v(idy,idx+1) - 2*v(idy,idx) + v(idy,idx-1)) / dx^2;
d2vdy2 = (v(idy+1,idx) - 2*v(idy,idx) + v(idy-1,idx)) / dy^2;
dpdy   = (p(idy+1,idx) - p(idy,idx)) / dy;
v_star(idy,idx) = v(idy,idx) + alpha * dt * (-duvdx - dv2dy - dpdy + nu*(d2vdx2 + d2vdy2));
end

% Vectorized Pressure Poisson
function p_prime = solve_pressure_poisson_vec(u_star, v_star, dx, dy, dt, tol, max_iter)
n = size(u_star,1);
p_prime = zeros(n);
idx = 2:n-1; idy = 2:n-1;
for iter = 1:max_iter
    p_old = p_prime;
    rhs = ((u_star(idy,idx) - u_star(idy,idx-1))/dx + (v_star(idy,idx) - v_star(idy-1,idx))/dy)/dt;
    p_prime(idy,idx) = 0.25 * (p_prime(idy,idx+1) + p_prime(idy,idx-1) ...
                             + p_prime(idy+1,idx) + p_prime(idy-1,idx) - dx^2 * rhs);
    if max(max(abs(p_prime - p_old))) < tol
        break;
    end
end
end

% Vectorized Corrector Step
function [u, v, p] = corrector_step_vec(u_star, v_star, p, p_prime, dx, dy, dt, alpha)
n = size(p,1);
u = u_star; v = v_star;
idx = 2:n-1; idy = 2:n-1;
u(idy,idx) = u_star(idy,idx) - alpha * dt * (p_prime(idy,idx+1) - p_prime(idy,idx)) / dx;
v(idy,idx) = v_star(idy,idx) - alpha * dt * (p_prime(idy+1,idx) - p_prime(idy,idx)) / dy;
p = p + alpha * p_prime;

% Boundary conditions
u(:,1) = 0; u(:,end) = 0; u(1,:) = 0; u(end,:) = 1;
v(:,1) = 0; v(:,end) = 0; v(1,:) = 0; v(end,:) = 0;
end