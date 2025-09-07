function SIMPLE2D_LidDrivenCavity()
%% Simulation parameters
clear; clc; close all;
Re = 100;               % Reynolds number
L = 1.0;                % Cavity length
n = 151;                % Grid size (n x n)
dx = L/(n-1); dy = dx;  % Square cells
dt = 0.0005;            % Time step
nu = 1/Re;              % Kinematic viscosity
alpha_u = 0.5;          % Under-relaxation for velocity
alpha_p = 0.2;          % Under-relaxation for pressure
tol = 1e-5;             % SIMPLE convergence tolerance
max_iter = 500;         % Max SIMPLE iterations per time step
total_time = 2;        % Total simulation time

%% Grid and variables
[X, Y] = meshgrid(0:dx:L, 0:dy:L);
u = zeros(n); v = zeros(n); p = zeros(n);
u(end,:) = 1;           % Lid moves to the right

time = 0; step = 0;
max_steps = ceil(total_time/dt);
res_u_arr = nan(max_steps,1); res_v_arr = nan(max_steps,1); res_p_arr = nan(max_steps,1);

%% Figure setup
hFig = figure('Name','SIMPLE Lid Driven Cavity','Units','normalized','Position',[0.05 0.1 0.9 0.8]);
subplot(2,2,1); title('Velocity Vectors');
subplot(2,2,2); title('Velocity Magnitude Contour');
subplot(2,2,3); title('Pressure Contour');
subplot(2,2,4); title('Residuals');

fprintf('Starting SIMPLE Lid Driven Cavity Simulation...\n');
fprintf('SIMPLE tolerance: %g\n', tol);

tic;  % Start measuring real time

%% Time stepping loop
while time < total_time
    step = step + 1;
    u_old = u; v_old = v; p_old = p;

    % --- SIMPLE Inner Iteration (with tolerance) ---
    for iter = 1:max_iter
        u_prev = u; v_prev = v; p_prev = p;
        [u_star, v_star] = predictor_step(u, v, p, dx, dy, dt, nu, alpha_u);
        p_prime = solve_pressure_poisson(u_star, v_star, dx, dy, dt, tol, max_iter);
        [u, v, p] = corrector_step(u_star, v_star, p, p_prime, dx, dy, dt, alpha_p);

        res_u = max(max(abs(u - u_prev)));
        res_v = max(max(abs(v - v_prev)));
        res_p = max(max(abs(p - p_prev)));

        if max([res_u, res_v, res_p]) < tol
            break;
        end
    end

    res_u_arr(step) = res_u; res_v_arr(step) = res_v; res_p_arr(step) = res_p;

    % --- Real-Time Animation: Update every time step ---
    % Velocity vector plot (Top-Left)
    subplot(2,2,1); cla;
    Xs = X(1:4:end,1:4:end); Ys = Y(1:4:end,1:4:end);
    Us = u(1:4:end,1:4:end); Vs = v(1:4:end,1:4:end);
    quiver(Xs, Ys, Us, Vs, 2, 'k');
    title('Velocity Vectors'); xlabel('X'); ylabel('Y'); axis equal tight;

    % Velocity magnitude contour (Top-Right)
    subplot(2,2,2); cla;
    velMag = sqrt(u.^2 + v.^2);
    contourf(X, Y, velMag, 20, 'LineColor','none'); colorbar;
    title('Velocity Magnitude Contour'); xlabel('X'); ylabel('Y'); axis equal tight;

    % Pressure contour (Bottom-Left)
    subplot(2,2,3); cla;
    contourf(X, Y, p, 20, 'LineColor','none'); colorbar;
    title('Pressure Contour'); xlabel('X'); ylabel('Y'); axis equal tight;

    % Residual subplot (Bottom-Right)
    subplot(2,2,4); cla;
    semilogy(1:step, res_u_arr(1:step), '-r', 'DisplayName', 'u-res'); hold on;
    semilogy(1:step, res_v_arr(1:step), '-g', 'DisplayName', 'v-res');
    semilogy(1:step, res_p_arr(1:step), '-b', 'DisplayName', 'p-res');
    hold off;
    xlabel('Step'); ylabel('Residual');
    title('Residuals');
    legend('Location','northeast');
    xlim([1 max_steps]); ylim([1e-8 1]); grid on;

    drawnow;   % Force figure window update

    fprintf('Simulation time: %.4f / %.2f\n', time, total_time);

    time = time + dt;
end

elapsedTime = toc;  % End measuring real time

fprintf('Simulation complete.\n');
fprintf('Elapsed real time: %.2f seconds.\n', elapsedTime);
end

%% Predictor step
function [u_star, v_star] = predictor_step(u, v, p, dx, dy, dt, nu, alpha)
n = size(u,1);
u_star = u; v_star = v;
for j = 2:n-1
    for i = 2:n-1
        du2dx = ((u(j,i)+u(j,i+1))^2 - (u(j,i-1)+u(j,i))^2)/(4*dx);
        duvdy = ((v(j,i)+v(j,i+1))*(u(j,i)+u(j+1,i)) - ...
                 (v(j-1,i)+v(j-1,i+1))*(u(j-1,i)+u(j,i)))/(4*dy);
        d2udx2 = (u(j,i+1)-2*u(j,i)+u(j,i-1))/dx^2;
        d2udy2 = (u(j+1,i)-2*u(j,i)+u(j-1,i))/dy^2;
        dpdx = (p(j,i+1) - p(j,i))/dx;
        u_star(j,i) = u(j,i) + alpha * dt * (-du2dx - duvdy - dpdx + nu*(d2udx2 + d2udy2));
    end
end
for j = 2:n-1
    for i = 2:n-1
        dv2dy = ((v(j,i)+v(j+1,i))^2 - (v(j-1,i)+v(j,i))^2)/(4*dy);
        duvdx = ((u(j+1,i)+u(j,i))*(v(j,i+1)+v(j,i)) - ...
                 (u(j+1,i-1)+u(j,i-1))*(v(j,i)+v(j,i-1)))/(4*dx);
        d2vdx2 = (v(j,i+1)-2*v(j,i)+v(j,i-1))/dx^2;
        d2vdy2 = (v(j+1,i)-2*v(j,i)+v(j-1,i))/dy^2;
        dpdy = (p(j+1,i) - p(j,i))/dy;
        v_star(j,i) = v(j,i) + alpha * dt * (-duvdx - dv2dy - dpdy + nu*(d2vdx2 + d2vdy2));
    end
end
end

%% Pressure Poisson solver
function p_prime = solve_pressure_poisson(u_star, v_star, dx, dy, dt, tol, max_iter)
n = size(u_star,1);
p_prime = zeros(n);
for iter = 1:max_iter
    p_old = p_prime;
    for j = 2:n-1
        for i = 2:n-1
            rhs = ((u_star(j,i) - u_star(j,i-1))/dx + (v_star(j,i) - v_star(j-1,i))/dy)/dt;
            p_prime(j,i) = 0.25 * (p_prime(j,i+1) + p_prime(j,i-1) + ...
                                   p_prime(j+1,i) + p_prime(j-1,i) - dx^2 * rhs);
        end
    end
    if max(max(abs(p_prime - p_old))) < tol
        break;
    end
end
end

%% Corrector step
function [u, v, p] = corrector_step(u_star, v_star, p, p_prime, dx, dy, dt, alpha)
n = size(p,1);
u = u_star; v = v_star;
for j = 2:n-1
    for i = 2:n-1
        u(j,i) = u_star(j,i) - alpha * dt * (p_prime(j,i+1) - p_prime(j,i)) / dx;
        v(j,i) = v_star(j,i) - alpha * dt * (p_prime(j+1,i) - p_prime(j,i)) / dy;
    end
end
p = p + alpha * p_prime;
% Boundary conditions
u(:,1) = 0; u(:,end) = 0; u(1,:) = 0; u(end,:) = 1;
v(:,1) = 0; v(:,end) = 0; v(1,:) = 0; v(end,:) = 0;
end