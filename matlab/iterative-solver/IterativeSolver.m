% IterativeSolver.m
%
% Corrected 2D Lid-Driven Cavity Pressure-Correction Solver
%
% Author: Ahmed Kandil (kandil.ahmed.amr@gmail.com)
% License: MIT (see LICENSE file in repository)
%
% Usage:
%   - Open in MATLAB
%   - Run IterativeSolver()
%   - Adjust parameters at the top of the file as needed
%   - GIFs and summary plots are generated in the working directory
%
% Features:
%   - Modular predictor / pressure-correction / corrector structure
%   - Separate GIF for each scene, saved with your original filenames
%   - Residual plots and final results summary
%   - Continuity residual added
%
% -------------------------------------------------------------------------

function IterativeSolver()

%% USER-ADJUSTABLE PARAMETERS
Re = 100;                % Reynolds number
L = 1.0;                 % Cavity length
n = 51;                  % Grid size (n x n)
dt = 0.001;              % Time step
total_time = 1.0;        % Total simulated time
alpha_u = 0.7;           % Under-relaxation for velocity
alpha_p = 0.3;           % Under-relaxation for pressure
tol = 1e-6;              % SIMPLE inner iteration tolerance
max_iter = 300;          % Max SIMPLE inner iterations per step
record_gif = true;       % Record GIFs? true/false
gif_stride = 1;          % Save every gif_stride time step
poisson_tol = 1e-6;      % Pressure Poisson tolerance
poisson_max = 500;       % Pressure Poisson max iterations

%% INITIALIZATION
clearvars -except Re L n dt total_time alpha_u alpha_p tol max_iter record_gif gif_stride poisson_tol poisson_max;
clc; close all;

nu = 1/Re;
dx = L/(n-1);
dy = dx;

x = linspace(0, L, n);
y = linspace(0, L, n);
[X, Y] = meshgrid(x, y);

max_steps = ceil(total_time/dt);

u = zeros(n, n);
v = zeros(n, n);
p = zeros(n, n);

% Lid-driven cavity top boundary
u(end, :) = 1.0;

res_u_arr    = zeros(max_steps, 1);
res_v_arr    = zeros(max_steps, 1);
res_p_arr    = zeros(max_steps, 1);
res_cont_arr = zeros(max_steps, 1);
simple_iter_arr = zeros(max_steps, 1);

% GIF names kept exactly as in your original code
gif_files = struct();
if record_gif
    gif_files.velocity_vectors = 'iterative_velocity_vectors.gif';
    gif_files.velocity_contour = 'iterative_velocity_contour.gif';
    gif_files.pressure_contour = 'iterative_pressure_contour.gif';
    gif_files.residuals        = 'iterative_residuals.gif';
    gif_files.streamlines      = 'iterative_streamlines.gif';

    delete_if_exists(gif_files.velocity_vectors);
    delete_if_exists(gif_files.velocity_contour);
    delete_if_exists(gif_files.pressure_contour);
    delete_if_exists(gif_files.residuals);
    delete_if_exists(gif_files.streamlines);
end

fprintf('Starting SIMPLE Lid Driven Cavity Simulation...\n');
fprintf('Grid size: %dx%d, Re: %d\n', n, n, Re);
fprintf('dt: %.4g, tol: %g, max_iter: %d\n', dt, tol, max_iter);
tic;

time = 0;
step = 0;

%% MAIN TIME STEPPING LOOP
while time < total_time + 1e-14
    step = step + 1;

    % SIMPLE Inner Iterations
    for iter = 1:max_iter
        u_prev = u;
        v_prev = v;
        p_prev = p;

        % Predictor
        [u_star, v_star] = predictor_step(u, v, p, dx, dy, dt, nu);

        % Pressure correction
        p_prime = solve_pressure_poisson(u_star, v_star, dx, dy, dt, poisson_tol, poisson_max);

        % Corrector
        [u_corr, v_corr, p_corr] = corrector_step(u_star, v_star, p, p_prime, dx, dy, dt, alpha_p);

        % Velocity under-relaxation
        u = alpha_u * u_corr + (1 - alpha_u) * u_prev;
        v = alpha_u * v_corr + (1 - alpha_u) * v_prev;
        p = p_corr;

        % Boundary conditions
        [u, v, p] = apply_boundary_conditions(u, v, p);

        % Residuals
        res_u = max(abs(u(:) - u_prev(:)));
        res_v = max(abs(v(:) - v_prev(:)));
        res_p = max(abs(p(:) - p_prev(:)));

        div = continuity_field(u, v, dx, dy);
        res_cont = max(abs(div(2:end-1, 2:end-1)), [], 'all');

        if max([res_u, res_v, res_p, res_cont]) < tol
            break;
        end
    end

    res_u_arr(step) = res_u;
    res_v_arr(step) = res_v;
    res_p_arr(step) = res_p;
    res_cont_arr(step) = res_cont;
    simple_iter_arr(step) = iter;

    % GIF capture and write directly to disk
    if record_gif && mod(step-1, gif_stride) == 0
        write_all_gifs(X, Y, u, v, p, ...
            res_u_arr(1:step), res_v_arr(1:step), res_p_arr(1:step), res_cont_arr(1:step), ...
            step, iter, time, L, gif_files);
    end

    % Advance time
    time = time + dt;

    if mod(step, 20) == 0 || step == 1
        fprintf('Step %4d / %4d | t = %.4f | SIMPLE iters = %3d | res = [%.3e %.3e %.3e %.3e]\n', ...
            step, max_steps, min(time, total_time), iter, res_u, res_v, res_p, res_cont);
    end
end

%% FINAL SUMMARY AND PLOTS
elapsedTime = toc;

res_u_arr    = res_u_arr(1:step);
res_v_arr    = res_v_arr(1:step);
res_p_arr    = res_p_arr(1:step);
res_cont_arr = res_cont_arr(1:step);
simple_iter_arr = simple_iter_arr(1:step);

fprintf('\nSimulation complete.\n');
fprintf('Elapsed real time: %.2f seconds (%.2f minutes).\n', elapsedTime, elapsedTime/60);
fprintf('Total time steps: %d, Final time: %.4f s\n', step, min(time, total_time));
fprintf('Average time per step: %.4f seconds\n', elapsedTime/step);

plot_final_results(X, Y, u, v, p, res_u_arr, res_v_arr, res_p_arr, res_cont_arr, Re);

end

% =========================================================================
function [u_star, v_star] = predictor_step(u, v, p, dx, dy, dt, nu)
n = size(u,1);

u_star = u;
v_star = v;

for j = 2:n-1
    for i = 2:n-1
        % u-momentum
        du2dx = ((u(j,i) + u(j,i+1))^2 - (u(j,i-1) + u(j,i))^2) / (4*dx);

        duvdy = ((v(j,i) + v(j,i+1))*(u(j,i) + u(j+1,i)) - ...
                 (v(j-1,i) + v(j-1,i+1))*(u(j-1,i) + u(j,i))) / (4*dy);

        d2udx2 = (u(j,i+1) - 2*u(j,i) + u(j,i-1)) / dx^2;
        d2udy2 = (u(j+1,i) - 2*u(j,i) + u(j-1,i)) / dy^2;

        dpdx = (p(j,i+1) - p(j,i-1)) / (2*dx);

        u_star(j,i) = u(j,i) + dt * (-du2dx - duvdy - dpdx + nu*(d2udx2 + d2udy2));

        % v-momentum
        dv2dy = ((v(j,i) + v(j+1,i))^2 - (v(j-1,i) + v(j,i))^2) / (4*dy);

        duvdx = ((u(j+1,i) + u(j,i))*(v(j,i+1) + v(j,i)) - ...
                 (u(j+1,i-1) + u(j,i-1))*(v(j,i) + v(j,i-1))) / (4*dx);

        d2vdx2 = (v(j,i+1) - 2*v(j,i) + v(j,i-1)) / dx^2;
        d2vdy2 = (v(j+1,i) - 2*v(j,i) + v(j-1,i)) / dy^2;

        dpdy = (p(j+1,i) - p(j-1,i)) / (2*dy);

        v_star(j,i) = v(j,i) + dt * (-duvdx - dv2dy - dpdy + nu*(d2vdx2 + d2vdy2));
    end
end
end

% =========================================================================
function p_prime = solve_pressure_poisson(u_star, v_star, dx, dy, dt, tol, max_iter)
n = size(u_star,1);
p_prime = zeros(n, n);

beta = 1.7; % SOR relaxation

for iter = 1:max_iter
    p_old = p_prime;

    for j = 2:n-1
        for i = 2:n-1
            rhs = ((u_star(j,i+1) - u_star(j,i-1)) / (2*dx) + ...
                   (v_star(j+1,i) - v_star(j-1,i)) / (2*dy)) / dt;

            p_new = ( ...
                (p_prime(j,i+1) + p_prime(j,i-1))*dy^2 + ...
                (p_prime(j+1,i) + p_prime(j-1,i))*dx^2 - ...
                rhs*dx^2*dy^2 ) / (2*(dx^2 + dy^2));

            p_prime(j,i) = (1 - beta)*p_prime(j,i) + beta*p_new;
        end
    end

    % Neumann BCs
    p_prime(1,:)   = p_prime(2,:);
    p_prime(end,:) = p_prime(end-1,:);
    p_prime(:,1)   = p_prime(:,2);
    p_prime(:,end) = p_prime(:,end-1);

    % Reference point
    p_prime(2,2) = 0;

    if max(abs(p_prime(:) - p_old(:))) < tol
        break;
    end
end
end

% =========================================================================
function [u, v, p] = corrector_step(u_star, v_star, p, p_prime, dx, dy, dt, alpha_p)
n = size(p,1);

u = u_star;
v = v_star;

for j = 2:n-1
    for i = 2:n-1
        u(j,i) = u_star(j,i) - dt * (p_prime(j,i+1) - p_prime(j,i-1)) / (2*dx);
        v(j,i) = v_star(j,i) - dt * (p_prime(j+1,i) - p_prime(j-1,i)) / (2*dy);
    end
end

p = p + alpha_p * p_prime;
end

% =========================================================================
function [u, v, p] = apply_boundary_conditions(u, v, p)
% Velocity BCs
u(:,1)   = 0;
u(:,end) = 0;
u(1,:)   = 0;
u(end,:) = 1;

v(:,1)   = 0;
v(:,end) = 0;
v(1,:)   = 0;
v(end,:) = 0;

% Pressure Neumann BCs
p(1,:)   = p(2,:);
p(end,:) = p(end-1,:);
p(:,1)   = p(:,2);
p(:,end) = p(:,end-1);

% Pressure reference
p(2,2) = 0;
end

% =========================================================================
function div = continuity_field(u, v, dx, dy)
n = size(u,1);
div = zeros(n,n);

for j = 2:n-1
    for i = 2:n-1
        div(j,i) = (u(j,i+1) - u(j,i-1)) / (2*dx) + ...
                   (v(j+1,i) - v(j-1,i)) / (2*dy);
    end
end
end

% =========================================================================
function write_all_gifs(X, Y, u, v, p, res_u, res_v, res_p, res_cont, step, iter, time, L, gif_files)

% Velocity vectors
h1 = figure('Visible','off', 'Color', 'w');
quiver(X(1:4:end,1:4:end), Y(1:4:end,1:4:end), ...
       u(1:4:end,1:4:end), v(1:4:end,1:4:end), 2, 'k');
title(sprintf('Velocity Vectors\nStep %d, SIMPLE iter %d, Time = %.3f s', step, iter, time));
xlabel('X'); ylabel('Y'); axis equal tight; grid on;
write_gif_frame(h1, gif_files.velocity_vectors, step == 1);
close(h1);

% Velocity magnitude contour
h2 = figure('Visible','off', 'Color', 'w');
velMag = sqrt(u.^2 + v.^2);
contourf(X, Y, velMag, 20, 'LineColor', 'none');
colorbar; axis equal tight;
title(sprintf('Velocity Magnitude\nStep %d, SIMPLE iter %d, Time = %.3f s', step, iter, time));
xlabel('X'); ylabel('Y');
write_gif_frame(h2, gif_files.velocity_contour, step == 1);
close(h2);

% Pressure contour
h3 = figure('Visible','off', 'Color', 'w');
contourf(X, Y, p, 20, 'LineColor', 'none');
colorbar; axis equal tight;
title(sprintf('Pressure Field\nStep %d, SIMPLE iter %d, Time = %.3f s', step, iter, time));
xlabel('X'); ylabel('Y');
write_gif_frame(h3, gif_files.pressure_contour, step == 1);
close(h3);

% Streamlines
h4 = figure('Visible','off', 'Color', 'w');
startx = linspace(0, L, 15);
starty = linspace(0, L, 15);
[sx, sy] = meshgrid(startx, starty);
streamline(X, Y, u, v, sx, sy);
axis equal tight; grid on;
title(sprintf('Streamlines\nStep %d, SIMPLE iter %d, Time = %.3f s', step, iter, time));
xlabel('X'); ylabel('Y');
write_gif_frame(h4, gif_files.streamlines, step == 1);
close(h4);

% Residual history
h5 = figure('Visible','off', 'Color', 'w');
semilogy(1:step, res_u, '-r', ...
         1:step, res_v, '-g', ...
         1:step, res_p, '-b', ...
         1:step, res_cont, '-k', 'LineWidth', 1.2);
xlabel('Time Step'); ylabel('Residual (log scale)'); grid on;
legend('u-res', 'v-res', 'p-res', 'cont-res', 'Location', 'northeast');
title(sprintf('Convergence History\nStep %d, SIMPLE iter %d, Time = %.3f s', step, iter, time));
write_gif_frame(h5, gif_files.residuals, step == 1);
close(h5);

end

% =========================================================================
function write_gif_frame(fig_handle, filename, is_first_frame)
drawnow;
frame = getframe(fig_handle);
img = frame2im(frame);
[A, map] = rgb2ind(img, 256);

if is_first_frame
    imwrite(A, map, filename, 'gif', 'LoopCount', Inf, 'DelayTime', 0.1);
else
    imwrite(A, map, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
end
end

% =========================================================================
function delete_if_exists(filename)
if exist(filename, 'file')
    delete(filename);
end
end

% =========================================================================
function plot_final_results(X, Y, u, v, p, res_u, res_v, res_p, res_cont, Re)
figure('Name','Final Results - Lid Driven Cavity', ...
       'Units','normalized','Position',[0.05 0.05 0.9 0.85], 'Color', 'w');

subplot(2,3,1);
quiver(X(1:3:end,1:3:end), Y(1:3:end,1:3:end), u(1:3:end,1:3:end), v(1:3:end,1:3:end), 1.5, 'k');
hold on;
startx = linspace(0, 1, 20);
starty = linspace(0, 1, 20);
[sx, sy] = meshgrid(startx, starty);
streamline(X, Y, u, v, sx, sy);
title('Velocity Vectors and Streamlines'); xlabel('X'); ylabel('Y'); axis equal tight; grid on;

subplot(2,3,2);
velMag = sqrt(u.^2 + v.^2);
contourf(X, Y, velMag, 20, 'LineColor', 'none'); colorbar;
title('Velocity Magnitude'); xlabel('X'); ylabel('Y'); axis equal tight;

subplot(2,3,3);
contourf(X, Y, p, 20, 'LineColor', 'none'); colorbar;
title('Pressure Field'); xlabel('X'); ylabel('Y'); axis equal tight;

subplot(2,3,4);
vorticity = compute_vorticity(u, v, X, Y);
contourf(X, Y, vorticity, 20, 'LineColor', 'none'); colorbar;
title('Vorticity Field'); xlabel('X'); ylabel('Y'); axis equal tight;

subplot(2,3,5);
mid = ceil(size(u,1)/2);
plot(Y(:,mid), u(:,mid), 'b-', 'LineWidth', 2); hold on;
plot(X(mid,:)', v(mid,:)', 'r-', 'LineWidth', 2);
title('Centerline Velocity Profiles');
xlabel('Position'); ylabel('Velocity');
legend('u on vertical centerline', 'v on horizontal centerline', 'Location', 'best');
grid on;

subplot(2,3,6);
semilogy(1:length(res_u), res_u, 'r-', ...
         1:length(res_v), res_v, 'g-', ...
         1:length(res_p), res_p, 'b-', ...
         1:length(res_cont), res_cont, 'k-', 'LineWidth', 1.2);
title('Convergence History'); xlabel('Time Step'); ylabel('Residual (log scale)');
legend('u-residual','v-residual','p-residual','cont-residual','Location','northeast'); grid on;

annotation('textbox', [0.02, 0.02, 0.3, 0.05], 'String', ...
    sprintf('Re = %d, Grid = %dx%d', Re, size(u,1), size(u,1)), ...
    'FitBoxToText', 'on', 'BackgroundColor', 'white');

saveas(gcf, 'final_results.png');
end

% =========================================================================
function vort = compute_vorticity(u, v, X, Y)
dx = X(1,2) - X(1,1);
dy = Y(2,1) - Y(1,1);

[dudy, ~] = gradient(u, dy, dx);
[~, dvdx] = gradient(v, dy, dx);

vort = dvdx - dudy;
end
