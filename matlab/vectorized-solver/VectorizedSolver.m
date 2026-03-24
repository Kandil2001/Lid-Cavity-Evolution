% VectorizedSolver.m
%
% True Staggered-Grid 2D Lid-Driven Cavity Solver
% SIMPLE-style pressure-correction solver (vectorized MATLAB version)
%
% Author: Ahmed Kandil (kandil.ahmed.amr@gmail.com)
% License: MIT (see LICENSE file in repository)
%
% Usage:
%   - Open in MATLAB
%   - Run VectorizedSolver()
%   - Adjust parameters at the top of the file as needed
%   - GIFs and summary plots are generated in the working directory
%
% Storage:
%   p : cell centers        -> (N x N)
%   u : vertical cell faces -> (N x (N+1))
%   v : horizontal faces    -> ((N+1) x N)
%
% Outputs:
%   vectorized_velocity_vectors.gif
%   vectorized_velocity_contour.gif
%   vectorized_pressure_contour.gif
%   vectorized_residuals.gif
%   vectorized_streamlines.gif
%   final_results_vectorized.png
%
% -------------------------------------------------------------------------

function VectorizedSolver()

%% USER-ADJUSTABLE PARAMETERS
Re = 100;                % Reynolds number
L = 1.0;                 % Cavity length
N = 51;                  % Number of pressure cells in each direction
dt = 5e-4;               % Time step
total_time = 1.0;        % Total simulated time

alpha_u = 0.7;           % Under-relaxation for velocity
alpha_p = 0.3;           % Under-relaxation for pressure

tol = 1e-6;              % Outer SIMPLE tolerance
max_iter = 200;          % Max SIMPLE iterations per time step

poisson_tol = 1e-6;      % Pressure Poisson tolerance
poisson_max = 800;       % Pressure Poisson max iterations

record_gif = true;       % Record GIFs? true/false
gif_stride = 1;          % Save every gif_stride step

%% INITIALIZATION
clearvars -except Re L N dt total_time alpha_u alpha_p tol max_iter poisson_tol poisson_max record_gif gif_stride;
clc; close all;

nu = 1 / Re;
dx = L / N;
dy = L / N;

max_steps = ceil(total_time / dt);

% Cell-center coordinates for pressure / plotting
xc = linspace(dx/2, L - dx/2, N);
yc = linspace(dy/2, L - dy/2, N);
[Xc, Yc] = meshgrid(xc, yc);

% True staggered storage
p = zeros(N, N);        % cell centers
u = zeros(N, N+1);      % vertical faces
v = zeros(N+1, N);      % horizontal faces

[u, v, p] = apply_boundary_conditions(u, v, p);

% Residual arrays
res_u_arr    = zeros(max_steps, 1);
res_v_arr    = zeros(max_steps, 1);
res_p_arr    = zeros(max_steps, 1);
res_cont_arr = zeros(max_steps, 1);
iter_arr     = zeros(max_steps, 1);

% GIF filenames
gif_files = struct();
if record_gif
    gif_files.velocity_vectors = 'vectorized_velocity_vectors.gif';
    gif_files.velocity_contour = 'vectorized_velocity_contour.gif';
    gif_files.pressure_contour = 'vectorized_pressure_contour.gif';
    gif_files.residuals        = 'vectorized_residuals.gif';
    gif_files.streamlines      = 'vectorized_streamlines.gif';

    delete_if_exists(gif_files.velocity_vectors);
    delete_if_exists(gif_files.velocity_contour);
    delete_if_exists(gif_files.pressure_contour);
    delete_if_exists(gif_files.residuals);
    delete_if_exists(gif_files.streamlines);
end

fprintf('Starting VECTORlZED true staggered-grid lid-driven cavity simulation...\n');
fprintf('Pressure cells: %dx%d\n', N, N);
fprintf('u size: %dx%d | v size: %dx%d | p size: %dx%d\n', size(u,1), size(u,2), size(v,1), size(v,2), size(p,1), size(p,2));
fprintf('Re = %g, dt = %g, total_time = %g\n', Re, dt, total_time);
fprintf('alpha_u = %g, alpha_p = %g\n', alpha_u, alpha_p);
fprintf('tol = %g, max_iter = %d\n\n', tol, max_iter);

tic;
time = 0.0;
step = 0;

%% MAIN TIME LOOP
while time < total_time - 1e-14
    step = step + 1;

    for iter = 1:max_iter
        u_prev = u;
        v_prev = v;
        p_prev = p;

        [u_star, v_star] = predictor_step_staggered_vectorized(u, v, p, dx, dy, dt, nu);
        p_prime = solve_pressure_poisson_staggered_vectorized(u_star, v_star, dx, dy, dt, poisson_tol, poisson_max);
        [u_corr, v_corr, p_corr] = corrector_step_staggered_vectorized(u_star, v_star, p, p_prime, dx, dy, dt, alpha_p);

        % Velocity under-relaxation
        u = alpha_u * u_corr + (1 - alpha_u) * u_prev;
        v = alpha_u * v_corr + (1 - alpha_u) * v_prev;
        p = p_corr;

        [u, v, p] = apply_boundary_conditions(u, v, p);

        % Residuals
        res_u = max(abs(u(:) - u_prev(:)));
        res_v = max(abs(v(:) - v_prev(:)));
        res_p = max(abs(p(:) - p_prev(:)));

        div = continuity_field_staggered(u, v, dx, dy);
        res_cont = max(abs(div(:)));

        if max([res_u, res_v, res_p, res_cont]) < tol
            break;
        end
    end

    res_u_arr(step)    = res_u;
    res_v_arr(step)    = res_v;
    res_p_arr(step)    = res_p;
    res_cont_arr(step) = res_cont;
    iter_arr(step)     = iter;

    [uc, vc] = staggered_to_cell_center(u, v);

    if record_gif && mod(step-1, gif_stride) == 0
        write_all_gifs(Xc, Yc, uc, vc, p, ...
            res_u_arr(1:step), res_v_arr(1:step), res_p_arr(1:step), res_cont_arr(1:step), ...
            step, iter, time, L, gif_files);
    end

    time = time + dt;

    if mod(step, 20) == 0 || step == 1
        fprintf('Step %4d / %4d | t = %.4f | iters = %3d | res = [%.3e %.3e %.3e %.3e]\n', ...
            step, max_steps, min(time, total_time), iter, res_u, res_v, res_p, res_cont);
    end
end

elapsedTime = toc;

% Trim histories
res_u_arr    = res_u_arr(1:step);
res_v_arr    = res_v_arr(1:step);
res_p_arr    = res_p_arr(1:step);
res_cont_arr = res_cont_arr(1:step);
iter_arr     = iter_arr(1:step);

[uc, vc] = staggered_to_cell_center(u, v);

fprintf('\nSimulation complete.\n');
fprintf('Elapsed real time: %.2f seconds (%.2f minutes)\n', elapsedTime, elapsedTime/60);
fprintf('Total time steps: %d\n', step);
fprintf('Final time: %.6f s\n', min(time, total_time));
fprintf('Average time per step: %.4f s\n', elapsedTime / step);

plot_final_results(Xc, Yc, uc, vc, p, res_u_arr, res_v_arr, res_p_arr, res_cont_arr, Re, N);

end

%% =========================================================================
function [u, v, p] = apply_boundary_conditions(u, v, p)
% u: N x (N+1)
% v: (N+1) x N
% p: N x N

% Normal velocity at walls
u(:,1)   = 0;   % left wall
u(:,end) = 0;   % right wall
v(1,:)   = 0;   % bottom wall
v(end,:) = 0;   % top wall

% Tangential velocity
u(1,:)   = 0;   % bottom wall
u(end,:) = 1;   % moving lid
u(end,1) = 0;
u(end,end) = 0;

v(:,1)   = 0;   % left wall tangential approx
v(:,end) = 0;   % right wall tangential approx

% Pressure Neumann BC
p(1,:)   = p(2,:);
p(end,:) = p(end-1,:);
p(:,1)   = p(:,2);
p(:,end) = p(:,end-1);

% Reference pressure
p(1,1) = 0;
end

%% =========================================================================
function [u_star, v_star] = predictor_step_staggered_vectorized(u, v, p, dx, dy, dt, nu)

N = size(p,1);

u_star = u;
v_star = v;

%% --- u-momentum on vertical faces: rows 2:N-1, cols 2:N
J = 2:N-1;
I = 2:N;

ue = 0.5 * (u(J,I)   + u(J,I+1));
uw = 0.5 * (u(J,I-1) + u(J,I));

un = 0.5 * (u(J,I)   + u(J+1,I));
us = 0.5 * (u(J-1,I) + u(J,I));

vn = 0.5 * (v(J+1,I-1) + v(J+1,I));
vs = 0.5 * (v(J,  I-1) + v(J,  I));

du2dx = (ue.^2 - uw.^2) / dx;
duvdy = (vn .* un - vs .* us) / dy;

d2udx2 = (u(J,I+1) - 2*u(J,I) + u(J,I-1)) / dx^2;
d2udy2 = (u(J+1,I) - 2*u(J,I) + u(J-1,I)) / dy^2;

dpdx = (p(J,I) - p(J,I-1)) / dx;

u_star(J,I) = u(J,I) + dt * (-du2dx - duvdy - dpdx + nu * (d2udx2 + d2udy2));

%% --- v-momentum on horizontal faces: rows 2:N, cols 2:N-1
J = 2:N;
I = 2:N-1;

vn = 0.5 * (v(J,I)   + v(J+1,I));
vs = 0.5 * (v(J-1,I) + v(J,I));

ve = 0.5 * (v(J,I)   + v(J,I+1));
vw = 0.5 * (v(J,I-1) + v(J,I));

ue = 0.5 * (u(J-1,I+1) + u(J,I+1));
uw = 0.5 * (u(J-1,I)   + u(J,I));

dv2dy = (vn.^2 - vs.^2) / dy;
duvdx = (ue .* ve - uw .* vw) / dx;

d2vdx2 = (v(J,I+1) - 2*v(J,I) + v(J,I-1)) / dx^2;
d2vdy2 = (v(J+1,I) - 2*v(J,I) + v(J-1,I)) / dy^2;

dpdy = (p(J,I) - p(J-1,I)) / dy;

v_star(J,I) = v(J,I) + dt * (-duvdx - dv2dy - dpdy + nu * (d2vdx2 + d2vdy2));

end

%% =========================================================================
function p_prime = solve_pressure_poisson_staggered_vectorized(u_star, v_star, dx, dy, dt, tol, max_iter)

N = size(u_star,1);
p_prime = zeros(N, N);

beta = 1.7;

rhs = continuity_field_staggered(u_star, v_star, dx, dy) / dt;

for iter = 1:max_iter
    p_old = p_prime;

    J = 2:N-1;
    I = 2:N-1;

    p_new = ( ...
        (p_prime(J,I+1) + p_prime(J,I-1)) * dy^2 + ...
        (p_prime(J+1,I) + p_prime(J-1,I)) * dx^2 - ...
        rhs(J,I) * dx^2 * dy^2 ) ./ (2 * (dx^2 + dy^2));

    p_prime(J,I) = (1 - beta) * p_prime(J,I) + beta * p_new;

    % Neumann BC
    p_prime(1,:)   = p_prime(2,:);
    p_prime(end,:) = p_prime(end-1,:);
    p_prime(:,1)   = p_prime(:,2);
    p_prime(:,end) = p_prime(:,end-1);

    % Reference point
    p_prime(1,1) = 0;

    if max(abs(p_prime(:) - p_old(:))) < tol
        break;
    end
end
end

%% =========================================================================
function [u, v, p] = corrector_step_staggered_vectorized(u_star, v_star, p, p_prime, dx, dy, dt, alpha_p)

N = size(p,1);

u = u_star;
v = v_star;

% Correct u on vertical faces: cols 2:N
u(:,2:N) = u_star(:,2:N) - dt * (p_prime(:,2:N) - p_prime(:,1:N-1)) / dx;

% Correct v on horizontal faces: rows 2:N
v(2:N,:) = v_star(2:N,:) - dt * (p_prime(2:N,:) - p_prime(1:N-1,:)) / dy;

p = p + alpha_p * p_prime;
end

%% =========================================================================
function div = continuity_field_staggered(u, v, dx, dy)
% Divergence at pressure cell centers
% u: N x (N+1), v: (N+1) x N

div = (u(:,2:end) - u(:,1:end-1)) / dx + ...
      (v(2:end,:) - v(1:end-1,:)) / dy;
end

%% =========================================================================
function [uc, vc] = staggered_to_cell_center(u, v)
uc = 0.5 * (u(:,1:end-1) + u(:,2:end));
vc = 0.5 * (v(1:end-1,:) + v(2:end,:));
end

%% =========================================================================
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

%% =========================================================================
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

%% =========================================================================
function delete_if_exists(filename)
if exist(filename, 'file')
    delete(filename);
end
end

%% =========================================================================
function plot_final_results(X, Y, u, v, p, res_u, res_v, res_p, res_cont, Re, N)

figure('Name','Final Results - Lid Driven Cavity (Vectorized)', ...
       'Units','normalized','Position',[0.05 0.05 0.9 0.85], 'Color', 'w');

subplot(2,3,1);
quiver(X(1:3:end,1:3:end), Y(1:3:end,1:3:end), ...
       u(1:3:end,1:3:end), v(1:3:end,1:3:end), 1.5, 'k');
hold on;
startx = linspace(min(X(:)), max(X(:)), 20);
starty = linspace(min(Y(:)), max(Y(:)), 20);
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

annotation('textbox', [0.02, 0.02, 0.35, 0.05], 'String', ...
    sprintf('Re = %d, Pressure grid = %dx%d', Re, N, N), ...
    'FitBoxToText', 'on', 'BackgroundColor', 'white');

saveas(gcf, 'final_results_vectorized.png');
end

%% =========================================================================
function vort = compute_vorticity(u, v, X, Y)
dx = X(1,2) - X(1,1);
dy = Y(2,1) - Y(1,1);

[dudy, ~] = gradient(u, dy, dx);
[~, dvdx] = gradient(v, dy, dx);

vort = dvdx - dudy;
end
