% VectorizedSolver.m
%
% SIMPLE 2D Lid-Driven Cavity Solver (Finite Volume, Staggered Grid)
% Vectorized version for enhanced MATLAB performance.
%
% Author: Ahmed Kandil (kandil.ahmed.amr@gmail.com)
% License: MIT
%
% Usage:
%   - Open in MATLAB
%   - Run VectorizedSolver()
%   - Adjust parameters at the top as needed
%   - GIFs and summary plots are generated in the working directory

function VectorizedSolver()
%% USER-ADJUSTABLE PARAMETERS
Re = 100;                % Reynolds number (try 100, 400, 1000)
L = 1.0;                 % Cavity length
n = 51;                  % Grid size (n x n, e.g., 31/51/101)
dt = 0.002;              % Time step
total_time = 1.0;        % Total simulated time
alpha_u = 0.7;           % Under-relaxation for velocity (0.5 - 0.8 typical)
alpha_p = 0.3;           % Under-relaxation for pressure (0.2 - 0.5 typical)
tol = 1e-6;              % SIMPLE inner iteration tolerance
max_iter = 300;          % Max SIMPLE inner iterations per step
record_gif = true;       % Record GIFs? true/false

%% INITIALIZATION
clearvars -except Re L n dt total_time alpha_u alpha_p tol max_iter record_gif;
clc; close all;

nu = 1/Re;
dx = L/(n-1); dy = dx;

max_steps = ceil(total_time/dt);
res_u_arr = zeros(1, max_steps);
res_v_arr = zeros(1, max_steps);
res_p_arr = zeros(1, max_steps);

[X, Y] = meshgrid(0:dx:L, 0:dy:L);
u = zeros(n); v = zeros(n); p = zeros(n);
u(end,:) = 1; % Lid moves right

if record_gif
    gif_scenes.velocity_vectors = struct('frames', [], 'filename', 'vectorized_velocity_vectors.gif');
    gif_scenes.velocity_contour = struct('frames', [], 'filename', 'vectorized_velocity_contour.gif');
    gif_scenes.pressure_contour = struct('frames', [], 'filename', 'vectorized_pressure_contour.gif');
    gif_scenes.residuals        = struct('frames', [], 'filename', 'vectorized_residuals.gif');
    gif_scenes.streamlines      = struct('frames', [], 'filename', 'vectorized_streamlines.gif');
end

fprintf('Starting Vectorized SIMPLE Lid Driven Cavity Simulation...\n');
fprintf('Grid size: %dx%d, Re: %d\n', n, n, Re);
fprintf('dt: %.4g, tol: %g, max_iter: %d\n', dt, tol, max_iter);
tic;
time = 0; step = 0;

%% MAIN TIME STEPPING LOOP
while time < total_time
    step = step + 1;
    u_old = u; v_old = v; p_old = p;

    % SIMPLE Inner Iterations (vectorized)
    for iter = 1:max_iter
        u_prev = u; v_prev = v; p_prev = p;
        [u_star, v_star] = predictor_step_vec(u, v, p, dx, dy, dt, nu, alpha_u);
        p_prime = solve_pressure_poisson_vec(u_star, v_star, dx, dy, dt, tol, max_iter);
        [u, v, p] = corrector_step_vec(u_star, v_star, p, p_prime, dx, dy, dt, alpha_p);

        % Boundary conditions
        u(:,1) = 0; u(:,end) = 0; u(1,:) = 0; u(end,:) = 1;
        v(:,1) = 0; v(:,end) = 0; v(1,:) = 0; v(end,:) = 0;
        p(1,:) = p(2,:);   p(end,:) = p(end-1,:);
        p(:,1) = p(:,2);   p(:,end) = p(:,end-1);

        % SIMPLE residuals
        res_u = max(max(abs(u - u_prev)));
        res_v = max(max(abs(v - v_prev)));
        res_p = max(max(abs(p - p_prev)));
        if max([res_u, res_v, res_p]) < tol
            break;
        end
    end

    res_u_arr(step) = res_u; res_v_arr(step) = res_v; res_p_arr(step) = res_p;

    % GIF CAPTURE (same as before, just change filenames)
    if record_gif
        h1 = figure('Visible','off');
        quiver(X(1:4:end,1:4:end), Y(1:4:end,1:4:end), u(1:4:end,1:4:end), v(1:4:end,1:4:end), 2, 'k');
        title(sprintf('Velocity Vectors\nStep %d, SIMPLE iter %d, Time = %.3f s', step, iter, time));
        xlabel('X'); ylabel('Y'); axis equal tight; grid on;
        gif_scenes.velocity_vectors.frames = [gif_scenes.velocity_vectors.frames, getframe(h1)];
        close(h1);

        h2 = figure('Visible','off');
        velMag = sqrt(u.^2 + v.^2);
        contourf(X, Y, velMag, 20, 'LineColor','none');
        colorbar; axis equal tight;
        title(sprintf('Velocity Magnitude\nStep %d, SIMPLE iter %d, Time = %.3f s', step, iter, time));
        xlabel('X'); ylabel('Y');
        gif_scenes.velocity_contour.frames = [gif_scenes.velocity_contour.frames, getframe(h2)];
        close(h2);

        h3 = figure('Visible','off');
        contourf(X, Y, p, 20, 'LineColor','none');
        colorbar; axis equal tight;
        title(sprintf('Pressure Field\nStep %d, SIMPLE iter %d, Time = %.3f s', step, iter, time));
        xlabel('X'); ylabel('Y');
        gif_scenes.pressure_contour.frames = [gif_scenes.pressure_contour.frames, getframe(h3)];
        close(h3);

        h4 = figure('Visible','off');
        startx = linspace(0, L, 15); starty = linspace(0, L, 15);
        [sx, sy] = meshgrid(startx, starty);
        streamline(X, Y, u, v, sx, sy);
        axis equal tight; grid on;
        title(sprintf('Streamlines\nStep %d, SIMPLE iter %d, Time = %.3f s', step, iter, time));
        xlabel('X'); ylabel('Y');
        gif_scenes.streamlines.frames = [gif_scenes.streamlines.frames, getframe(h4)];
        close(h4);

        h5 = figure('Visible','off');
        semilogy(1:step, res_u_arr(1:step), '-r', 1:step, res_v_arr(1:step), '-g', 1:step, res_p_arr(1:step), '-b');
        xlabel('Time Step'); ylabel('Residual (log scale)'); grid on;
        legend('u-res','v-res','p-res','Location','northeast');
        title(sprintf('Convergence History\nStep %d, SIMPLE iter %d, Time = %.3f s', step, iter, time));
        gif_scenes.residuals.frames = [gif_scenes.residuals.frames, getframe(h5)];
        close(h5);
    end

    time = time + dt;
end

%% CREATE GIFs FROM FRAMES
if record_gif
    fprintf('Creating GIF files...\n');
    create_gifs(gif_scenes);
end

%% FINAL SUMMARY AND PLOTS
elapsedTime = toc;
fprintf('\nSimulation complete.\n');
fprintf('Elapsed real time: %.2f seconds (%.2f minutes).\n', elapsedTime, elapsedTime/60);
fprintf('Total time steps: %d, Final time: %.4f s\n', step, time);
fprintf('Average time per step: %.4f seconds\n', elapsedTime/step);

plot_final_results(X, Y, u, v, p, res_u_arr(1:step), res_v_arr(1:step), res_p_arr(1:step));
end

function create_gifs(gif_scenes)
    scenes = fieldnames(gif_scenes);
    for i = 1:length(scenes)
        scene_name = scenes{i};
        scene_data = gif_scenes.(scene_name);
        if ~isempty(scene_data.frames)
            fprintf('Creating %s...\n', scene_data.filename);
            for j = 1:length(scene_data.frames)
                [A, map] = rgb2ind(scene_data.frames(j).cdata, 256);
                if j == 1
                    imwrite(A, map, scene_data.filename, 'gif', ...
                           'LoopCount', Inf, 'DelayTime', 0.1);
                else
                    imwrite(A, map, scene_data.filename, 'gif', ...
                           'WriteMode', 'append', 'DelayTime', 0.1);
                end
            end
            fprintf('Saved: %s\n', scene_data.filename);
        end
    end
end

function [u_star, v_star] = predictor_step_vec(u, v, p, dx, dy, dt, nu, alpha)
n = size(u,1);
inv_4dx = 1/(4*dx); inv_4dy = 1/(4*dy);
inv_dx_sq = 1/dx^2; inv_dy_sq = 1/dy^2;
alpha_dt = alpha * dt;
u_star = u; v_star = v;

% Interior indices
ii = 2:n-1; jj = 2:n-1;

% For u
du2dx = ((u(jj,ii) + u(jj,ii+1)).^2 - (u(jj,ii-1) + u(jj,ii)).^2) * inv_4dx;
duvdy = ((v(jj,ii) + v(jj,ii+1)) .* (u(jj,ii) + u(jj+1,ii)) ...
        - (v(jj-1,ii) + v(jj-1,ii+1)) .* (u(jj-1,ii) + u(jj,ii))) * inv_4dy;
d2udx2 = (u(jj,ii+1) - 2*u(jj,ii) + u(jj,ii-1)) * inv_dx_sq;
d2udy2 = (u(jj+1,ii) - 2*u(jj,ii) + u(jj-1,ii)) * inv_dy_sq;
dpdx = (p(jj,ii+1) - p(jj,ii)) / dx;
u_star(jj,ii) = u(jj,ii) + alpha_dt * (-du2dx - duvdy - dpdx + nu*(d2udx2 + d2udy2));

% For v
dv2dy = ((v(jj,ii) + v(jj+1,ii)).^2 - (v(jj-1,ii) + v(jj,ii)).^2) * inv_4dy;
duvdx = ((u(jj+1,ii) + u(jj,ii)) .* (v(jj,ii+1) + v(jj,ii)) ...
        - (u(jj+1,ii-1) + u(jj,ii-1)) .* (v(jj,ii) + v(jj,ii-1))) * inv_4dx;
d2vdx2 = (v(jj,ii+1) - 2*v(jj,ii) + v(jj,ii-1)) * inv_dx_sq;
d2vdy2 = (v(jj+1,ii) - 2*v(jj,ii) + v(jj-1,ii)) * inv_dy_sq;
dpdy = (p(jj+1,ii) - p(jj,ii)) / dy;
v_star(jj,ii) = v(jj,ii) + alpha_dt * (-duvdx - dv2dy - dpdy + nu*(d2vdx2 + d2vdy2));
end

function p_prime = solve_pressure_poisson_vec(u_star, v_star, dx, dy, dt, tol, max_iter)
n = size(u_star,1);
p_prime = zeros(n);
inv_dx = 1/dx; inv_dy = 1/dy;
dt_rhs_factor = 1/dt;
laplacian_factor = 0.25;
ii = 2:n-1; jj = 2:n-1;
for iter = 1:max_iter
    p_old = p_prime;
    rhs = ((u_star(jj,ii) - u_star(jj,ii-1))*inv_dx + (v_star(jj,ii) - v_star(jj-1,ii))*inv_dy) * dt_rhs_factor;
    % Jacobi update (vectorized)
    p_prime(jj,ii) = laplacian_factor * (p_prime(jj,ii+1) + p_prime(jj,ii-1) + ...
                                         p_prime(jj+1,ii) + p_prime(jj-1,ii) - dx^2 * rhs);
    % Homogeneous Neumann BCs
    p_prime(1,:) = p_prime(2,:);
    p_prime(end,:) = p_prime(end-1,:);
    p_prime(:,1) = p_prime(:,2);
    p_prime(:,end) = p_prime(:,end-1);

    if max(max(abs(p_prime - p_old))) < tol
        break;
    end
end
end

function [u, v, p] = corrector_step_vec(u_star, v_star, p, p_prime, dx, dy, dt, alpha)
n = size(p,1);
alpha_dt_dx = alpha * dt / dx; alpha_dt_dy = alpha * dt / dy;
ii = 2:n-1; jj = 2:n-1;
u = u_star; v = v_star;
u(jj,ii) = u_star(jj,ii) - alpha_dt_dx * (p_prime(jj,ii+1) - p_prime(jj,ii));
v(jj,ii) = v_star(jj,ii) - alpha_dt_dy * (p_prime(jj+1,ii) - p_prime(jj,ii));
p = p + alpha * p_prime;
end

function plot_final_results(X, Y, u, v, p, res_u, res_v, res_p)
figure('Name','Final Results - Lid Driven Cavity',...
       'Units','normalized','Position',[0.05 0.05 0.9 0.85], 'Color', 'w');
subplot(2,3,1);
quiver(X(1:3:end,1:3:end), Y(1:3:end,1:3:end), u(1:3:end,1:3:end), v(1:3:end,1:3:end), 1.5, 'k');
hold on;
startx = linspace(0, 1, 20); starty = linspace(0, 1, 20);
[sx, sy] = meshgrid(startx, starty);
streamline(X, Y, u, v, sx, sy);
title('Velocity Vectors and Streamlines'); xlabel('X'); ylabel('Y'); axis equal tight; grid on;
subplot(2,3,2);
velMag = sqrt(u.^2 + v.^2);
contourf(X, Y, velMag, 20, 'LineColor','none'); colorbar;
title('Velocity Magnitude'); xlabel('X'); ylabel('Y'); axis equal tight;
subplot(2,3,3);
contourf(X, Y, p, 20, 'LineColor','none'); colorbar;
title('Pressure Field'); xlabel('X'); ylabel('Y'); axis equal tight;
subplot(2,3,4);
vorticity = curl(X, Y, u, v);
contourf(X, Y, vorticity, 20, 'LineColor','none'); colorbar;
title('Vorticity Field'); xlabel('X'); ylabel('Y'); axis equal tight;
subplot(2,3,5);
plot(u(ceil(size(u,1)/2),:), Y(ceil(size(Y,1)/2),:), 'b-', 'LineWidth', 2);
hold on; plot(u(:,ceil(size(u,2)/2)), X(:,ceil(size(X,2)/2)), 'r-', 'LineWidth', 2);
title('Centerline Velocity Profiles'); xlabel('Velocity'); ylabel('Position');
legend('Vertical Centerline', 'Horizontal Centerline'); grid on;
subplot(2,3,6);
semilogy(1:length(res_u), res_u, 'r-', 1:length(res_v), res_v, 'g-', 1:length(res_p), res_p, 'b-');
title('Convergence History'); xlabel('Time Step'); ylabel('Residual (log scale)');
legend('u-residual','v-residual','p-residual','Location','northeast'); grid on;
annotation('textbox', [0.02, 0.02, 0.3, 0.05], 'String', ...
    sprintf('Re = %d, Grid = %dx%d', 100, size(u,1), size(u,1)), ...
    'FitBoxToText', 'on', 'BackgroundColor', 'white');
saveas(gcf, 'final_results.png');
end

function vort = curl(X, Y, u, v)
[dudy, dudx] = gradient(u, Y(1,2)-Y(1,1), X(2,1)-X(1,1));
[dvdy, dvdx] = gradient(v, Y(1,2)-Y(1,1), X(2,1)-X(1,1));
vort = dvdx - dudy;
end
