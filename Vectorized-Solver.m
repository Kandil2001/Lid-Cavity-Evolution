function SIMPLE2D_LidDrivenCavity_Vectorized()
%% SIMPLE2D_LidDrivenCavity_Vectorized
% 2D Lid-Driven Cavity using SIMPLE algorithm (Vectorized Version)
% Structure and outputs are matched to the iterative implementation for fair comparison.
% - Modular function breakdown
% - Real-time visualization and GIF creation
% - Identical boundary conditions, plotting, and output metrics

% --- Simulation parameters ---
clear; clc; close all;
Re = 100;                  % Reynolds number
L = 1.0;                   % Cavity size
n = 151;                   % Grid size (n x n)
dx = L/(n-1); dy = dx;
dt = 0.0005;               % Time step
nu = 1/Re;                 % Kinematic viscosity
alpha_u = 0.5;             % Under-relaxation for velocity
alpha_p = 0.2;             % Under-relaxation for pressure
tol = 1e-5;                % SIMPLE convergence tolerance
max_iter = 500;            % Max SIMPLE iterations per time step
total_time = 2;            % Total simulation time

record_gif = true; gif_frame_interval = 5;

% --- Grid generation and variable initialization ---
[X, Y] = meshgrid(0:dx:L, 0:dy:L);
u = zeros(n); v = zeros(n); p = zeros(n);
u(end,:) = 1;              % Moving lid BC

time = 0; step = 0;
max_steps = ceil(total_time/dt);
res_u_arr = zeros(1000,1); res_v_arr = zeros(1000,1); res_p_arr = zeros(1000,1);

% --- GIF structures ---
if record_gif
    gif_scenes = struct();
    gif_scenes.velocity_vectors = struct('frames', [], 'filename', 'velocity_vectors_vec.gif');
    gif_scenes.velocity_contour = struct('frames', [], 'filename', 'velocity_contour_vec.gif');
    gif_scenes.pressure_contour = struct('frames', [], 'filename', 'pressure_contour_vec.gif');
    gif_scenes.residuals = struct('frames', [], 'filename', 'residuals_convergence_vec.gif');
    gif_scenes.streamlines = struct('frames', [], 'filename', 'streamlines_vec.gif');
end

% --- Figure setup for real-time visualization (matched layout) ---
hFig = figure('Name','SIMPLE Lid Driven Cavity (Vectorized)','Units','normalized',...
    'Position',[0.05 0.1 0.9 0.8], 'Color', 'w');
subplot(2,3,1); h_quiver = quiver(X(1:4:end,1:4:end), Y(1:4:end,1:4:end), u(1:4:end,1:4:end), v(1:4:end,1:4:end), 2, 'k');
title('Velocity Vectors'); xlabel('X'); ylabel('Y'); axis equal tight; grid on;
subplot(2,3,2); velMag = sqrt(u.^2 + v.^2); [~, h_contour_vel] = contourf(X, Y, velMag, 20, 'LineColor','none'); colorbar;
title('Velocity Magnitude Contour'); xlabel('X'); ylabel('Y'); axis equal tight;
subplot(2,3,3); [~, h_contour_p] = contourf(X, Y, p, 20, 'LineColor','none'); colorbar;
title('Pressure Contour'); xlabel('X'); ylabel('Y'); axis equal tight;
subplot(2,3,4); startx = linspace(0, L, 15); starty = linspace(0, L, 15); [sx, sy] = meshgrid(startx, starty);
h_stream = streamline(X, Y, u, v, sx, sy); set(h_stream, 'Color', 'b', 'LineWidth', 1);
title('Streamlines'); xlabel('X'); ylabel('Y'); axis equal tight; grid on;
subplot(2,3,5);
h_res_u = semilogy(1, 1e-1, '-r', 'LineWidth', 1.5); hold on;
h_res_v = semilogy(1, 1e-1, '-g', 'LineWidth', 1.5);
h_res_p = semilogy(1, 1e-1, '-b', 'LineWidth', 1.5); hold off;
xlabel('Time Step'); ylabel('Residual (log scale)');
title('Convergence History'); legend('u','v','p','Location','northeast');
xlim([1 max_steps]); ylim([1e-8 1]); grid on;
subplot(2,3,6); h_info = text(0.1, 0.9, sprintf('Time: %.3f s\nStep: %d\nRe: %d', time, step, Re), 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold');
axis off; title('Simulation Info');

% --- Simulation execution ---
fprintf('Starting SIMPLE Lid Driven Cavity Simulation (Vectorized)...\n');
fprintf('Grid size: %dx%d, Reynolds number: %d\n', n, n, Re);
fprintf('SIMPLE tolerance: %g, Max iterations: %d\n', tol, max_iter);

tic;
while time < total_time
    step = step + 1;
    if step > length(res_u_arr)
        new_size = 2 * length(res_u_arr);
        res_u_arr(new_size) = 0; res_v_arr(new_size) = 0; res_p_arr(new_size) = 0;
    end

    u_old = u; v_old = v; p_old = p;

    % --- SIMPLE Inner Iteration ---
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

    % --- Real-Time Animation: Identical subplot updates ---
    subplot(2,3,1); set(h_quiver, 'UData', u(1:4:end,1:4:end), 'VData', v(1:4:end,1:4:end));
    title(sprintf('Velocity Vectors (Time = %.3f s)', time));
    subplot(2,3,2); cla;
    velMag = sqrt(u.^2 + v.^2);
    contourf(X, Y, velMag, 20, 'LineColor','none'); colorbar;
    title(sprintf('Velocity Magnitude (Max = %.2f)', max(velMag(:))));
    subplot(2,3,3); cla;
    contourf(X, Y, p, 20, 'LineColor','none'); colorbar;
    title(sprintf('Pressure Field (Max = %.2f)', max(p(:))));
    subplot(2,3,4); delete(findobj(gca, 'Type', 'line'));
    h_stream = streamline(X, Y, u, v, sx, sy); set(h_stream, 'Color', 'b', 'LineWidth', 1);
    title('Streamlines');
    subplot(2,3,5);
    set(h_res_u, 'XData', 1:step, 'YData', res_u_arr(1:step));
    set(h_res_v, 'XData', 1:step, 'YData', res_v_arr(1:step));
    set(h_res_p, 'XData', 1:step, 'YData', res_p_arr(1:step));
    title(sprintf('Residuals (Iteration %d)', step)); xlim([1 max(step, 10)]);
    subplot(2,3,6); cla;
    text(0.1, 0.9, sprintf('Time: %.3f s\nStep: %d\nRe: %d\nResiduals:\n  u: %.2e\n  v: %.2e\n  p: %.2e', ...
        time, step, Re, res_u, res_v, res_p), 'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold');
    axis off; title('Simulation Info');

    % --- GIF capture ---
    if record_gif && mod(step, gif_frame_interval) == 0
        subplot(2,3,1); gif_scenes.velocity_vectors.frames = [gif_scenes.velocity_vectors.frames, getframe(gcf)];
        subplot(2,3,2); gif_scenes.velocity_contour.frames = [gif_scenes.velocity_contour.frames, getframe(gcf)];
        subplot(2,3,3); gif_scenes.pressure_contour.frames = [gif_scenes.pressure_contour.frames, getframe(gcf)];
        subplot(2,3,5); gif_scenes.residuals.frames = [gif_scenes.residuals.frames, getframe(gcf)];
        subplot(2,3,4); gif_scenes.streamlines.frames = [gif_scenes.streamlines.frames, getframe(gcf)];
    end

    drawnow;
    if mod(step, 10) == 0
        fprintf('Step: %d, Simulation time: %.4f/%.2f, Residuals: u=%.2e, v=%.2e, p=%.2e\n', ...
            step, time, total_time, res_u, res_v, res_p);
    end
    time = time + dt;
end

elapsedTime = toc;
fprintf('\nSimulation complete.\n');
fprintf('Elapsed real time: %.2f seconds (%.2f minutes).\n', elapsedTime, elapsedTime/60);
fprintf('Total time steps: %d, Final time: %.4f s\n', step, time);
fprintf('Average time per step: %.4f seconds\n', elapsedTime/step);

% --- GIF creation ---
if record_gif
    fprintf('Creating GIF files (Vectorized)...\n');
    create_gifs(gif_scenes);
end

% --- Final results and post-processing ---
plot_final_results(X, Y, u, v, p, res_u_arr(1:step), res_v_arr(1:step), res_p_arr(1:step));
end

% ---------------------- VECTORIZED FUNCTIONS ----------------------
function [u_star, v_star] = predictor_step_vec(u, v, p, dx, dy, dt, nu, alpha)
n = size(u,1);
u_star = u; v_star = v;
idx = 2:n-1; idy = 2:n-1;

du2dx  = ((u(idy,idx) + u(idy,idx+1)).^2 - (u(idy,idx-1) + u(idy,idx)).^2) / (4*dx);
duvdy  = ((v(idy,idx) + v(idy,idx+1)) .* (u(idy,idx) + u(idy+1,idx)) ...
        - (v(idy-1,idx) + v(idy-1,idx+1)) .* (u(idy-1,idx) + u(idy,idx))) / (4*dy);
d2udx2 = (u(idy,idx+1) - 2*u(idy,idx) + u(idy,idx-1)) / dx^2;
d2udy2 = (u(idy+1,idx) - 2*u(idy,idx) + u(idy-1,idx)) / dy^2;
dpdx   = (p(idy,idx+1) - p(idy,idx)) / dx;
u_star(idy,idx) = u(idy,idx) + alpha * dt * (-du2dx - duvdy - dpdx + nu*(d2udx2 + d2udy2));

dv2dy  = ((v(idy,idx) + v(idy+1,idx)).^2 - (v(idy-1,idx) + v(idy,idx)).^2) / (4*dy);
duvdx  = ((u(idy+1,idx) + u(idy,idx)) .* (v(idy,idx+1) + v(idy,idx)) ...
        - (u(idy+1,idx-1) + u(idy,idx-1)) .* (v(idy,idx) + v(idy,idx-1))) / (4*dx);
d2vdx2 = (v(idy,idx+1) - 2*v(idy,idx) + v(idy,idx-1)) / dx^2;
d2vdy2 = (v(idy+1,idx) - 2*v(idy,idx) + v(idy-1,idx)) / dy^2;
dpdy   = (p(idy+1,idx) - p(idy,idx)) / dy;
v_star(idy,idx) = v(idy,idx) + alpha * dt * (-duvdx - dv2dy - dpdy + nu*(d2vdx2 + d2vdy2));
end

function p_prime = solve_pressure_poisson_vec(u_star, v_star, dx, dy, dt, tol, max_iter)
n = size(u_star,1);
p_prime = zeros(n);
idx = 2:n-1; idy = 2:n-1;
for iter = 1:max_iter
    p_old = p_prime;
    rhs = ((u_star(idy,idx) - u_star(idy,idx-1))/dx + (v_star(idy,idx) - v_star(idy-1,idx))/dy)/dt;
    p_prime(idy,idx) = 0.25 * (p_prime(idy,idx+1) + p_prime(idy,idx-1) ...
                             + p_prime(idy+1,idx) + p_prime(idy-1,idx) - dx^2 * rhs);
    % Neumann BCs
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
u = u_star; v = v_star;
idx = 2:n-1; idy = 2:n-1;
u(idy,idx) = u_star(idy,idx) - alpha * dt * (p_prime(idy,idx+1) - p_prime(idy,idx)) / dx;
v(idy,idx) = v_star(idy,idx) - alpha * dt * (p_prime(idy+1,idx) - p_prime(idy,idx)) / dy;
p = p + alpha * p_prime;
% Velocity BCs
u(:,1) = 0; u(:,end) = 0; u(1,:) = 0; u(end,:) = 1;
v(:,1) = 0; v(:,end) = 0; v(1,:) = 0; v(end,:) = 0;
% Pressure BCs
p(1,:) = p(2,:); p(end,:) = p(end-1,:); p(:,1) = p(:,2); p(:,end) = p(:,end-1);
end

function create_gifs(gif_scenes)
    scenes = fieldnames(gif_scenes);
    for i = 1:length(scenes)
        scene_name = scenes{i};
        scene_data = gif_scenes.(scene_name);
        if ~isempty(scene_data.frames)
            for j = 1:length(scene_data.frames)
                [A, map] = rgb2ind(scene_data.frames(j).cdata, 256);
                if j == 1
                    imwrite(A, map, scene_data.filename, 'gif', 'LoopCount', Inf, 'DelayTime', 0.1);
                else
                    imwrite(A, map, scene_data.filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
                end
            end
        end
    end
end

function plot_final_results(X, Y, u, v, p, res_u, res_v, res_p)
figure('Name','Final Results - Lid Driven Cavity (Vectorized)',...
       'Units','normalized','Position',[0.05 0.05 0.9 0.85], 'Color', 'w');
subplot(2,3,1); quiver(X(1:3:end,1:3:end), Y(1:3:end,1:3:end), u(1:3:end,1:3:end), v(1:3:end,1:3:end), 1.5, 'k'); hold on;
startx = linspace(0, 1, 20); starty = linspace(0, 1, 20); [sx, sy] = meshgrid(startx, starty);
streamline(X, Y, u, v, sx, sy); title('Velocity Vectors and Streamlines'); xlabel('X'); ylabel('Y'); axis equal tight; grid on;
subplot(2,3,2); velMag = sqrt(u.^2 + v.^2); contourf(X, Y, velMag, 20, 'LineColor','none'); colorbar; title('Velocity Magnitude'); xlabel('X'); ylabel('Y'); axis equal tight;
subplot(2,3,3); contourf(X, Y, p, 20, 'LineColor','none'); colorbar; title('Pressure Field'); xlabel('X'); ylabel('Y'); axis equal tight;
subplot(2,3,4); vorticity = vorticity_field(X, Y, u, v); contourf(X, Y, vorticity, 20, 'LineColor','none'); colorbar; title('Vorticity Field'); xlabel('X'); ylabel('Y'); axis equal tight;
subplot(2,3,5); plot(u(ceil(end/2),:), Y(ceil(end/2),:), 'b-', 'LineWidth', 2); hold on; plot(u(:,ceil(end/2)), X(:,ceil(end/2)), 'r-', 'LineWidth', 2); title('Centerline Velocity Profiles'); xlabel('Velocity'); ylabel('Position'); legend('Vertical Centerline', 'Horizontal Centerline'); grid on;
subplot(2,3,6); semilogy(1:length(res_u), res_u, 'r-', 1:length(res_v), res_v, 'g-', 1:length(res_p), res_p, 'b-', 'LineWidth', 1.5); title('Convergence History'); xlabel('Time Step'); ylabel('Residual (log scale)'); legend('u','v','p','Location','northeast'); grid on;
annotation('textbox', [0.02, 0.02, 0.3, 0.05], 'String', sprintf('Re = %d, Grid = %dx%d', 100, size(u,1), size(u,2)), 'FitBoxToText', 'on', 'BackgroundColor', 'white');
saveas(gcf, 'final_results_vec.png');
end

function vort = vorticity_field(X, Y, u, v)
[dudy, dudx] = gradient(u, Y(1,2)-Y(1,1), X(2,1)-X(1,1));
[dvdy, dvdx] = gradient(v, Y(1,2)-Y(1,1), X(2,1)-X(1,1));
vort = dvdx - dudy;
end
