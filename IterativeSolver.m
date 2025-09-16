function IterativeSolver()
%% SIMPLE2D_LidDrivenCavity
% Solves the 2D Lid-Driven Cavity problem using the SIMPLE algorithm.
% This implementation uses triple-nested loops for educational clarity,
% establishing a performance baseline for future optimized versions.
%
% Key Features:
%   - Finite Volume discretization on a staggered grid
%   - Predictor-Corrector method for pressure-velocity coupling
%   - Iterative solution of the Pressure Poisson equation
%   - Silent visualization for GIF creation
%
% OUTPUTS:
%   Elapsed real time is printed to the console, providing the benchmark metric.
%   GIF files are saved for each visualization scene.
%
% BENCHMARK: On a 101x101 grid, Re=100, this implementation serves as a
%            performance baseline for comparison with optimized versions.
%
% -------------------------------------------------------------------------

%% Simulation parameters
clear; clc; close all;

% Physical parameters
Re = 100;               % Reynolds number
L = 1.0;                % Cavity length

% Numerical parameters
n = 101;                % Grid size (101x101 - good balance of accuracy/speed)
dx = L/(n-1); 
dy = dx;                % Square cells
dt = 0.002;             % Time step (more stable for Re=100)
nu = 1/Re;              % Kinematic viscosity
alpha_u = 0.7;          % Under-relaxation for velocity (increased for better convergence)
alpha_p = 0.3;          % Under-relaxation for pressure (increased for better convergence)
tol = 1e-6;             % Stricter tolerance for better convergence
max_iter = 1000;        % More iterations for pressure convergence
total_time = 1.0;       % Longer simulation time for better development

% GIF recording parameters
record_gif = true;      % Set to true to record GIFs

%% Optimized Calculations and Pre-allocation 🚀
% Pre-calculate constant terms for efficiency
inv_4dx = 1/(4*dx);
inv_4dy = 1/(4*dy);
inv_dx_sq = 1/dx^2;
inv_dy_sq = 1/dy^2;
inv_dt_dx_dy = 1/(dt * dx * dy);
inv_dt_dx = 1/(dt*dx);
inv_dt_dy = 1/(dt*dy);
alpha_dt = alpha_u * dt;
alpha_dt_dx = alpha_u * dt / dx;
alpha_dt_dy = alpha_u * dt / dy;

% Pre-allocate arrays to their final size to avoid re-allocation inside loops
max_steps = ceil(total_time/dt);
res_u_arr = zeros(1, max_steps); 
res_v_arr = zeros(1, max_steps); 
res_p_arr = zeros(1, max_steps);

%% Grid generation and variable initialization
[X, Y] = meshgrid(0:dx:L, 0:dy:L);
u = zeros(n);           % x-velocity component
v = zeros(n);           % y-velocity component
p = zeros(n);           % Pressure field

% Boundary conditions: No-slip on all walls except moving lid
u(end,:) = 1;           % Lid moves to the right with unit velocity

% Initialize time tracking and residuals storage
time = 0; 
step = 0;

% Initialize GIF structures if recording
if record_gif
    gif_scenes = struct();
    gif_scenes.velocity_vectors = struct('frames', [], 'filename', 'iterative_velocity_vectors.gif');
    gif_scenes.velocity_contour = struct('frames', [], 'filename', 'iterative_velocity_contour.gif');
    gif_scenes.pressure_contour = struct('frames', [], 'filename', 'iterative_pressure_contour.gif');
    gif_scenes.residuals = struct('frames', [], 'filename', 'iterative_residuals.gif');
    gif_scenes.streamlines = struct('frames', [], 'filename', 'iterative_streamlines.gif');
end

%% Simulation execution
fprintf('Starting SIMPLE Lid Driven Cavity Simulation...\n');
fprintf('Grid size: %dx%d, Reynolds number: %d\n', n, n, Re);
fprintf('SIMPLE tolerance: %g, Max iterations: %d\n', tol, max_iter);
fprintf('Starting time integration. This may take a while...\n');

tic;  % Start measuring real time

%% Time stepping loop
while time < total_time
    step = step + 1;
    
    % Store previous time step values
    u_old = u; 
    v_old = v; 
    p_old = p;

    % --- SIMPLE Inner Iteration (with tolerance) ---
    for iter = 1:max_iter
        u_prev = u; 
        v_prev = v; 
        p_prev = p;
        
        % Predictor step: Solve momentum equations for intermediate velocities
        [u_star, v_star] = predictor_step(u, v, p, dx, dy, dt, nu, alpha_u);
        
        % Pressure correction: Solve pressure Poisson equation
        p_prime = solve_pressure_poisson(u_star, v_star, dx, dy, dt, tol, max_iter);
        
        % Corrector step: Update velocities and pressure
        [u, v, p] = corrector_step(u_star, v_star, p, p_prime, dx, dy, dt, alpha_p);
        
        % Calculate residuals
        res_u = max(max(abs(u - u_prev)));
        res_v = max(max(abs(v - v_prev)));
        res_p = max(max(abs(p - p_prev)));

        % Check for convergence
        if max([res_u, res_v, res_p]) < tol
            break;
        end
    end

    % Store residuals for plotting
    res_u_arr(step) = res_u;
    res_v_arr(step) = res_v;
    res_p_arr(step) = res_p;

    % Capture frames for GIFs at each time step
    if record_gif
        % Create off-screen figure for capturing frames
        hFig = figure('Visible', 'off', 'Position', [100, 100, 1200, 800], 'Color', 'w');
        
        % Velocity vectors subplot
        subplot(2,3,1);
        quiver(X(1:4:end,1:4:end), Y(1:4:end,1:4:end), ...
                     u(1:4:end,1:4:end), v(1:4:end,1:4:end), 2, 'k');
        title(sprintf('Velocity Vectors (Time = %.3f s)', time)); 
        xlabel('X'); ylabel('Y'); 
        axis equal tight;
        grid on;
        
        % Velocity magnitude contour subplot
        subplot(2,3,2);
        velMag = sqrt(u.^2 + v.^2);
        if max(velMag(:)) - min(velMag(:)) > eps
            contourf(X, Y, velMag, 20, 'LineColor','none');
        else
            imagesc(velMag);
            set(gca, 'YDir', 'normal');
        end
        colorbar;
        title(sprintf('Velocity Magnitude (Time = %.3f s)', time)); 
        xlabel('X'); ylabel('Y'); 
        axis equal tight;
        
        % Pressure contour subplot
        subplot(2,3,3);
        if max(p(:)) - min(p(:)) > eps
            contourf(X, Y, p, 20, 'LineColor','none');
        else
            imagesc(p);
            set(gca, 'YDir', 'normal');
        end
        colorbar;
        title(sprintf('Pressure Field (Time = %.3f s)', time)); 
        xlabel('X'); ylabel('Y'); 
        axis equal tight;
        
        % Streamlines subplot
        subplot(2,3,4);
        startx = linspace(0, L, 15);
        starty = linspace(0, L, 15);
        [sx, sy] = meshgrid(startx, starty);
        h_stream = streamline(X, Y, u, v, sx, sy);
        set(h_stream, 'Color', 'b', 'LineWidth', 1);
        title(sprintf('Streamlines (Time = %.3f s)', time));
        xlabel('X'); ylabel('Y');
        axis equal tight;
        grid on;
        
        % Residuals subplot
        subplot(2,3,5);
        semilogy(1:step, res_u_arr(1:step), '-r', 'DisplayName', 'u-res', 'LineWidth', 1.5);
        hold on;
        semilogy(1:step, res_v_arr(1:step), '-g', 'DisplayName', 'v-res', 'LineWidth', 1.5);
        semilogy(1:step, res_p_arr(1:step), '-b', 'DisplayName', 'p-res', 'LineWidth', 1.5);
        hold off;
        xlabel('Time Step'); 
        ylabel('Residual (log scale)');
        title(sprintf('Convergence History (Time = %.3f s)', time));
        legend('Location','northeast');
        xlim([1 max(step, 10)]); 
        ylim([1e-8 1]);  % Show residuals down to 1e-8
        grid on;
        
        % Info subplot
        subplot(2,3,6);
        text(0.1, 0.9, sprintf('Time: %.3f s\nStep: %d\nRe: %d\nResiduals:\n  u: %.2e\n  v: %.2e\n  p: %.2e', ...
              time, step, Re, res_u, res_v, res_p), ...
             'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold');
        axis off;
        title('Simulation Info');
        
        % Capture frames for each scene
        % Velocity vectors scene
        subplot(2,3,1);
        frame = getframe(hFig);
        gif_scenes.velocity_vectors.frames = [gif_scenes.velocity_vectors.frames, frame];
        
        % Velocity contour scene
        subplot(2,3,2);
        frame = getframe(hFig);
        gif_scenes.velocity_contour.frames = [gif_scenes.velocity_contour.frames, frame];
        
        % Pressure contour scene
        subplot(2,3,3);
        frame = getframe(hFig);
        gif_scenes.pressure_contour.frames = [gif_scenes.pressure_contour.frames, frame];
        
        % Residuals scene
        subplot(2,3,5);
        frame = getframe(hFig);
        gif_scenes.residuals.frames = [gif_scenes.residuals.frames, frame];
        
        % Streamlines scene
        subplot(2,3,4);
        frame = getframe(hFig);
        gif_scenes.streamlines.frames = [gif_scenes.streamlines.frames, frame];
        
        % Close the figure
        close(hFig);
    end

    % Advance time
    time = time + dt;
end

%% Create GIFs from captured frames
if record_gif
    fprintf('Creating GIF files...\n');
    create_gifs(gif_scenes);
end

%% Post-processing and results
elapsedTime = toc;  % End measuring real time

fprintf('\nSimulation complete.\n');
fprintf('Elapsed real time: %.2f seconds (%.2f minutes).\n', elapsedTime, elapsedTime/60);
fprintf('Total time steps: %d, Final time: %.4f s\n', step, time);

% Calculate and display performance metrics
fprintf('Average time per step: %.4f seconds\n', elapsedTime/step);

% Plot final results in a new figure
plot_final_results(X, Y, u, v, p, res_u_arr(1:step), res_v_arr(1:step), res_p_arr(1:step));

end

%% Function to create GIF files
function create_gifs(gif_scenes)
    scenes = fieldnames(gif_scenes);
    
    for i = 1:length(scenes)
        scene_name = scenes{i};
        scene_data = gif_scenes.(scene_name);
        
        if ~isempty(scene_data.frames)
            fprintf('Creating %s...\n', scene_data.filename);
            
            % Create GIF
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

%% Predictor step: Solve momentum equations for intermediate velocities
function [u_star, v_star] = predictor_step(u, v, p, dx, dy, dt, nu, alpha)
n = size(u,1);
u_star = u; 
v_star = v;

% Pre-calculate constant terms for efficiency
inv_4dx = 1/(4*dx);
inv_4dy = 1/(4*dy);
inv_dx_sq = 1/dx^2;
inv_dy_sq = 1/dy^2;
alpha_dt = alpha * dt;

% Calculate intermediate u-velocity
for j = 2:n-1
    for i = 2:n-1
        % Convective terms (using divergence form)
        du2dx = ((u(j,i)+u(j,i+1))^2 - (u(j,i-1)+u(j,i))^2) * inv_4dx;
        duvdy = ((v(j,i)+v(j,i+1))*(u(j,i)+u(j+1,i)) - ...
                 (v(j-1,i)+v(j-1,i+1))*(u(j-1,i)+u(j,i))) * inv_4dy;
        
        % Diffusion terms
        d2udx2 = (u(j,i+1)-2*u(j,i)+u(j,i-1)) * inv_dx_sq;
        d2udy2 = (u(j+1,i)-2*u(j,i)+u(j-1,i)) * inv_dy_sq;
        
        % Pressure gradient
        dpdx = (p(j,i+1) - p(j,i)) / dx;
        
        % Update u_star with under-relaxation
        u_star(j,i) = u(j,i) + alpha_dt * (-du2dx - duvdy - dpdx + nu*(d2udx2 + d2udy2));
    end
end

% Calculate intermediate v-velocity
for j = 2:n-1
    for i = 2:n-1
        % Convective terms (using divergence form)
        dv2dy = ((v(j,i)+v(j+1,i))^2 - (v(j-1,i)+v(j,i))^2) * inv_4dy;
        duvdx = ((u(j+1,i)+u(j,i))*(v(j,i+1)+v(j,i)) - ...
                 (u(j+1,i-1)+u(j,i-1))*(v(j,i)+v(j,i-1))) * inv_4dx;
        
        % Diffusion terms
        d2vdx2 = (v(j,i+1)-2*v(j,i)+v(j,i-1)) * inv_dx_sq;
        d2vdy2 = (v(j+1,i)-2*v(j,i)+v(j-1,i)) * inv_dy_sq;
        
        % Pressure gradient
        dpdy = (p(j+1,i) - p(j,i)) / dy;
        
        % Update v_star with under-relaxation
        v_star(j,i) = v(j,i) + alpha_dt * (-duvdx - dv2dy - dpdy + nu*(d2vdx2 + d2vdy2));
    end
end
end

%% Pressure Poisson solver
function p_prime = solve_pressure_poisson(u_star, v_star, dx, dy, dt, tol, max_iter)
n = size(u_star,1);
p_prime = zeros(n);

% Pre-calculate constant terms for efficiency
inv_dx = 1/dx;
inv_dy = 1/dy;
dt_rhs_factor = 1/(dt);
laplacian_factor = 0.25;

for iter = 1:max_iter
    p_old = p_prime;
    
    % Iterate over interior points
    for j = 2:n-1
        for i = 2:n-1
            % Source term (divergence of intermediate velocity field)
            rhs = ( (u_star(j,i) - u_star(j,i-1))*inv_dx + (v_star(j,i) - v_star(j-1,i))*inv_dy ) * dt_rhs_factor;
            
            % Jacobi iteration for pressure correction
            p_prime(j,i) = laplacian_factor * (p_prime(j,i+1) + p_prime(j,i-1) + ...
                                   p_prime(j+1,i) + p_prime(j-1,i) - dx^2 * rhs);
        end
    end
    
    % Apply homogeneous Neumann boundary conditions
    p_prime(1,:) = p_prime(2,:);     % Bottom
    p_prime(end,:) = p_prime(end-1,:); % Top
    p_prime(:,1) = p_prime(:,2);     % Left
    p_prime(:,end) = p_prime(:,end-1); % Right
    
    % Check for convergence
    if max(max(abs(p_prime - p_old))) < tol
        break;
    end
end
end

%% Corrector step: Update velocities and pressure
function [u, v, p] = corrector_step(u_star, v_star, p, p_prime, dx, dy, dt, alpha)
n = size(p,1);
u = u_star; 
v = v_star;

% Pre-calculate constant terms for efficiency
alpha_dt_dx = alpha * dt / dx;
alpha_dt_dy = alpha * dt / dy;

% Velocity correction
for j = 2:n-1
    for i = 2:n-1
        u(j,i) = u_star(j,i) - alpha_dt_dx * (p_prime(j,i+1) - p_prime(j,i));
        v(j,i) = v_star(j,i) - alpha_dt_dy * (p_prime(j+1,i) - p_prime(j,i));
    end
end

% Pressure correction
p = p + alpha * p_prime;

% Apply boundary conditions
% Velocity BCs
u(:,1) = 0; u(:,end) = 0; u(1,:) = 0; u(end,:) = 1; % Moving lid
v(:,1) = 0; v(:,end) = 0; v(1,:) = 0; v(end,:) = 0; % No-slip on walls

% Pressure BCs (homogeneous Neumann)
p(1,:) = p(2,:);       % Bottom
p(end,:) = p(end-1,:); % Top
p(:,1) = p(:,2);       % Left
p(:,end) = p(:,end-1); % Right
end

%% Post-processing function
function plot_final_results(X, Y, u, v, p, res_u, res_v, res_p)
% Create a comprehensive figure with final results

figure('Name','Final Results - Lid Driven Cavity',...
       'Units','normalized','Position',[0.05 0.05 0.9 0.85], 'Color', 'w');

% Velocity vectors and streamlines
subplot(2,3,1);
quiver(X(1:3:end,1:3:end), Y(1:3:end,1:3:end),...
       u(1:3:end,1:3:end), v(1:3:end,1:3:end), 1.5, 'k');
hold on;
startx = linspace(0, 1, 20);
starty = linspace(0, 1, 20);
[sx, sy] = meshgrid(startx, starty);
streamline(X, Y, u, v, sx, sy);
title('Velocity Vectors and Streamlines');
xlabel('X'); ylabel('Y');
axis equal tight;
grid on;

% Velocity magnitude contour
subplot(2,3,2);
velMag = sqrt(u.^2 + v.^2);
if max(velMag(:)) - min(velMag(:)) > eps
    contourf(X, Y, velMag, 20, 'LineColor','none');
else
    imagesc(velMag);
    set(gca, 'YDir', 'normal');
end
colorbar;
title('Velocity Magnitude');
xlabel('X'); ylabel('Y');
axis equal tight;

% Pressure contour
subplot(2,3,3);
if max(p(:)) - min(p(:)) > eps
    contourf(X, Y, p, 20, 'LineColor','none');
else
    imagesc(p);
    set(gca, 'YDir', 'normal');
end
colorbar;
title('Pressure Field');
xlabel('X'); ylabel('Y');
axis equal tight;

% Vorticity contour
subplot(2,3,4);
vorticity = curl(X, Y, u, v);
if max(vorticity(:)) - min(vorticity(:)) > eps
    contourf(X, Y, vorticity, 20, 'LineColor','none');
else
    imagesc(vorticity);
    set(gca, 'YDir', 'normal');
end
colorbar;
title('Vorticity Field');
xlabel('X'); ylabel('Y');
axis equal tight;

% Centerline velocity profiles
subplot(2,3,5);
plot(u(ceil(size(u,1)/2),:), Y(ceil(size(Y,1)/2),:), 'b-', 'LineWidth', 2);
hold on;
plot(u(:,ceil(size(u,2)/2)), X(:,ceil(size(X,2)/2)), 'r-', 'LineWidth', 2);
title('Centerline Velocity Profiles');
xlabel('Velocity');
ylabel('Position');
legend('Vertical Centerline', 'Horizontal Centerline');
grid on;

% Residual history
subplot(2,3,6);
semilogy(1:length(res_u), res_u, 'r-', 'DisplayName', 'u-residual', 'LineWidth', 1.5);
hold on;
semilogy(1:length(res_v), res_v, 'g-', 'DisplayName', 'v-residual', 'LineWidth', 1.5);
semilogy(1:length(res_p), res_p, 'b-', 'DisplayName', 'p-residual', 'LineWidth', 1.5);
title('Convergence History');
xlabel('Time Step');
ylabel('Residual (log scale)');
legend('Location','northeast');
grid on;

% Add annotation with simulation parameters
annotation('textbox', [0.02, 0.02, 0.3, 0.05], 'String', ...
           sprintf('Re = %d, Grid = %dx%d', 100, size(u,1), size(u,1)), ...
           'FitBoxToText', 'on', 'BackgroundColor', 'white');

% Save final results figure
saveas(gcf, 'final_results.png');
end

%% Helper function to calculate vorticity
function vort = curl(X, Y, u, v)
% Calculate vorticity field from velocity components
[dudy, dudx] = gradient(u, Y(1,2)-Y(1,1), X(2,1)-X(1,1));
[dvdy, dvdx] = gradient(v, Y(1,2)-Y(1,1), X(2,1)-X(1,1));
vort = dvdx - dudy;
end
