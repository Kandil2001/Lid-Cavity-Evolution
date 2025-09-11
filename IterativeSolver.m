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
%   - Real-time visualization of convergence
%
% OUTPUTS:
%   Elapsed real time is printed to the console, providing the benchmark metric.
%   Figures for velocity, pressure, and residuals are animated during the run.
%   GIF files are saved for each visualization scene.
%
% BENCHMARK: On a 151x151 grid, Re=100, this implementation serves as a
%            performance baseline for comparison with optimized versions.
%
% -------------------------------------------------------------------------

%% Simulation parameters
clear; clc; close all;

% Physical parameters
Re = 100;               % Reynolds number
L = 1.0;                % Cavity length

% Numerical parameters
n = 151;                % Grid size (n x n) - 151x151 provides good resolution
dx = L/(n-1); 
dy = dx;                % Square cells
dt = 0.0005;            % Time step - stability constraint for explicit methods
nu = 1/Re;              % Kinematic viscosity
alpha_u = 0.5;          % Under-relaxation for velocity (0 < alpha_u <= 1)
alpha_p = 0.2;          % Under-relaxation for pressure (0 < alpha_p <= 1)
tol = 1e-5;             % SIMPLE convergence tolerance
max_iter = 500;         % Max SIMPLE iterations per time step
total_time = 2;         % Total simulation time

% GIF recording parameters
record_gif = true;      % Set to true to record GIFs
gif_frame_interval = 5; % Save every Nth frame to GIF

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
max_steps = ceil(total_time/dt);

% Pre-allocate arrays for residual history with reasonable initial size
res_u_arr = zeros(1000,1); 
res_v_arr = zeros(1000,1); 
res_p_arr = zeros(1000,1);

% Initialize GIF structures if recording
if record_gif
    gif_scenes = struct();
    gif_scenes.velocity_vectors = struct('frames', [], 'filename', 'velocity_vectors.gif');
    gif_scenes.velocity_contour = struct('frames', [], 'filename', 'velocity_contour.gif');
    gif_scenes.pressure_contour = struct('frames', [], 'filename', 'pressure_contour.gif');
    gif_scenes.residuals = struct('frames', [], 'filename', 'residuals_convergence.gif');
    gif_scenes.streamlines = struct('frames', [], 'filename', 'streamlines.gif');
end

%% Figure setup for real-time visualization
hFig = figure('Name','SIMPLE Lid Driven Cavity','Units','normalized',...
             'Position',[0.05 0.1 0.9 0.8], 'Color', 'w');

% Velocity vectors subplot
subplot(2,3,1);
h_quiver = quiver(X(1:4:end,1:4:end), Y(1:4:end,1:4:end), ...
                 u(1:4:end,1:4:end), v(1:4:end,1:4:end), 2, 'k');
title('Velocity Vectors'); 
xlabel('X'); ylabel('Y'); 
axis equal tight;
grid on;

% Velocity magnitude contour subplot
subplot(2,3,2);
velMag = sqrt(u.^2 + v.^2);
% Avoid constant ZData warning by checking if data varies
if max(velMag(:)) - min(velMag(:)) > eps
    [~, h_contour_vel] = contourf(X, Y, velMag, 20, 'LineColor','none');
else
    % If data is constant, create a uniform plot
    imagesc(velMag);
    set(gca, 'YDir', 'normal');
end
colorbar;
title('Velocity Magnitude Contour'); 
xlabel('X'); ylabel('Y'); 
axis equal tight;

% Pressure contour subplot
subplot(2,3,3);
% Avoid constant ZData warning by checking if data varies
if max(p(:)) - min(p(:)) > eps
    [~, h_contour_p] = contourf(X, Y, p, 20, 'LineColor','none');
else
    % If data is constant, create a uniform plot
    imagesc(p);
    set(gca, 'YDir', 'normal');
end
colorbar;
title('Pressure Contour'); 
xlabel('X'); ylabel('Y'); 
axis equal tight;

% Streamlines subplot
subplot(2,3,4);
startx = linspace(0, L, 15);
starty = linspace(0, L, 15);
[sx, sy] = meshgrid(startx, starty);
h_stream = streamline(X, Y, u, v, sx, sy);
set(h_stream, 'Color', 'b', 'LineWidth', 1);
title('Streamlines');
xlabel('X'); ylabel('Y');
axis equal tight;
grid on;

% Residuals subplot
subplot(2,3,5);
h_res_u = semilogy(1, res_u_arr(1), '-r', 'DisplayName', 'u-res', 'LineWidth', 1.5);
hold on;
h_res_v = semilogy(1, res_v_arr(1), '-g', 'DisplayName', 'v-res', 'LineWidth', 1.5);
h_res_p = semilogy(1, res_p_arr(1), '-b', 'DisplayName', 'p-res', 'LineWidth', 1.5);
hold off;
xlabel('Time Step'); 
ylabel('Residual (log scale)');
title('Convergence History');
legend('Location','northeast');
xlim([1 max_steps]); 
ylim([1e-8 1]); 
grid on;

% Info subplot
subplot(2,3,6);
h_info = text(0.1, 0.9, sprintf('Time: %.3f s\nStep: %d\nRe: %d', time, step, Re), ...
             'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold');
axis off;
title('Simulation Info');

%% Simulation execution
fprintf('Starting SIMPLE Lid Driven Cavity Simulation...\n');
fprintf('Grid size: %dx%d, Reynolds number: %d\n', n, n, Re);
fprintf('SIMPLE tolerance: %g, Max iterations: %d\n', tol, max_iter);
fprintf('Starting time integration. This may take a while...\n');

tic;  % Start measuring real time

%% Time stepping loop
while time < total_time
    step = step + 1;
    
    % Check if we need to expand residual storage
    if step > length(res_u_arr)
        % Double the storage to avoid growing arrays every step
        new_size = 2 * length(res_u_arr);
        res_u_arr(new_size) = 0;
        res_v_arr(new_size) = 0;
        res_p_arr(new_size) = 0;
    end
    
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

    % --- Real-Time Animation: Update every time step ---
    % Update velocity vectors
    subplot(2,3,1);
    set(h_quiver, 'UData', u(1:4:end,1:4:end), 'VData', v(1:4:end,1:4:end));
    title(sprintf('Velocity Vectors (Time = %.3f s)', time));
    
    % Update velocity magnitude contour
    subplot(2,3,2);
    velMag = sqrt(u.^2 + v.^2);
    cla;
    if max(velMag(:)) - min(velMag(:)) > eps
        contourf(X, Y, velMag, 20, 'LineColor','none');
    else
        imagesc(velMag);
        set(gca, 'YDir', 'normal');
    end
    colorbar;
    title(sprintf('Velocity Magnitude (Max = %.2f)', max(velMag(:))));
    
    % Update pressure contour
    subplot(2,3,3);
    cla;
    if max(p(:)) - min(p(:)) > eps
        contourf(X, Y, p, 20, 'LineColor','none');
    else
        imagesc(p);
        set(gca, 'YDir', 'normal');
    end
    colorbar;
    title(sprintf('Pressure Field (Max = %.2f)', max(p(:))));
    
    % Update streamlines
    subplot(2,3,4);
    delete(findobj(gca, 'Type', 'line'));
    h_stream = streamline(X, Y, u, v, sx, sy);
    set(h_stream, 'Color', 'b', 'LineWidth', 1);
    title('Streamlines');
    
    % Update residuals plot
    subplot(2,3,5);
    set(h_res_u, 'XData', 1:step, 'YData', res_u_arr(1:step));
    set(h_res_v, 'XData', 1:step, 'YData', res_v_arr(1:step));
    set(h_res_p, 'XData', 1:step, 'YData', res_p_arr(1:step));
    title(sprintf('Residuals (Iteration %d)', step));
    xlim([1 max(step, 10)]); % Adjust x-axis as simulation progresses
    
    % Update info panel
    subplot(2,3,6);
    cla;
    text(0.1, 0.9, sprintf('Time: %.3f s\nStep: %d\nRe: %d\nResiduals:\n  u: %.2e\n  v: %.2e\n  p: %.2e', ...
          time, step, Re, res_u, res_v, res_p), ...
         'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold');
    axis off;
    title('Simulation Info');
    
    % Capture frames for GIFs
    if record_gif && mod(step, gif_frame_interval) == 0
        % Velocity vectors scene
        subplot(2,3,1);
        frame = getframe(gcf);
        gif_scenes.velocity_vectors.frames = [gif_scenes.velocity_vectors.frames, frame];
        
        % Velocity contour scene
        subplot(2,3,2);
        frame = getframe(gcf);
        gif_scenes.velocity_contour.frames = [gif_scenes.velocity_contour.frames, frame];
        
        % Pressure contour scene
        subplot(2,3,3);
        frame = getframe(gcf);
        gif_scenes.pressure_contour.frames = [gif_scenes.pressure_contour.frames, frame];
        
        % Residuals scene
        subplot(2,3,5);
        frame = getframe(gcf);
        gif_scenes.residuals.frames = [gif_scenes.residuals.frames, frame];
        
        % Streamlines scene
        subplot(2,3,4);
        frame = getframe(gcf);
        gif_scenes.streamlines.frames = [gif_scenes.streamlines.frames, frame];
    end
    
    % Force figure update
    drawnow;
    
    % Display progress every 10 steps
    if mod(step, 10) == 0
        fprintf('Step: %d, Simulation time: %.4f/%.2f, Residuals: u=%.2e, v=%.2e, p=%.2e\n', ...
                step, time, total_time, res_u, res_v, res_p);
    end

    % Advance time
    time = time + dt;
end

%% Create GIFs from captured frames
if record_gif
    fprintf('Creating GIF files for GitHub...\n');
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

% Calculate intermediate u-velocity
for j = 2:n-1
    for i = 2:n-1
        % Convective terms (using divergence form)
        du2dx = ((u(j,i)+u(j,i+1))^2 - (u(j,i-1)+u(j,i))^2)/(4*dx);
        duvdy = ((v(j,i)+v(j,i+1))*(u(j,i)+u(j+1,i)) - ...
                 (v(j-1,i)+v(j-1,i+1))*(u(j-1,i)+u(j,i)))/(4*dy);
        
        % Diffusion terms
        d2udx2 = (u(j,i+1)-2*u(j,i)+u(j,i-1))/dx^2;
        d2udy2 = (u(j+1,i)-2*u(j,i)+u(j-1,i))/dy^2;
        
        % Pressure gradient
        dpdx = (p(j,i+1) - p(j,i))/dx;
        
        % Update u_star with under-relaxation
        u_star(j,i) = u(j,i) + alpha * dt * (-du2dx - duvdy - dpdx + nu*(d2udx2 + d2udy2));
    end
end

% Calculate intermediate v-velocity
for j = 2:n-1
    for i = 2:n-1
        % Convective terms (using divergence form)
        dv2dy = ((v(j,i)+v(j+1,i))^2 - (v(j-1,i)+v(j,i))^2)/(4*dy);
        duvdx = ((u(j+1,i)+u(j,i))*(v(j,i+1)+v(j,i)) - ...
                 (u(j+1,i-1)+u(j,i-1))*(v(j,i)+v(j,i-1)))/(4*dx);
        
        % Diffusion terms
        d2vdx2 = (v(j,i+1)-2*v(j,i)+v(j,i-1))/dx^2;
        d2vdy2 = (v(j+1,i)-2*v(j,i)+v(j-1,i))/dy^2;
        
        % Pressure gradient
        dpdy = (p(j+1,i) - p(j,i))/dy;
        
        % Update v_star with under-relaxation
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
    
    % Iterate over interior points
    for j = 2:n-1
        for i = 2:n-1
            % Source term (divergence of intermediate velocity field)
            rhs = ((u_star(j,i) - u_star(j,i-1))/dx + (v_star(j,i) - v_star(j-1,i))/dy)/dt;
            
            % Jacobi iteration for pressure correction
            p_prime(j,i) = 0.25 * (p_prime(j,i+1) + p_prime(j,i-1) + ...
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

% Velocity correction
for j = 2:n-1
    for i = 2:n-1
        u(j,i) = u_star(j,i) - alpha * dt * (p_prime(j,i+1) - p_prime(j,i)) / dx;
        v(j,i) = v_star(j,i) - alpha * dt * (p_prime(j+1,i) - p_prime(j,i)) / dy;
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
