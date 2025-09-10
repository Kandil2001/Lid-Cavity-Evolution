function SIMPLE2D_LidDrivenCavity_Optimized()
% SIMPLE2D Lid-Driven Cavity Simulation (Performance-Optimized, No Vectorization)
% -------------------------------------------------------------------------
% - Finite Volume, Staggered Grid, SIMPLE Algorithm
% - Optimized triple-nested loops, precomputed constants, minimal memory overhead
% - Real-time visualization, GIF saving (optional)
% -------------------------------------------------------------------------
clear; clc; close all;

%% Parameters
Re = 100;               % Reynolds number
L = 1.0;                % Cavity length (square)
n = 151;                % Grid size
dx = L/(n-1); dy = dx;  % Uniform grid
dt = 0.001;            % Time step
nu = 1/Re;              % Kinematic viscosity
alpha_u = 0.5;          % Underrelaxation (velocity)
alpha_p = 0.2;          % Underrelaxation (pressure)
tol = 1e-5;             % SIMPLE tolerance
max_iter = 500;         % Max SIMPLE iterations/step
total_time = 2;         % Total simulation time
max_steps = ceil(total_time/dt);

% GIF recording
record_gif = true; gif_frame_interval = 5;

%% Grid + Initialization
[X, Y] = meshgrid(0:dx:L, 0:dy:L);
u = zeros(n); v = zeros(n); p = zeros(n);
u(end,:) = 1; % Lid BC

time = 0; step = 0;
res_u_arr = zeros(1000,1); res_v_arr = zeros(1000,1); res_p_arr = zeros(1000,1);

if record_gif
    gif_scenes = struct();
    scnames = {'velocity_vectors','velocity_contour','pressure_contour','residuals','streamlines'};
    for i = 1:length(scnames)
        gif_scenes.(scnames{i}) = struct('frames',[],'filename',[scnames{i},'.gif']);
    end
end

%% Figure Setup
hFig = figure('Name','SIMPLE Lid Driven Cavity','Units','normalized','Position',[0.05 0.1 0.9 0.8],'Color','w');
subplot(2,3,1); h_quiver = quiver(X(1:4:end,1:4:end), Y(1:4:end,1:4:end), u(1:4:end,1:4:end), v(1:4:end,1:4:end),2,'k');
title('Velocity Vectors'); xlabel('X'); ylabel('Y'); axis equal tight; grid on;
subplot(2,3,2); contourf(X,Y,sqrt(u.^2+v.^2),20,'LineColor','none'); colorbar; title('Velocity Magnitude Contour'); axis equal tight;
subplot(2,3,3); contourf(X,Y,p,20,'LineColor','none'); colorbar; title('Pressure Contour'); axis equal tight;
subplot(2,3,4); [sx,sy]=meshgrid(linspace(0,L,15),linspace(0,L,15)); h_stream = streamline(X,Y,u,v,sx,sy); set(h_stream,'Color','b','LineWidth',1); title('Streamlines'); axis equal tight; grid on;
subplot(2,3,5); h_res_u=semilogy(1,1e-1,'-r','LineWidth',1.5); hold on; h_res_v=semilogy(1,1e-1,'-g','LineWidth',1.5); h_res_p=semilogy(1,1e-1,'-b','LineWidth',1.5); hold off; xlabel('Time Step'); ylabel('Residual'); title('Convergence History'); legend('u','v','p'); grid on; xlim([1 max_steps]); ylim([1e-8 1]);
subplot(2,3,6); h_info = text(0.1,0.9,sprintf('Time: %.3f s\nStep: %d\nRe: %d',time,step,Re),'Units','normalized','FontSize',12,'FontWeight','bold'); axis off; title('Simulation Info');

%% Main Time-Stepping Loop
fprintf('Starting SIMPLE Lid Driven Cavity Simulation...\n');
tic;
while time < total_time
    step = step + 1;
    if step > numel(res_u_arr)
        res_u_arr(2*step) = 0; res_v_arr(2*step) = 0; res_p_arr(2*step) = 0;
    end

    u_old = u; v_old = v; p_old = p;

    % SIMPLE Inner Iteration
    for iter = 1:max_iter
        u_prev = u; v_prev = v; p_prev = p;
        [u_star, v_star] = predictor_step_fast(u, v, p, dx, dy, dt, nu, alpha_u);
        p_prime = solve_pressure_poisson_fast(u_star, v_star, dx, dy, dt, tol, max_iter);
        [u, v, p] = corrector_step_fast(u_star, v_star, p, p_prime, dx, dy, dt, alpha_p);

        % Residuals
        res_u = max(max(abs(u-u_prev)));
        res_v = max(max(abs(v-v_prev)));
        res_p = max(max(abs(p-p_prev)));
        if max([res_u,res_v,res_p]) < tol, break; end
    end
    res_u_arr(step) = res_u; res_v_arr(step) = res_v; res_p_arr(step) = res_p;

    % --- Visualization ---
    % (Update only every N steps for speed if desired)
    subplot(2,3,1); set(h_quiver,'UData',u(1:4:end,1:4:end),'VData',v(1:4:end,1:4:end)); title(sprintf('Velocity Vectors (t=%.3f)',time));
    subplot(2,3,2); cla; contourf(X,Y,sqrt(u.^2+v.^2),20,'LineColor','none'); colorbar; title(sprintf('Velocity Magnitude (max=%.2f)',max(sqrt(u(:).^2+v(:).^2))));
    subplot(2,3,3); cla; contourf(X,Y,p,20,'LineColor','none'); colorbar; title(sprintf('Pressure Field (max=%.2f)',max(p(:))));
    subplot(2,3,4); delete(findobj(gca,'Type','line')); h_stream = streamline(X,Y,u,v,sx,sy); set(h_stream,'Color','b','LineWidth',1); title('Streamlines');
    subplot(2,3,5); set(h_res_u,'XData',1:step,'YData',res_u_arr(1:step)); set(h_res_v,'XData',1:step,'YData',res_v_arr(1:step)); set(h_res_p,'XData',1:step,'YData',res_p_arr(1:step)); title(sprintf('Residuals (Step %d)',step)); xlim([1 max(step,10)]);
    subplot(2,3,6); cla; text(0.1,0.9,sprintf('Time: %.3f s\nStep: %d\nRe: %d\nResiduals:\n  u: %.2e\n  v: %.2e\n  p: %.2e',time,step,Re,res_u,res_v,res_p),'Units','normalized','FontSize',10,'FontWeight','bold'); axis off; title('Simulation Info');
    drawnow;

    % GIF frame capture
    if record_gif && mod(step,gif_frame_interval)==0
        subplot(2,3,1); gif_scenes.velocity_vectors.frames = [gif_scenes.velocity_vectors.frames, getframe(gcf)];
        subplot(2,3,2); gif_scenes.velocity_contour.frames = [gif_scenes.velocity_contour.frames, getframe(gcf)];
        subplot(2,3,3); gif_scenes.pressure_contour.frames = [gif_scenes.pressure_contour.frames, getframe(gcf)];
        subplot(2,3,5); gif_scenes.residuals.frames = [gif_scenes.residuals.frames, getframe(gcf)];
        subplot(2,3,4); gif_scenes.streamlines.frames = [gif_scenes.streamlines.frames, getframe(gcf)];
    end

    if mod(step,25)==0
        fprintf('Step: %d, t=%.3f/%.2f, Residuals: u=%e, v=%e, p=%e\n',step,time,total_time,res_u,res_v,res_p);
    end
    time = time + dt;
end
elapsedTime = toc;
fprintf('\nSimulation complete. Elapsed time: %.2f seconds.\n',elapsedTime);

% Save GIFs
if record_gif, create_gifs(gif_scenes); end

% Final Post-Processing
plot_final_results_fast(X,Y,u,v,p,res_u_arr(1:step),res_v_arr(1:step),res_p_arr(1:step));

end

% ------------------- FAST LOOP FUNCTIONS (NO VECTORIZATION) --------------------

function [u_star, v_star] = predictor_step_fast(u, v, p, dx, dy, dt, nu, alpha)
n = size(u,1);
u_star = u; v_star = v;
inv4dx = 0.25/dx; inv4dy = 0.25/dy; inv_dx2 = 1/(dx*dx); inv_dy2 = 1/(dy*dy); dt_alpha = dt * alpha;

for j=2:n-1
    for i=2:n-1
        uji=u(j,i); ujip1=u(j,i+1); ujim1=u(j,i-1); ujp1i=u(j+1,i); ujm1i=u(j-1,i);
        vji=v(j,i); vjip1=v(j,i+1); vjm1i=v(j-1,i); vjm1ip1=v(j-1,i+1); vjp1i=v(j+1,i);

        du2dx = ((uji+ujip1)^2 - (ujim1+uji)^2) * inv4dx;
        duvdy = ((vji+vjip1)*(uji+ujp1i) - (vjm1i+vjm1ip1)*(ujm1i+uji)) * inv4dy;
        d2udx2 = (ujip1-2*uji+ujim1)*inv_dx2; d2udy2 = (ujp1i-2*uji+ujm1i)*inv_dy2;
        dpdx = (p(j,i+1)-p(j,i))/dx;

        u_star(j,i) = uji + dt_alpha*(-du2dx-duvdy-dpdx + nu*(d2udx2+d2udy2));
    end
end

for j=2:n-1
    for i=2:n-1
        vji=v(j,i); vjp1i=v(j+1,i); vjm1i=v(j-1,i);
        uji=u(j,i); ujp1i=u(j+1,i); ujim1=u(j,i-1); ujp1im1=u(j+1,i-1);
        vjip1=v(j,i+1); vjim1=v(j,i-1);
        dv2dy = ((vji+vjp1i)^2 - (vjm1i+vji)^2)*inv4dy;
        duvdx = ((ujp1i+uji)*(vjip1+vji) - (ujp1im1+ujim1)*(vji+vjim1))*inv4dx;
        d2vdx2 = (vjip1-2*vji+vjim1)*inv_dx2; d2vdy2 = (vjp1i-2*vji+vjm1i)*inv_dy2;
        dpdy = (p(j+1,i)-p(j,i))/dy;

        v_star(j,i) = vji + dt_alpha*(-duvdx-dv2dy-dpdy + nu*(d2vdx2+d2vdy2));
    end
end
end

function p_prime = solve_pressure_poisson_fast(u_star, v_star, dx, dy, dt, tol, max_iter)
n = size(u_star,1); p_prime = zeros(n);
fact = 0.25; dx2=dx*dx; dy2=dy*dy; invdt=1.0/dt;
for iter=1:max_iter
    p_old = p_prime;
    for j=2:n-1
        for i=2:n-1
            rhs = ((u_star(j,i)-u_star(j,i-1))/dx + (v_star(j,i)-v_star(j-1,i))/dy)*invdt;
            p_prime(j,i) = fact*(p_prime(j,i+1)+p_prime(j,i-1)+p_prime(j+1,i)+p_prime(j-1,i)-dx2*rhs);
        end
    end
    for i=1:n, p_prime(1,i)=p_prime(2,i); p_prime(n,i)=p_prime(n-1,i); p_prime(i,1)=p_prime(i,2); p_prime(i,n)=p_prime(i,n-1); end
    if max(max(abs(p_prime-p_old)))<tol, break; end
end
end

function [u, v, p] = corrector_step_fast(u_star, v_star, p, p_prime, dx, dy, dt, alpha)
n = size(p,1); u=u_star; v=v_star; dt_alpha=dt*alpha;
for j=2:n-1
    for i=2:n-1
        u(j,i) = u_star(j,i) - dt_alpha*(p_prime(j,i+1)-p_prime(j,i))/dx;
        v(j,i) = v_star(j,i) - dt_alpha*(p_prime(j+1,i)-p_prime(j,i))/dy;
    end
end
p = p + alpha*p_prime;
for i=1:n, u(i,1)=0; u(i,n)=0; v(i,1)=0; v(i,n)=0; end
for j=1:n, u(1,j)=0; u(n,j)=1; v(1,j)=0; v(n,j)=0; end
for i=1:n, p(1,i)=p(2,i); p(n,i)=p(n-1,i); p(i,1)=p(i,2); p(i,n)=p(i,n-1); end
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

function plot_final_results_fast(X,Y,u,v,p,res_u,res_v,res_p)
figure('Name','Final Results','Units','normalized','Position',[0.05 0.05 0.9 0.85],'Color','w');
subplot(2,3,1); quiver(X(1:3:end,1:3:end),Y(1:3:end,1:3:end),u(1:3:end,1:3:end),v(1:3:end,1:3:end),1.5,'k'); hold on;
[sx,sy]=meshgrid(linspace(0,1,20),linspace(0,1,20)); streamline(X,Y,u,v,sx,sy); title('Velocity Vectors & Streamlines'); axis equal tight; grid on;
subplot(2,3,2); contourf(X,Y,sqrt(u.^2+v.^2),20,'LineColor','none'); colorbar; title('Velocity Magnitude'); axis equal tight;
subplot(2,3,3); contourf(X,Y,p,20,'LineColor','none'); colorbar; title('Pressure Field'); axis equal tight;
subplot(2,3,4); vorticity = vorticity_field(X,Y,u,v); contourf(X,Y,vorticity,20,'LineColor','none'); colorbar; title('Vorticity Field'); axis equal tight;
subplot(2,3,5); plot(u(ceil(end/2),:),Y(ceil(end/2),:),'b-','LineWidth',2); hold on; plot(u(:,ceil(end/2)),X(:,ceil(end/2)),'r-','LineWidth',2); title('Centerline Velocity'); legend('Vertical','Horizontal'); grid on;
subplot(2,3,6); semilogy(1:length(res_u),res_u,'r-',1:length(res_v),res_v,'g-',1:length(res_p),res_p,'b-','LineWidth',1.5); title('Residual History'); legend('u','v','p'); grid on;
annotation('textbox',[0.02,0.02,0.3,0.05],'String',sprintf('Re=%d, Grid=%dx%d',100,size(u,1),size(u,2)),'FitBoxToText','on','BackgroundColor','white');
saveas(gcf,'final_results.png');
end

function vort = vorticity_field(X,Y,u,v)
[dudy, dudx] = gradient(u, Y(1,2)-Y(1,1), X(2,1)-X(1,1));
[dvdy, dvdx] = gradient(v, Y(1,2)-Y(1,1), X(2,1)-X(1,1));
vort = dvdx - dudy;
end
