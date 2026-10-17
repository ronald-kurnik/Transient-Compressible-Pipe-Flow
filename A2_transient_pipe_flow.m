%% 02_transient_pipe_flow.m
% 1-D unsteady compressible / acoustic pipe flow (method of characteristics style)
% + simple heat transfer wall model
%
% Models a short pipe with prescribed pressure at both ends.
% Uses 1-D Euler equations with isothermal wall heat loss.

clear; close all; clc;

%% Gas properties (air-like)
gamma = 1.4;
R     = 287;          % J/kg-K
cv    = R/(gamma-1);
cp    = gamma*cv;

%% Geometry
L     = 0.5;          % m  (characteristic length scale of generator duct)
D     = 0.05;         % m
A     = pi*(D/2)^2;
Pwall = pi*D;         % perimeter

%% Grid
Nx    = 100;
dx    = L/(Nx-1);
x     = linspace(0,L,Nx)';

%% Initial conditions (rest)
p0    = 1e5;          % Pa
T0    = 300;          % K
rho0  = p0/(R*T0);
u0    = zeros(Nx,1);
rho   = rho0 * ones(Nx,1);
u     = u0;
p     = p0 * ones(Nx,1);
T     = T0 * ones(Nx,1);
e     = cv*T;         % internal energy

%% Boundary forcing (simulate valve / piston-like pressure pulse)
t_pulse = 0.002;      % s
dp      = 0.2e5;      % Pa amplitude

%% Time integration (explicit, CFL-limited)
c0    = sqrt(gamma*R*T0);
CFL   = 0.5;
dt    = CFL * dx / (c0 + 50);   % safety
Tfinal= 0.010;
Nt    = ceil(Tfinal/dt);
dt    = Tfinal/Nt;

% Storage for animation / post
rho_hist = zeros(Nx,Nt+1);
u_hist   = zeros(Nx,Nt+1);
p_hist   = zeros(Nx,Nt+1);
rho_hist(:,1) = rho; u_hist(:,1) = u; p_hist(:,1) = p;

%% Simple wall heat transfer coefficient
h_wall = 50;          % W/m^2-K  (order-of-magnitude)

for n = 1:Nt
    t = n*dt;
    
    % Primitive -> conservative
    mom = rho.*u;
    E   = rho.*(e + 0.5*u.^2);   % total energy
    
    % Fluxes (central + artificial viscosity for stability demo)
    % Better: use HLLC or Rusanov; here simple Rusanov-like
    F_rho = rho.*u;
    F_mom = rho.*u.^2 + p;
    F_E   = u.*(E + p);
    
    % Rusanov numerical flux
    c     = sqrt(gamma*R*T);
    lambda= abs(u)+c;
    alpha = max(lambda);                 % global for simplicity
    
    F_rho_hat = zeros(Nx-1,1);
    F_mom_hat = zeros(Nx-1,1);
    F_E_hat   = zeros(Nx-1,1);
    for i=1:Nx-1
        amax = max(lambda(i),lambda(i+1));
        F_rho_hat(i) = 0.5*(F_rho(i)+F_rho(i+1)) - 0.5*amax*(rho(i+1)-rho(i));
        F_mom_hat(i) = 0.5*(F_mom(i)+F_mom(i+1)) - 0.5*amax*(mom(i+1)-mom(i));
        F_E_hat(i)   = 0.5*(F_E(i)+F_E(i+1))     - 0.5*amax*(E(i+1)-E(i));
    end
    
    % Update interior
    dF_rho = (F_rho_hat(2:end) - F_rho_hat(1:end-1))/dx;
    dF_mom = (F_mom_hat(2:end) - F_mom_hat(1:end-1))/dx;
    dF_E   = (F_E_hat(2:end)   - F_E_hat(1:end-1))/dx;
    
    rho(2:end-1) = rho(2:end-1) - dt * dF_rho;
    mom(2:end-1) = mom(2:end-1) - dt * dF_mom;
    E(2:end-1)   = E(2:end-1)   - dt * dF_E;
    
    % Wall heat loss (source term on energy)
    Qdot = -h_wall * Pwall/A * (T(2:end-1) - T0);   % W/m^3
    E(2:end-1) = E(2:end-1) + dt * Qdot;
    
    % Recover primitives
    u   = mom ./ rho;
    e   = E./rho - 0.5*u.^2;
    T   = e / cv;
    p   = rho * R .* T;
    
    % Boundary conditions (pressure prescribed, characteristic)
    % Left: p = p0 + dp*sin^2 pulse
    if t < t_pulse
        p_left = p0 + dp * sin(pi*t/t_pulse)^2;
    else
        p_left = p0;
    end
    p(1)   = p_left;
    T(1)   = T0;          % isothermal inlet for simplicity
    rho(1) = p(1)/(R*T(1));
    u(1)   = u(2);        % zero-gradient velocity (approx)
    
    % Right: fixed pressure
    p(end)   = p0;
    T(end)   = T0;
    rho(end) = p0/(R*T0);
    u(end)   = u(end-1);
    
    % Store
    rho_hist(:,n+1) = rho;
    u_hist(:,n+1)   = u;
    p_hist(:,n+1)   = p;
end

%% Visualization
figure('Position',[50 50 1000 700], 'Name','Transient Pipe Flow');
t_vec = (0:Nt)*dt;

subplot(2,2,1);
imagesc(t_vec*1000, x, p_hist/1e5);
colorbar; xlabel('t [ms]'); ylabel('x [m]');
title('Pressure [bar]'); set(gca,'YDir','normal');

subplot(2,2,2);
imagesc(t_vec*1000, x, u_hist);
colorbar; xlabel('t [ms]'); ylabel('x [m]');
title('Velocity [m/s]'); set(gca,'YDir','normal');

subplot(2,2,3);
plot(x, p_hist(:,1)/1e5, 'k--', x, p_hist(:,round(Nt/4))/1e5, 'b', ...
     x, p_hist(:,round(Nt/2))/1e5, 'r', x, p_hist(:,end)/1e5, 'g', 'LineWidth',1.3);
legend('t=0','t=T/4','t=T/2','t=T','Location','best');
xlabel('x [m]'); ylabel('p [bar]'); title('Pressure profiles'); grid on;

subplot(2,2,4);
plot(x, u_hist(:,round(Nt/4)), 'b', x, u_hist(:,round(Nt/2)), 'r', ...
     x, u_hist(:,end), 'g', 'LineWidth',1.3);
xlabel('x [m]'); ylabel('u [m/s]'); title('Velocity profiles'); grid on;

fprintf('\n=== Transient pipe flow demo complete ===\n');
fprintf('This is a minimal 1-D Euler solver with wall heat transfer.\n');
fprintf('In production you would use a proper Riemann solver, adaptive mesh,\n');
fprintf('and couple to detailed chemical kinetics or system-level models.\n');
