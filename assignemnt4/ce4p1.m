%% Hyperbolic PDEs
%% Part 1: Model problem

% given is the advection equation u_t + a u_x = 0
% for a > 0, x in [0, D] and t in [0, Tmax]

%% constants
tau = 2.2;
D = 4.5;
a = 1;
Tmax = 6;

% gsin(t, tau)
% gsq(t, tau)
% are defined externally

%% (a) 
% compare Upwind, Lax-Friedrich and Lax-Wendroff

% number of grid points right of the zero
N = 99;
dx = D/N;

cfl = 0.9;
dt = cfl*dx;

% round to integer and redefine dt, cfl
t = 0:dt:Tmax;
% M = round(Tmax/dt);
% dt = Tmax/M;
% cfl = dt/dx;

fprintf('CFL: %f, dt = %f, dx = %f, N = %d, M = %d, T = %f\n', cfl, dt, dx, N, length(t), t(end));

% initial condition
u0 = zeros(1, N+1);

% spacial grid
x = 0:dx:D;
    
%% plotting
linewidth = 2;

figure();
subplot(1, 2, 1);

[uuw, ~] = upwind(cfl, a, dt, Tmax, N, u0, @(t) gsin(t, tau));
[ulxf, ~] = lxf(cfl, a, dt, Tmax, N, u0, @(t) gsin(t, tau));
[ulw, ~] = lw(cfl, a, dt, Tmax, N, u0, @(t) gsin(t, tau));

% plot at Tmax
plot(x, gsin(Tmax - x, tau), 'linewidth', linewidth);
hold on;
plot(x, ulxf(end, :), 'linewidth', linewidth);
plot(x, uuw(end, :), 'linewidth', linewidth);
plot(x, ulw(end, :), 'linewidth', linewidth);
xlabel('x');
ylabel('u');
legend('exact', 'LxF', 'Upwind', 'LW');
title('Sine wave');
grid on;

subplot(1, 2, 2);
[uuw, ~] = upwind(cfl, a, dt, Tmax, N, u0, @(t) gsq(t, tau));
[ulxf, ~] = lxf(cfl, a, dt, Tmax, N, u0, @(t) gsq(t, tau));
[ulw, ~] = lw(cfl, a, dt, Tmax, N, u0, @(t) gsq(t, tau));

% plot for at Tmax
plot(x, gsq(Tmax - x, tau), 'linewidth', linewidth);
hold on;
plot(x, ulxf(end, :), 'linewidth', linewidth);
plot(x, uuw(end, :), 'linewidth', linewidth);
plot(x, ulw(end, :), 'linewidth', linewidth);
xlabel('x');
ylabel('u');
legend('exact', 'LxF', 'Upwind', 'LW');
title('Square wave');
grid on;


% %% 3D plots
% [X, T] = meshgrid(x, t);

% % upwind
% [uuw, t] = upwind(cfl, a, dt, Tmax, N, u0, @(t) gsin(t, tau));
% % LxF
% [ulxf, ~] = lxf(cfl, a, dt, Tmax, N, u0, @(t) gsin(t, tau));
% % LW
% [ulw, ~] = lw(cfl, a, dt, Tmax, N, u0, @(t) gsin(t, tau));

% figure();
% mesh(T, X, ulxf);
% xlabel('t');
% ylabel('x');
% title('Lax-Friedrich');

% figure();
% mesh(T, X, uuw);
% xlabel('t');
% ylabel('x');
% title('Upwind');

% figure();
% mesh(T, X, ulw);
% xlabel('t');
% ylabel('x');
% title('Lax-Wendroff');
