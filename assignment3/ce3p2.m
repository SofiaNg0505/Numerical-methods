%%ce3p2
%%solve the 1D heat equation with semi-discretization and explicit euler

set(0,'defaultTextInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');


%% (a)
% solve from tau in [0, 2]
n = 99;
dx = 1/(n+1);
u0 = zeros(n+1, 1);

%% unstable

dt = 2.1*dx*dx;
tauMax = 2;
m = ceil(tauMax/dt);
dt = tauMax/m;

fprintf('dx = %f, dt = %f, dt/dx^2 = %f\n', dx, dt, dt/dx/dx);

a = 1/4;

% get the matrix A and vector b for the system with parameter a and n points
[A, b, x] = mol(a, dx, n);

pulse = @(t) sin(pi*t)*unitstep(1-t);

f = @(t, u) A*u + b*pulse(t);

% solve with explicit euler
[t, u] = euler(f, dt, m, u0);

% 3D plot
u = [arrayfun(pulse, 0:dt:tauMax); u'];
[X, T] = meshgrid(x, t);

fig3d = figure();
mesh(T, X, u');
xlabel('$\tau$');
ylabel('$x$');
zlabel('$u$');
title("Unstable solution");

% exportgraphics(fig3d, "ce3p2_unstable_3d.pdf");


%% stable
dt = 1.9*dx*dx;
tauMax = 2;
m = 2*ceil(tauMax/dt/2);
dt = tauMax/m;

fprintf('dx = %f, dt = %f, dt/dx^2 = %f\n', dx, dt, dt/dx/dx);

a = 1/4;

% get the matrix A and vector b for the system with parameter a and n points
[A, b, x] = mol(a, dx, n);

pulse = @(t) sin(pi*t*unitstep(1-t));

f = @(t, u) A*u + b*pulse(t);

% solve with explicit euler
[t, u] = euler(f, dt, m, u0);

% 3D plot
u = [arrayfun(pulse, 0:dt:tauMax); u'];
[X, T] = meshgrid(x, t);

fig3d = figure();
mesh(T, X, u');
xlabel('$\tau$');
ylabel('$x$');
zlabel('$u$');
title("Stable solution");

% exportgraphics(fig3d, "ce3p2_stable_3d.pdf");


%%
fig = figure();
% plot u(x, 1)
plot(x, u(:, m/2));
xlabel('$x$');
ylabel('$u$');
title("Spacial temperature profile at $\tau = 1$");
% exportgraphics(fig, "ce3p2_spacial.pdf");

fig = figure();
% plot u(0, tau)
plot(t, u(1, :));
hold on;
% plot u(1, tau)
plot(t, u(end, :));
xlabel('$\tau$');
ylabel('$u$');
legend("$x = 0$", "$x = 1$");
title("Temperature profile over time");
% exportgraphics(fig, "ce3p2_time.pdf");

%% (b)

% dt = 5*dx*dx;
tauMax = 2;
% m = ceil(tauMax/dt);
% dt = tauMax/m;
tspan = [0 tauMax];

for n = [99, 199, 399]
    dx = 1/(n+1);
    u0 = zeros(n+1, 1);
    
    [A, b, x] = mol(a, dx, n);
    f = @(t, u) A*u + b*pulse(t);
    
    opt = odeset('RelTol', 1e-3, 'AbsTol', 1e-6);
    
    %ode23
    tic; [t, u] = ode23(f, tspan, u0, opt); toc
    fprintf('ode23: N = %d, steps = %d\n', n+1, size(t, 1));
    
    %ode23s
    tic; [t, u] = ode23s(f, tspan, u0, opt); toc
    fprintf('ode23s: N = %d, steps = %d\n', n+1, size(t, 1));
    
    %ode23sJ
    opt = odeset('RelTol', 1e-3, 'AbsTol', 1e-3, 'Jacobian', A);
    tic; [t, u] = ode23s(f, tspan, u0, opt); toc
    fprintf('ode23sJ: N = %d, steps = %d\n', n+1, size(t, 1));
end

