%% Constants

Tcool = 50;
Thot = 200;
alpha = 0.5;
L = 5;
a = 1;
v = 1;
Tmax = 6;

% spacial grid
N = 100;
dx = L/N;


cfl = 0.9;
dt = cfl * dx;

M = round(Tmax/dt);
dt = Tmax / M;
cfl = dt/dx;

%% 3d plots
figure();
T0 = zeros(1, N+1) + Tcool;

x = 0:dx:L;
subplot(1,2,1);
[uuw, t] = upwind_2(cfl, a, dt, Tmax, N, M, @(t) T_bc(t, Tcool, Thot), alpha, Tcool, Thot,  T0);
[X,T] = meshgrid(x,t);
mesh(X,T,uuw, 'facecolor', 'flat');
xlabel('x');
ylabel('t');
zlabel('T');
title('Upwind');

subplot(1,2,2);
[lw, t] = lw_2(cfl, a, dt, Tmax, N, M, @(t) T_bc(t, Tcool, Thot), alpha, Tcool, Thot,  T0 , dx);
[X,T] = meshgrid(x,t);
mesh(X,T,lw, 'facecolor', 'flat');
xlabel('x');
ylabel('t');
zlabel('T');
title('Lax-Wendroff');

%% 2d plots
%b)
x_acc = 0:(dx/4):L;
[lw_acc, t_acc] = lw_2(cfl, a, dt/4, Tmax, 4*N, 4*M, @(t) T_bc(t, Tcool, Thot), alpha, Tcool, Thot,  zeros(1, 4*N+1) + Tcool, dx/4);
[uuw_acc, ~] = upwind_2(cfl, a, dt/4, Tmax, 4*N, 4*M, @(t) T_bc(t, Tcool, Thot), alpha, Tcool, Thot,  zeros(1, 4*N+1) + Tcool);

figure();
subplot(1,2,1);
plot(x, uuw(round(length(t)/2),:), 'linewidth', 2);
hold on;
plot(x, lw(round(length(t)/2),:), 'linewidth', 2);
plot(x_acc, lw_acc(round(size(lw_acc, 1)/2), :), 'linewidth', 2);
plot(x_acc, uuw_acc(round(size(uuw_acc, 1)/2), :), 'linewidth', 2);
xlabel('x');
ylabel('T');
title('2D - at t = 3');
legend('Upwind', 'Lax-Wendroff', 'LW ref', 'Upwind ref');

subplot(1,2,2)
plot(x, uuw(length(t),:), 'linewidth', 2);
hold on;
plot(x, lw(length(t),:), 'linewidth', 2);
plot(x_acc, lw_acc(end, :), 'linewidth', 2);
plot(x_acc, uuw_acc(end, :), 'linewidth', 2);
xlabel('x');
ylabel('T');
title('2D - at t = 6');
legend('Upwind', 'Lax-Wendroff', 'LW ref', 'Upwind ref');





