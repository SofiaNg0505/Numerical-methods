%% intro
set(0,'defaultTextInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');

r1 = 0.04;
r2 = 10^4;
r3 = 3*10^7;

f = @(t, y) [-r1*y(1)+r2*y(2)*y(3); r1*y(1)-r2*y(2)*y(3)-r3*y(2)*y(2); r3*y(2)*y(2)];

y0 = [1; 0; 0];

%% solve with Runge-Kutta
T = 1;

N = 125*2^3
h = T/N

[t, u] = rungeKutta(f, h, N, y0);

fig = figure();

% tiledlayout(1, 3);

subplot(1, 3, 1);
loglog(t, u(:, 1), 'linewidth', 1.2);
xlabel('$t$');
ylabel('$x_1$');
% title('$x_1$');

subplot(1, 3, 2);
loglog(t, u(:, 2), 'linewidth', 1.2);
xlabel('$t$');
ylabel('$x_2$');

subplot(1, 3, 3);
loglog(t, u(:, 3), 'linewidth', 1.2);
xlabel('$t$');
ylabel('$x_3$');

sgtitle('Runge-Kutta solution, $N = 1000$');
exportgraphics(fig, 'ce1p3_runge_kutta.pdf');

%% ode23
T = 1;

for k = [3, 4, 5, 11]
    relTol = 10^(-k);
    opt = odeset('RelTol', relTol, 'AbsTol', relTol/1000);

    [t, y] = ode23(f, [0, T], y0, opt);
    fprintf("RelTol = 10^(-%d) has %d steps\n", k, size(t, 1));
end

% 10^(-5)
fig = figure();
subplot(1, 2, 1);
relTol = 10^(-5);
opt = odeset('RelTol', relTol, 'AbsTol', relTol/1000);
[t, y] = ode23(f, [0, T], y0, opt);
loglog(t(1:end-2), diff(t(1:end-1)));
xlabel('$t$');
ylabel('$h$');
title('RelTol = $10^{-5}$');
grid on;

% 10^(-11)
% fig = figure();
subplot(1, 2, 2);
relTol = 10^(-11);
opt = odeset('RelTol', relTol, 'AbsTol', relTol/1000);
[t, y] = ode23(f, [0, T], y0, opt);
loglog(t(1:end-2), diff(t(1:end-1)));
xlabel('$t$');
ylabel('$h$');
title('RelTol = $10^{-11}$');
grid on;

sgtitle('Solver \texttt{ode23}');

exportgraphics(fig, 'ce1p3_ode23.pdf')


%% ode23s

T = 1000;

for k = [3, 4, 5, 11]
    relTol = 10^(-k);
    opt = odeset('RelTol', relTol, 'AbsTol', relTol/1000);

    [t, y] = ode23s(f, [0, T], y0, opt);
    fprintf("RelTol = 10^(-%d) has %d steps\n", k, size(t, 1));
end

% 10^(-5)
fig = figure();
subplot(1, 2, 1);
relTol = 10^(-5);
opt = odeset('RelTol', relTol, 'AbsTol', relTol/1000);
[t, y] = ode23s(f, [0, T], y0, opt);
loglog(t(1:end-2), diff(t(1:end-1)));
xlabel('$t$');
ylabel('$h$');
title('RelTol = $10^{-5}$');
grid on;

% 10^(-11)
% fig = figure();
subplot(1, 2, 2);
relTol = 10^(-11);
opt = odeset('RelTol', relTol, 'AbsTol', relTol/1000);
[t, y] = ode23s(f, [0, T], y0, opt);
loglog(t(1:end-2), diff(t(1:end-1)));
xlabel('$t$');
ylabel('$h$');
title('RelTol = $10^{-11}$');
grid on;

sgtitle('Solver \texttt{ode23s}');

exportgraphics(fig, 'ce1p3_ode23s.pdf')
