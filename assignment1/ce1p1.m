%% computer exercise 1
%% part 1

%%  init
set(0,'defaultTextInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');

% fixed values to the problem
a = [1/2, 1/2, sqrt(2)/2]';
alpha = 0.1;
m0 = [1 0 0]';
T = 40;

% vector product matrix
B = [0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0];
% system matrix
A = B*(eye(3) + alpha*B);
% right side of the ODE y'(t) = f(t, y)
F = @(t, y) A*y;

%% runge kutta once
% choose medium number of points for the plots
N = 80;
h = T/N;

[t, m] = rungeKutta(F, h, N, m0);

%% time evolution plots
fig = figure();
plot(t, m, '-o', 'linewidth', 1.2); 
hold on;
grid on;
xlabel('$t$');
ylabel('$m$');
title('Components of the magnetization vector $m$');
legend('$m_1$', '$m_2$', '$m_3$');
exportgraphics(fig, 'ce1p1_components_m.pdf');

%% 3d vector plot

% normalize the vector m
u = m ./ vecnorm(m, 2, 2);
% vector from origin to a
oa = [0 0 0; a'];


fig = figure();
% plot reference vector a
plot3(oa(:, 1), oa(:, 2), oa(:, 3), 'linewidth', 1.2);
hold on;
% plot unit vectors of m
plot3(u(:, 1), u(:, 2), u(:, 3), '-o', 'linewidth', 1.2);
axis equal;
grid on;
set(gca, 'Box', 'on');
view(-54.4, 8.6);
% 
% xlim([-1, 1]);
% ylim([-1, 1]);
% zlim([-1, 1]);

xlabel('x');
ylabel('y');
zlabel('z');

title('Trajectory of the direction of the magnetization vector $m$');
legend('$a$', '$m/\|m\|$');
exportgraphics(fig, 'ce1p1_trajectory_m.pdf');

%% order of accuracy
% lists for h and error
hs = [];
uT = [];
% run runge kutta for differnt step sizes
for N=[20 40 80 160 320 640]
    h = T/N;
    [t, m] = rungeKutta(F, h, N, m0);
    uT = [uT; m(end, :)];
    hs = [hs; h];
end

hs = hs(1:end -1);
err = vecnorm(diff(uT, 1), 2, 2);

 
fig = figure();
loglog(hs, err, '-o');
hold on;
loglog(hs, hs.^3);
xlabel('$h$');
ylabel('error');
grid on;
legend('$\| \tilde{m}_N(T) - \tilde{m}_{2N}(T) \|$', '$h^3$', 'location', 'best');
title('Order of accuracy estimation');
% looks like h^3

exportgraphics(fig, 'ce1p1_order_accuracy.pdf');

%% stability region

p = @(z) z.^3/6 + z.^2/2 + z + 1;
dp = @(z) z.^2/3 + z + 1;

for l = [0, -alpha+1i, -alpha-1i]
    g = @(h) abs(p(h*l))^2 - 1;
    h0 = 2;
    h = fzero(g, h0)
    abs(p(h*l))
end
% h= 2.1420

%% plot unstable/stable
T = 120;
N = floor(T/2.1420);
% unstable
h = T/N;

[t, m] = rungeKutta(F, h, N, m0);

fig = figure();
plot(t, m, '-o', 'linewidth', 1.2); 
hold on;
grid on;
xlabel('$t$');
ylabel('Components of $m$');
title(sprintf('Unstable solution, $h = %.4f$', h));
legend('$m_1$', '$m_2$', '$m_3$');
exportgraphics(fig, 'ce1p1_unstable.pdf');

% stable
h = T/(N+1);

[t, m] = rungeKutta(F, h, N, m0);

fig = figure();
plot(t, m, '-o', 'linewidth', 1.2); 
hold on;
grid on;
xlabel('$t$');
ylabel('Components of $m$');
title(sprintf('Stable solution $h = %.4f$', h));
legend('$m_1$', '$m_2$', '$m_3$');
exportgraphics(fig, 'ce1p1_stable.pdf');
