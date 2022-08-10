% part 2 fdm in 2D

%%  init
set(0,'defaultTextInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');

%% constants
Lx = 5;
Ly = 2;

% dirichlet boundary conditions on the left/ride site (x=const)
Tleft = 40;
Tright = 400;
% neumann condiditons on top and bottow (y=const)

%% part (a)
% step size
h = 0.1;
N = 50; % x steps
M = 20; % y steps

f = @(x, y) 100;

[T, X, Y] = fd2d(Lx, N, Ly, M, Tleft, Tright, h, f);
[t, x, y] = getPoint(T, X, Y, h);
fprintf("h = %f, (%f, %f): %f\n", h, x, y, t)

fig = figure();
mesh(X, Y, T, 'facecolor', 'flat');
xlabel('$x$');
ylabel('$y$');
zlabel('$T$');
title('Temperature distribution in a rectangle with $f = 100$');
% exportgraphics(fig, 'ce2p2_constant.pdf');

%% part (c)

hs = [];
ts = [];

g = @(x, y) 6000*exp(-5*(x-1)^2-10*(y-1.5)^2);

for k = 0:6
% step size
h = 0.1/2^k;
N = 50*2^k; % x steps
M = 20*2^k; % y steps

[T, X, Y] = fd2d(Lx, N, Ly, M, Tleft, Tright, h, g);

[t, x, y] = getPoint(T, X, Y, h);
fprintf("k = %d, h = %f, (%f, %f): %f\n", k, h, x, y, t);

hs = [hs; h];
ts = [ts; t];

end

%%
fig = figure();
loglog(hs(2:end), diff(ts), '-x');
hold on;
loglog(hs(2:end), hs(2:end).^2);
grid on;
legend('$T_{N}(3, 1) - T_{2N}(3, 1)$', '$h^2$', 'location', 'northwest');
xlabel('$h$');
ylabel('estimated error');
% exportgraphics(fig, 'ce2p2_order_accuracy.pdf');


%% visualize part (c)
h = 0.05;
N = 100; % x steps
M = 40; % y steps

[T, X, Y] = fd2d(Lx, N, Ly, M, Tleft, Tright, h, g);

fig1 = figure();
mesh(X, Y, T,'facecolor', 'flat');
xlabel('$x$');
ylabel('$y$');
zlabel('$T$');
% exportgraphics(fig1, 'ce2p2_c_mesh.pdf');

% hold on;
fig2 = figure();
imagesc([0, 5], [0, 2], T');
xlabel('$x$');
ylabel('$y$');
zlabel('$T$');
axis xy;
colorbar
% exportgraphics(fig2, 'ce2p2_c_imagesc.pdf');

fig3 = figure();
contour(X, Y, T, 'showtext', 'on', 'linewidth', 2);
xlabel('$x$');
ylabel('$y$');
zlabel('$T$');
grid on;
% exportgraphics(fig3, 'ce2p2_c_contour.pdf');
