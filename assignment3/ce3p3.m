%%ce3p2
%%solve the 2D heat equation with semi-discretization and cranck-nicolson
close all
clear all
clc
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

% the heating source term
f = @(x, y) 6000*exp(-5*(x-1)^2-10*(y-1.5)^2);

%%

% spacial step size
h = 0.1;

x = 0:h:Lx;
y = 0:h:Ly;

n = size(x, 2);
m = size(y, 2);

%%
% initial condition
u0 = meshgrid(40 + 72.*x(2:end-1), y);
u0 = reshape(u0', m*(n-2), 1);

[A, b] = mol_2d(Lx, n-1, Ly, m-1, Tleft, Tright, h, f);

dt = 0.2*h*h
tMax = 10;
t = 0:dt:tMax;

%u innehåller tidssteg. Extrahera dem 
u = crankNicolson(dt, size(t, 2), A, b, u0);

%% for u fixed time step
fig = figure();
plot(t, u(:, 619))
hold on
plot(t, 782.43*ones(size(t)),'r-')
xlabel('$\tau \rightarrow$')
ylabel('$T(x=3,y=1,\tau) \rightarrow$')
grid on
exportgraphics(fig, "ce3p3_u_of_tau.pdf");


%% For fixed time steps

[Y,X]=createGrid(Lx, Ly, h);
%time 
for i = [1, 126, 501, 5001]
    figure();
    tx = u(i, :);
    new_u0help = reshape(tx, (n-2), m);
    new_u0left = zeros(1, size(X,2));
    new_u0left(1,:) = 40;
    new_u0right = zeros(1, size(X,2));
    new_u0right(1,:) = 400;
    % add left/right BC
    new_u0= [new_u0left; new_u0help; new_u0right];
    mesh(X, Y, new_u0) 
    xlabel('$y$');
    ylabel('$x$');
    zlabel('$u(x,y,10)$');
end
















