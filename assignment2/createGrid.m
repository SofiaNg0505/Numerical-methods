function [X, Y] = createGrid(Lx, Ly, h)
%%[X, Y] = createGrid(Lx, Ly, h)
%%creates a mesh grid of the area (Lx, Ly) with step size h
%%indexing of X,Y is of kind (x, y)

x = 0:h:Lx;
y = 0:h:Ly;

[Y, X] = meshgrid(y, x);
end