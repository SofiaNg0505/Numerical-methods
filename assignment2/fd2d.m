function [T, X, Y] = fd2d(Lx, N, Ly, M, Tleft, Tright, h, f)

% dimension of the matrix system
dim = (N-1)*(M+1);

% sizes are (N+1, M+1)
[X, Y] = createGrid(Lx, Ly, h);


eY = ones(M+1, 1);
% tridiag matrix (M+1)x(M+1) (-1, 2, 1)
SY = spdiags([-eY, 2*eY, -eY], -1:1, M+1, M+1);
% (M+1)x(M+1) identity matrix
IY = spdiags(eY, 0, M+1, M+1);

eX = ones(N-1, 1);
% tridiag (N-1)x(N-1) matrix (-1, 2, 1)
SX = spdiags([-eX, 2*eX, -eX], -1:1, N-1, N-1);
% NxN identity matrix
IX = spdiags(eX, 0, N-1, N-1);

% change two entries to match the neumann BC on top and bottom
SY(1, 2) = -2;
SY(end, end-1) = -2;

% kron is the tensor product of matrices
A = kron(IY, SX) + kron(SY, IX);
% left part is creating the center diagonals (blocks of SX)
% right part is creating the far off diagonals

A = 1/h^2 * A;

% boundary conditions
% index (x, y)
F = zeros(N-1, M+1);
for i = 1:(N-1)
    for j = 1:(M+1)
        F(i, j) = f(X(i+1, j), Y(i+1, j));
    end
end
F(1, :) = F(1, :) + Tleft/h/h;
F(N-1, :) = F(N-1, :) + Tright/h/h;
F = reshape(F, dim, 1);

% solve the system
T = A\F;

T = reshape(T, N-1, M+1);

T = [Tleft*ones(1, M+1); T; Tright*ones(1, M+1)];

end