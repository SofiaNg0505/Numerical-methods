function T = fd1d(k1, k2, k3, d0, d1, d2, Q, T0, N)
%%fd1d(k1, k2, k3, d0, d1, d2, Q, T0, N)
%%Finite Differences in 1 dimension
%
% k1, k2, k3 are constants for the diagonals
% d0, d1, d2 co
% me from the Robin boundary conditions
% Q is the right hand side of the linear system
% T0 is the initial value of the temperature
% N is the number of sampling points, 
% so make sure the length of Q is N!

e = ones(N-1, 1);

A = spdiags([k1*e, k2*e, k3*e], -1:1, N-1, N-1);
A(N-1, N-2) = A(N-1, N-2) + k3*d2;
A(N-1, N-1) = A(N-1, N-1) + k3*d1;

if size(Q, 1) == 1
    b = Q(1:N-1)';
else 
    b = Q(1:N-1);
end

b(1) = b(1) - k1 * T0;
b(N-1) = b(N-1) - k3 * d0;

% Ax = b
x = A\b;

T = [T0; x; d0 + d1*x(N-1) + d2*x(N-2)];
end