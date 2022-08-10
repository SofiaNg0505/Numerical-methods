function [A, b, x] = mol(a, dx, n)

e = ones(n+1, 1);

x = (0:dx:1)';

k = a / dx / dx;

A = spdiags([k*e, -k*2*e, k*e], -1:1, n+1, n+1);
% addapt to neumann condition as the right boundary
A(n+1, n) =  2* A(n+1, n);

b = zeros(n+1, 1);
b(1) = k;

end