function [t, u] = rungeKutta(f, h, N, y0)
%%perform a runge kutta method to the ODE y'(t) = f(t, y)
%%with step size h, number of steps N and initial condition y0

% time length
T = N*h;
% time points
t = [0:h:T]';
% create matrix for values of u
u = zeros([N+1, size(y0, 1)]);
% set initial condition
u(1, :) = y0';
% run loop to iterate over the steps
for n = 1:N
    % short cut for the current appoximation value
    un = u(n, :)';
    % auxiliary values
    k1 = f(t(n), un);
    k2 = f(t(n) + h, un + h*k1);
    k3 = f(t(n) + h/2, un + h*k1/4 + h*k2/4);
    % next approximation
    u(n + 1, :) = transpose(un + h/6*(k1 + k2 + 4*k3));
end

end