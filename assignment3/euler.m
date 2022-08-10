function [t, u] = euler(f, h, N, y0)
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
    % next approximation
    u(n + 1, :) = un + h*f(t(n),un);
end

end