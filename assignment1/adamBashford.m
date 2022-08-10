function [t, u] = adamBashford(f, h, N, y0)
%%perform a adams bashford method to the ODE y'(t) = f(t, y)
%%with step size h, number of steps N and initial condition y0
if(size(y0, 1)) ~= 4
    error("The number of initial values for the adams bashford method is not correct!");
end

% time length
T = N*h;
% time points
t = [0:h:T]';
% create matrix for values of u
u = zeros([N+1, size(y0, 2)]);

% set initial condition
u(1:4, :) = y0;

for n = 4:N
    u(n+1, :) = transpose(u(n,:)' + (55*f(t(n), u(n,:)) - 59*f(t(n-1), u(n-1,:)) + 37*f(t(n-2), u(n-2,:)) - 9*f(t(n-3), u(n-3,:)))/24*h);
end

end