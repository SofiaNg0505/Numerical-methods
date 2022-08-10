function [u, t] = upwind(cfl, a, dt, T, N, u0, g)
%%[u, t] = upwind(cfl, a, dt, N, M, u0, g)
%%Upwind method
%%cfl number, constant a, time step dt, finish time T, N spacial points, 
%%initial condition u0, left boundary g


% init time points
t = 0:dt:T;
M = length(t);

% init array with dimension (time x space)
u = zeros(M, N+1);

if (size(u0) == [1, N])
    u(1, :) = u0;
else
    u(1, :) = u0';
end
u(1, 1) = g(0);


if (a > 0)
    
    % loop thru time
    for n = 1:(M-1)
        % BC on the left
        u(n+1, 1) = g(t(n+1));
        % loop thru space
        for j = 2:(N+1)
            u(n+1, j) = u(n, j) - a * cfl * ( u(n, j) - u(n, j-1) );
        end
    end
    
elseif (a < 0)
    
    % loop thru time
    for n = 1:(M-1)
        % BC on the left
        u(n+1, 1) = g(t(n+1));
        % loop thru space
        for j = 2:N
            u(n+1, j) = u(n, j) - a * cfl ( u(n, j+1) - u(n, j) );
        end
        
        % numerical BC on the right
        u(n+1, N+1) = 2*u(n+1, N) - u(n+1, N-1);
    end
    
else
    error('The variable a should not be zero!');
end

end