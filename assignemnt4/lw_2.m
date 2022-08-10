function [T, t] = lw_2(cfl, a, dt, Tmax, N, M, T_bc, alpha, Tcool, Thot, T0, dx)
%%[u, t] = lw(cfl, a, dt, N, M, u0, g)
%%Lax-Wendrof
%%cfl number, constant a, time step dt, finish time T, N spacial points, 
%%initial condition u0, left boundary g
sigma = a*dt/dx;
% init time points
t = 0:dt:Tmax;
M = length(t);

% init array with dimension (time x space)
T = zeros(M, N+1);

if (size(T0) == [1, N])
    T(1, :) = T0;
else
    T(1, :) = T0';
end
T(1, 1) = T_bc(t(1));


% loop thru time
for n = 1:(M-1)
    % BC on the left
    T(n+1, 1) = T_bc(t(n+1));
    % loop thru space
    for j = 2:N
        T(n+1, j) =  T(n, j) - ( sigma*(1-alpha*dt)/2 )*( T(n ,j+1) - T(n,j-1) ) + sigma*sigma * ( T(n ,j+1) - 2*T(n , j) + T(n, j-1) )/2 - dt*(1-alpha*dt/2)*( alpha*T(n ,j) - alpha*Tcool );
    end
    
    % numerical BC on the right
    T(n+1, N+1) = 2*T(n+1, N) - T(n+1, N-1);
end

end