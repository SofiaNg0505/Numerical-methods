function [T, t] = upwind_2(cfl, a, dt, Tmax, N, M, T_bc, ALPHA, Tcool, Thot, T0)
%%[u, t] = upwind(cfl, a, dt, Tmax,N,M, u0, g)
%%Upwind method
%%cfl number, constant a, time step dt, finish time T, N spacial points, 
%%initial condition u0, left boundary g

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


if (a > 0)
    
    % loop thru time
    for n = 1:(M-1)
        % BC on the left
        T(n+1, 1) = T_bc(t(n+1));
        % loop thru space
        for j = 2:(N+1)
            T(n+1, j) = T(n, j) - a * cfl * ( T(n, j) - T(n, j-1) )-dt*(ALPHA*T(n, j)-ALPHA*Tcool);
        end
    end
    
elseif (a < 0)
    
    % loop thru time
    for n = 1:(M-1)
        % BC on the left
        T(n+1, 1) = T_bc(t(n+1));
        % loop thru space
        for j = 2:N
            T(n+1, j) = T(n, j) - a * cfl (T(n, j+1) - T(n, j) -dt*( ALPHA * T(n, j)- ALPHA * Tcool)) ;
        end
        
        % numerical BC on the right
        T(n+1, N+1) = 2*T(n+1, N) - T(n+1, N-1);
    end
else
    error('The variable a should not be zero!');
end

end