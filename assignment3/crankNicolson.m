function u = crankNicolson(dt, k, A, b, u0)
    u = zeros(k, size(u0, 1));
    u(1, :) = u0';
    
    dim = size(A, 1);
    
    spI = spdiags(ones(dim), 0, dim, dim);
    Ml = spI - dt/2 * A;
    Mr = spI + dt/2 * A;
    
    for i = 2:k
        u(i, :) =  Ml\(Mr*(u(i-1, :))' + dt*b);
    end
end