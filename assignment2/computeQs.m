function Q = computeQs(h, a, b, L, Q0, N)
    Q = zeros(N+1, 1);
    z = 0:h:L;
    
    idx = ((z > a) & (z < b));
    
    Q(idx) = Q0*sin((z(idx) - a)/(b-a)*pi);
    
    Q = Q(1:N);
end