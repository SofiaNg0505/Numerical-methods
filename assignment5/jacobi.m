function [x, flag, relres, iter, resvec] = jacobi(M, T, b, x0, eps, max_iter)
x = x0;
r = b - (M - T)*x; %residual
iter = 0;
flag = 0;
resvec = [];

while(norm(r)/norm(b) >= eps && iter < max_iter)
    r = b - (M - T)*x;
    x = M\(T*x + b);
    resvec(end+1) = norm(r);
    iter = iter + 1;
end

if(iter >= max_iter)
    flag = 1;
end

relres = norm(r)/norm(b); %relative residual
end