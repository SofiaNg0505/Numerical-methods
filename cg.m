function [x, flag, relres, iter, resvec] = cg(A, b, x0, eps, max_iter)

% init values
x = x0;
p = x0;
ap = A*p;
r = b;
beta = 0;
iter = 0;
flag = 0;
resvec = [];

while(norm(r)/norm(b) >= eps && iter < max_iter)
  resvec(end+1) = norm(r);
  p = r + beta * p;
  ap = A * p;
  alpha = (p' * r) / (p' * ap);
  x = x + alpha * p;
  r = r - alpha * ap;
  beta = (r' * r) / resvec(end)^2; %should be resvec()'* resvec()
  iter = iter + 1;
end

if( iter >= max_iter)
  flag = 1;
end

relres = norm(r)/norm(b);

end