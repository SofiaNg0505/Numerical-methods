%% CE5 part 1

% k x 3 matrix
% (n, d, method)
% method: (1) Jacobi, (2) CG
ndm_pairs = [
    20000, 1, 1;
    20000, 1, 2;
    30000, 1, 1;
    30000, 1, 2;
    40000, 1, 1;
    40000, 1, 2;
    150, 2, 1;
    150, 2, 2;
    200, 2, 1;
    200, 2, 2;
    27, 3, 1;
    27, 3, 2;
    34, 3, 1;
    34, 3, 2
    ];

%% (a)
fprintf("*** (a) ***");
max_iter = 2000;
rng(43);
eps = 1e-9;
labels = {};

for i = 1:size(ndm_pairs, 1)
    
    n = ndm_pairs(i, 1);
    d = ndm_pairs(i, 2);
    method = ndm_pairs(i, 3);
    N = n^d;
    
    A = lap(n, d);
    
    M = spdiags([diag(A)], 0, N, N);
    T = M - A;
    
    b = rand(N, 1);
    x0 = zeros(N, 1);
    
    if (method == 1)
        % jacobi
        start = tic;
        [x_j, flag, relres, iter, resvec_j] = jacobi(M, T, b, x0, eps, max_iter);
        stop = toc(start);
        fprintf("N: %d, Jacobi after %d iterations in time %f: flag %d and relres: %f\n", N, iter, stop, flag, relres);
        
        semilogy(1:iter, resvec_j./norm(b));
        hold on;
        labels{end+1} = sprintf('Jacobi, n=%d, d=%d', n, d);
        
    elseif (method == 2)
        % CG
        start = tic;
        [x_cg, flag, relres, iter, resvec_cg] = cg(A, b, x0, eps, max_iter);
        stop = toc(start);
        fprintf("N: %d, CG after %d iterations in time %f: flag %d and relres: %f\n", N, iter, stop, flag, relres);
        
        semilogy(1:iter, resvec_cg./norm(b));
        hold on;
        labels{end+1} = sprintf('CG, N=%d, n=%d, d=%d', N, n, d);
    end
end

ylabel("relative residual error");
xlabel("number of iterations");
l = legend(labels);
l.ItemHitFcn = @hitcallback_ex1;
title("Convergence rate");
grid on;

%% (b)
fprintf("*** (b) ***\n");
eps = 1e-10;
max_iter = 40000;
rng(44);

for i = 1:size(ndm_pairs, 1)
    
    n = ndm_pairs(i, 1);
    d = ndm_pairs(i, 2);
    method = ndm_pairs(i, 3);
    if (method ~= 2)
        continue;
    end
    N = n^d;
    
    A = lap(n, d);
    b = rand(N, 1);
    x0 = zeros(N, 1);

    start = tic;
    [x_cg, flag, relres, iter, resvec_cg] = cg(A, b, x0, eps, max_iter);
    time_cg = toc(start);

    start = tic;
    x_backslash  = A\b;
    time_backslash = toc(start);

    reldiff = norm(x_cg - x_backslash) / norm(x_backslash);

    fprintf("N=%d, n=%d, d=%d, time CG(%d): %f, time Backslash: %f, rel. diff. %e\n", N, n, d, iter, time_cg, time_backslash, reldiff);

end

%%
function hitcallback_ex1(src,evnt)

if strcmp(evnt.Peer.Visible,'on')
    evnt.Peer.Visible = 'off';
else 
    evnt.Peer.Visible = 'on';
end

end