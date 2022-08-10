%% CE5 part 2(a)
fprintf('*** (a) ***\n');

load('cooling_flange.mat', 'A')
% spy(A)

n = length(A);
max_iter = 800;

tol = 1e-4;
rng(42);
% random b vector
b = rand([n,1]);

% pcg
eval1 = tic;
[x, flag, relres, iter, resvec] = pcg(A, b, tol, max_iter);
evaluate_pcg = toc(eval1);
if (flag ~= 0)
    fprintf("The pcg didn't converge, increase the number of iterations!\n");
end

% backslash
eval2 = tic;
x = A\b;
evaluate_backslash = toc(eval2);


figure();
semilogy(0:iter, resvec./norm(b));
hold on;
ylabel("relative residual error");
xlabel("number of iterations");
title(sprintf("Convergence rate (|b|=%f)", norm(b)));
grid on;

fprintf("The pcg after %d iterations with final relative residual %f\nin %f sec while backslash needed %f sec.\n", iter, relres, evaluate_pcg, evaluate_backslash);

%% (b) preconditioning
fprintf('*** (b) ***\n');

% diagonal of A
M = spdiags([diag(A)], 0, n, n);
% incomplete Cholesky factorization
L = ichol(A);

eval1 = tic;
[x, flag, relres, iter, resvec] = pcg(A, b, 0.0001, max_iter, M);
time_for_pcg_M = toc(eval1);
semilogy(0:iter, resvec./norm(b));

eval2 = tic;
[x, flag, relres, iter, resvec] = pcg(A, b, 0.0001, max_iter, L, L');
time_for_pcg_ichol = toc(eval2);
semilogy(0:iter, resvec./norm(b));
legend('pcg', 'pcg (diag)', 'pcg (ichol)');

eval3 = tic;
x = A\b;
time_for_backslash = toc(eval3);

eval4 = tic;
x = A.*inv(M)\(b.*inv(M));
time_for_backslash_M = toc(eval4);

eval5= tic;
x = A.*inv(L.*L')\(b.*inv(L.*L'));
time_for_backslash_ichol = toc(eval5);

fprintf("pcg (prec diag): %f sec\n", time_for_pcg_M);
fprintf("pcg (prec ichol): %f sec\n", time_for_pcg_ichol);
fprintf("backslash: %f sec\n", time_for_backslash);


%% (c) Cooling Flange 
fprintf('*** (c) ***\n');

load('convdiff.mat', 'A')

n = length(A);
max_iter = 800;
tol = 1e-4;

rng(43);
b = rand([n,1]);

% Plot convergence with pcg
[x, flag, relres, iter, resvec] = pcg(A, b, tol, 100);

figure()
plot(1:length(resvec), resvec./norm(b));
hold on;
ylabel("relative residual error");
xlabel("number of iterations");
title(sprintf("Convergence rate (|b|=%f)", norm(b)));
grid on;
legend('pcg');

% Compare time for backslash and GMRES
start = tic;
[x, flag, relres, iter, resvec] = gmres(A, b, [], tol, max_iter);
time_gmres = toc(start);

figure();
semilogy(0:iter(2), resvec./norm(b));
hold on;
xlabel("number of iterations");
ylabel("relative residual error");
title(sprintf("Convergence rate (|b|=%f)", norm(b)));
grid on;

start = tic;
x = A\b;
time_backslash = toc(start);

fprintf("The gmres after %d iterations with final relative residual %f\nin %f sec while backslash needed %f sec.\n", iter(2), relres, time_gmres, time_backslash);

% preconditioners
M = spdiags([diag(A)], 0, n, n);
[L, U]= ilu(A);

% Timing
start = tic;
[x, flag, relres, iter, resvec] = gmres(A, b, [], tol , max_iter, M);
time_gmres_M  = toc(start);
semilogy(0:iter(2), resvec./norm(b));

% timing for gmres with LU as precondition
start = tic;
[x, flag, relres, iter, resvec] = gmres(A, b, [], tol, max_iter, L, U);
time_gmres_LU  = toc(start);
semilogy(0:iter(2), resvec./norm(b));
legend('gmres', 'gmres (diag)', 'gmres (ilu)');

 
fprintf("gmres: %f sec\n", time_gmres);
fprintf("gmres (prec diag): %f sec\n", time_gmres_M);
fprintf("gmres (prec ilu): %f sec\n", time_gmres_LU);
fprintf("backslash: %f sec\n", time_backslash);


%%
% eval3= tic;
% x = A.*inv(M)\(b.*inv(M));
% gm_res_backslash_M  = toc(eval3)
% 
% 
% eval4= tic;
% x = A.*inv(L.*U)\(b.*inv(L.*U));
% gm_res_backslash_LU  = toc(eval4)
% 
% eval5= tic;
% x = A\b;
% gmres_Backslash  = toc(eval5)

