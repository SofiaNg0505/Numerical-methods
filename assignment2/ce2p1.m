%% part 1 : finite difference approximation in 1D

%%  init
set(0,'defaultTextInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');

%% constants
L = 1;
a = 0.2;
b = 0.3;
Q0 = 4000;
alpha0 = 100;
Tout = 20;
T0 = 50;

%% (a)
v = 10;

hs = [0.05, 0.025, 0.0125, 0.00625];
TN = [];

figa = figure();

for i = 1:4

    h = hs(i);
    N = 1/h;
    z = (0:h:L)';
    
    Q = computeQs(h, a, b, L, Q0, N);
    k = computeKs(h, v);
    d = computeDs(h, v, alpha0, Tout);

    T = fd1d(k(1), k(2), k(3), d(1), d(2), d(3), Q, T0, N);
   
    % save the values of T(L) for the error estimate later
    TN = [TN; T(N+1)];

    plot(z, T);
    hold on;    
end

xlabel('$z$');
ylabel('$T$');
legend('$h=0.05$','$h=0.025$','$h=0.0125$','$h=0.00625$', 'location', 'southwest');
grid on;
% exportgraphics(figa, 'ce2p1_a.pdf');


% decreasing h for every iteration, giving better approximations of T for
% every iterations. Taking the temperature at z = L for every iteration and
% making a log plot of those differences. Comparing with the different very
% small steplengths

figerror = figure();
loglog(hs(2:end),diff(TN));
hold on;
% plot h and h^2 for reference 
loglog(hs(2:end), hs(2:end));
loglog(hs(2:end), hs(2:end).^2);
legend('data', '$h$', '$h^2$', 'location', 'northwest');
grid on;
xlabel('$h$');
ylabel('error estimate');
% exportgraphics(figerror, 'ce2p1_order_accuaracy.pdf')


%% (b)

N = 320;
h = L/N; % h = 0.00625
fprintf("h = %f\n", h);
% does not depend on v
Q = computeQs(h, a, b, L, Q0, N);
z = 0:h:L;

lgd = [];
fig = figure();

for v = [1, 10, 30, 100]
    % depends on v
    k = computeKs(h, v);
    d = computeDs(h, v, alpha0, Tout);

    % compute solution
    T = fd1d(k(1), k(2), k(3), d(1), d(2), d(3), Q, T0, N);
    fprintf("v = %d, T(L) = %f\n", v, T(N+1));
    plot(z, T);
    hold on;
    
    lgd = [lgd, sprintf("$v=%d$", v)];
end


xlabel("horizontal position $z$");
ylabel("temperature $T$");
title("Temperature along a pipe");
legend(lgd)

grid on;
% exportgraphics(fig, "ce2p1_vs.pdf");