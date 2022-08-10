clear all,close all, clc
format short
set(groot,'DefaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')
warning('off','all')
warning

global mi R rho title_ptr
mi = 4*10^-5;
R = 0.1;
rho = 1000;
u0 = 0.15;
P0 = 10^4;

Z_SD=[];
N_1 = [100:100:1000];
for nn = N_1
    N = nn-1;
    u_0 = ones(N+1,1)*u0;
    sigma_0 = 0;
    v_0 = sparse(N,1);
    h_zMax = 0.05;
    z_max = 70 + 0.04;
    TOL = 2*10^-14;
    TOL_u = 4/3*(1/(N+1))^2;
    r = (0:R/(N+1):R-R/(N+1))';
    u_r = @(r) 2*u0*(1-(r/R).^2);
    
    disp("For N+1 = " + string(nn)+":")
    tic
    [X, P, z_SD]=  newton(u_0, sigma_0, v_0, P0, N, h_zMax, z_max, TOL, u_r(r), TOL_u);
    toc
    
    Z_SD = [Z_SD,z_SD];
end
Z_SD
title_ptr = [];
figure;
plot(N_1, Z_SD)
xlabel("N+1")
ylabel("z_{SD}")
title("z-value for Stationary Distribution, $z_{SD}$");
title_ptr = [title_ptr,get(gca, 'title')];

size(X)
z = [0:0.001:0.005, (0.005+0.005):0.005:0.04, (0.04+h_zMax):h_zMax:z_SD];
size(z)

r = (0:R/(N+1):R)';
u = [X(1:N+1,:);sparse(1,size(X,2))];
sigma = X(N+2,:);
v = [sparse(1,size(X,2));X(N+3:end,:);sparse(1,size(X,2))];

TOL_str = "$\frac{|X_{k+1} - X_{k}|}{|X_{k}|} <$ " + string(TOL);
TOL_u_str = "$\frac{|u^{n} - u|}{|u|} <\frac{4}{3}$" + string((1/(N+1))^2);
hr_hz_str =  "$h_r =$ " + string(R/(N+1)) + " and $h_z=$ " + string(h_zMax)+"(for $z \geq 0.04$)";
title_str = " with " + TOL_str + ", " + TOL_u_str + ", " + newline + hr_hz_str;
figure;
plot_rz(u,r,z, "u", title_str)
figure;
plot_rz(v,r,z, "v",title_str)


figure;
plot(r,u_r(r), '*')
hold on
plot(r, u(:,end), 'LineWidth', 1)
xlabel("r")
ylabel("u")
u_str = "$u = 2u_{0}(1-(\frac{r}{R})^{2})$";
leg_ptr = legend(u_str, "Approximated u(r, " +string(z_SD)+")");
set(leg_ptr, 'Interpreter','latex', 'FontSize', 12)
title("Stationary Distribution of u" +title_str)
title_ptr = [title_ptr,get(gca, 'title')];


figure;
plot(z,P)
xlabel("z")
ylabel("P")
title("Approximation of P" + title_str)
title_ptr = [title_ptr,get(gca, 'title')];

figure;
plot(z,sigma)
xlabel("z")
ylabel("\sigma")
title("Approximation of $\sigma$" + title_str)
title_ptr = [title_ptr,get(gca, 'title')];

set(title_ptr, 'Interpreter','latex','FontSize', 12);
function plot_rz(X,r,z, x_name, title_str)
    global title_ptr

    [Z,RR] = meshgrid(z,r);
    mesh(Z, RR,X)
    xlabel("z")
    ylabel("r")
    zlabel(x_name)

    title("Approximation of "+ x_name + title_str)
    title_ptr = [title_ptr,get(gca, 'title')];
end


function [X, P, z_SD] = newton(u_0, sigma_0, v_0, P0, N, h_zMax, z_max, TOL,u_r, TOL_u)
    global R rho
    
    h_r = R/(N+1);
    M = round((z_max-0.04)/h_zMax);
    X = sparse(2*(N+1),M+7+5+1);
    P = sparse(1,M+7+5+1);
    P(1) = P0;
    X(:,1) = [u_0; sigma_0; v_0];
    x_k = [u_0; sigma_0; v_0];
    x_k_1 = x_k;
    z_SD = 0;
    
    for i = 1:5+7+M
        if i < 6
            h_z = 0.001;
        elseif i < 7+5+1
            h_z = 0.005;
        else
            h_z = h_zMax;
        end
        J = jacobian(x_k(1:N+1), x_k(N+3:end), X(1:N+1,i), h_r,N, h_z);
        g = g_x(x_k(1:N+1), x_k(N+2), x_k(N+3:end), X(1:N+1,i), h_r, N, h_z);
       
        x_k_1 =  x_k -J\g;
        while norm(x_k_1-x_k)/norm(x_k)> TOL
            x_k = x_k_1;
            J = jacobian(x_k(1:N+1), x_k(N+3:end), X(1:N+1,i), h_r,N, h_z);
            g = g_x(x_k(1:N+1), x_k(N+2), x_k(N+3:end), X(1:N+1,i), h_r, N, h_z);
            
            x_k_1 =  x_k - J\g; 
        end
        
        P(i+1) = P(i) + h_z*rho*x_k(N+2);
        X(:,i+1) = x_k;
        if norm(x_k(1:N+1)-u_r)/norm(u_r) < TOL_u
            z_SD = (i-(7+5))*h_zMax + 0.04;
            P = P(1:1:i+1);
            X = X(:,1:i+1);            
            break;
        end
    end   
end

function g = g_x(u_n1, sigma_n1, v_n1, u_n, h_r, N, h_z)
    global mi 
    a_j = mi/h_r^2*(1 - (1./(2*(1:N)'))) + (1/2/h_r)*v_n1;

    b_j = -2*mi/h_r^2*ones(N,1);

    c_j = mi/h_r^2*(1 + 1./(2*(1:N-1)')) - (1/2/h_r)*v_n1(1:N-1);
    
    g1= -3*u_n1(1)+4*u_n1(2)-u_n1(3);

    g2 = u_n1(2:N+1).*(u_n1(2:N+1) - u_n(2:N+1))/h_z + ...
     + (-1)*(a_j.*u_n1(1:N) + b_j.*u_n1(2:N+1) + [c_j;0].*[u_n1(3:end);0])+sigma_n1;

    g3 = (u_n1(1)-u_n(1))/h_z + (4*v_n1(1)-v_n1(2))/h_r;

    g4 = (u_n1(2:N+1) - u_n(2:N+1))/h_z + ...
        + ([v_n1(2:N);0] - [0;v_n1(1:N-1)])/(2*h_r) + v_n1.*(1./(1:N)')/h_r;
    g = [g1;g2;g3;g4];
end

function J = jacobian(u_n1, v_n1, u_n, h_r, N, h_z)
    global mi
    a_j = mi/h_r^2*(1 - (1./(2*(1:N)'))) + (1/2/h_r)*v_n1;

    b_j = -2*mi/h_r^2*ones(N,1);

    c_j = mi/h_r^2*(1 + 1./(2*(1:N-1)')) - (1/2/h_r)*v_n1(1:N-1);
    
    I =  sparse(eye(N+1,N+1));
    J1 = [sparse(1,N+1);I(1:N,1:N).*(-a_j) sparse(N,1)] +...
    + [sparse(1,N+1);sparse(N,1) I(1:N,1:N).*(2*u_n1(2:N+1)/h_z - u_n(2:N+1)/h_z - b_j)] +...
    + [-3, 4, -1, sparse(1,N-2);sparse(N-1,2) I(1:N-1,1:N-1).*(-c_j); sparse(1,N+1)];

    J2 = [sparse(1,N+1);ones(N-1,1), I(1:N-1,1:N-1).*((-1/2/h_r)*u_n1(1:N-1) + (1/2/h_r)*u_n1(3:N+1)), sparse(N-1,1); 1, sparse(1,N-1), (-1/2/h_r)*u_n1(N)];

    J3 = I(1:N+1,1:N+1).*(1/h_z);

    J4 = [0, 4/h_r, -1/h_r, sparse(1,N-2); sparse(N,1), I(1:N,1:N).*((1/h_r)./[1:N])] +...
        + [sparse(2,N+1); sparse(N-1,1), I(1:N-1,1:N-1).*(-1/2/h_r),sparse(N-1,1)] +...
        + [sparse(1,N+1); sparse(N-1,2), I(1:N-1,1:N-1).*(1/2/h_r);sparse(1,N+1)];

    J = [J1, J2; J3, J4];
end
