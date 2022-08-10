clear all

h=0.001;
T=20;
% choose N by T and h
N=ceil(T/h);
mu=1/82.45;

position_earth=[-mu;0];
position_moon=[1-mu;0];

% initial condition
r0 = [1.2, 0, 0, -1]';

% function for ODE y'(t) = F(t, y)
F = @(t, r) [r(3), ...
            r(4) ...
            -(1-mu)*(r(1)+mu)/((r(1)+mu)^2+(r(2))^2)^(3/2)-mu*(r(1)+mu-1)/((r(1)+mu-1)^2+(r(2))^2)^(3/2)+2*r(4)+r(1), ...
            -(1-mu)*r(2)/((r(1)+mu)^2+(r(2))^2)^(3/2)-mu*r(2)/((r(1)+mu-1)^2+(r(2))^2)^(3/2)-2*r(3)+r(2)]';


% Runge-Kutta just for 3 steps
[~, r] = rungeKutta(F, h, 3, r0);

% Adam-Bashford for the rest
[t, r] = adamBashford(F, h, N, r);


% plot(position_earth(1),position_earth(2),'b*')
% hold on
% plot(position_moon(1),position_moon(2),'r*')
% hold on

% % Plotting the trajectory
% plot(r(:,1), r(:,2))
% legend('Position of the earth','Position of the moon','Trajectory of the satellite during accurate time')
% xlabel('x')
% ylabel('y')
% axis equal
% figure()

% % Plotting the velocity and/or the position coordinates
% plot(t,r(:,3:4))
% legend('x''-component', 'y''-component')
% xlabel('t')
% ylabel('r''(t)')

% % % Solving the ODE with the different methods
TOL = 0.25;
[t , r_E]= euler(F, h, N, r0);
[t, r_Ehalf]= euler(F,h/2, 2*N, r0);

[t, r_RK] = rungeKutta(F, h, N, r0);
[t, r_RKhalf] = rungeKutta(F, h/2, 2*N, r0);
% 
[t, r_AB] = adamBashford(F, h, N, r_RK(1:4,:));
[t, r_ABhalf] = adamBashford(F, h/2, 2*N, r_RK(1:4,:));
% 
r_method=[r_AB(:,1) r_AB(:,2)];
r_methodhalf_i=[r_ABhalf(:,1) r_ABhalf(:,2)];
r_methodhalf = r_methodhalf_i(1:2:end,:);

for tlimit=1:N 
    if abs(r_method(tlimit)-r_methodhalf(tlimit))> TOL
    limit=abs(r_method(tlimit-1)-r_methodhalf(tlimit-1)) 
    T_acc=(tlimit-1)*h
        break
    end  
end

t=[0:h:T_acc]
plot(r_method(1:length(t),1),r_method(1:length(t),2),'color',[1, 0.5, 0])
hold on
plot(position_earth(1),position_earth(2),'b*')
hold on 
plot(position_moon(1),position_moon(2),'r*')
xlabel('x')
ylabel('y')
axis equal
% % % % % % % % Analytical solution to the ODE  % % % % % % % %

options= odeset(RelTol= 10^-4)
[t,r] = ode23(F,[0 T_acc],r0,options)
hold on
plot(r(:,1),r(:,2))
hold on
% plot(position_earth(1),position_earth(2),'b*')
% hold on 
% grid on
% plot(position_moon(1),position_moon(2),'r*')
legend('Solution with Adam Bashford','Position of the moon','Position of the earth','Soltuion with ode23')
xlabel('x')
ylabel('y')
figure()

% % % %Plot timesteps % % % % %

hh=diff(t)
Hmax=max(hh)
Hmin=min(hh)
plot(t(1:size(hh)),hh)
xlabel('Time(t)')
ylabel('Size of time steps')

