set(0,'defaultTextInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');


%% general constants

tau = 2.2;
D = 4.5;
a = 1;
linewidth = 2;
Tmax = 6;

%% plotting

figure();
% plot();
xlabel('x');
ylabel('u');
xlim([0, D]);
ylim([-1.5, 1.5]);
labels = {};
l = legend();
l.ItemHitFcn = @hitcallback_ex1;
grid on;
hold on;

%%

while (1)
    mode = input("Which do you want to specify? Enter a number:\n(1) CFL number\n(2) dx\n(3) Custom\n");
    switch (mode)
        case 1
            cfl = input("Enter your CFL number:");
            N = 99;
            dx = D/N;
            x = 0:dx:D;
            dt = cfl*dx;
            t = 0:dt:Tmax;
            M = length(t);
            
        case 2
            dx = input("Enter the spacial step dx:");
            x = 0:dx:D;
            N = length(x) - 1;
            cfl = 1;
            dt = cfl * dx;
            t = 0:dt:Tmax;
            M = length(t);
        
            case 3
            dx = input("Enter the spacial step dx:");
            cfl = input("Enter your CFL number:");
            x = 0:dx:D;
            N = length(x) - 1;
            dt = cfl * dx;
            t = 0:dt:Tmax;
            M = length(t);
        
        otherwise
            warning("Your case does not exist! Aborting.");
            break;
    end % mode
    
    % initial condition
    u0 = zeros(1, N+1);
    
    fprintf('CFL: %f, dt = %f, T = %f, dx = %f, N = %d, M = %d\n', cfl, dt, t(end), dx, N, M);
    
    algo = input("Which method do you want to use? Enter a number:\n(1) Lax-Friedrich\n(2) Upwind\n(3) Lax-Wendroff\n(4) exact solution\n");
    
    bc = input("Choose the boundary condition:\n(1) sine wave\n(2) square wave\n");
    
    switch (bc)
        case 1
            g = @(t) gsin(t,tau);
        case 2
            g= @(t) gsq(t,tau);
        otherwise
            g = @(t) gsin(t,tau);
            warning("You selected the wrong boundary condition! Continuing with sine wave.")
    end % bc
    
    switch (algo)
        case 1
            % lax-friedrich
            [u, t] = lxf(cfl, a, dt, t(end), N, u0, g);
            plot(x, u(end, :), 'linewidth', linewidth);
            hold on;
            labels{end+1} = sprintf("LxF (cfl=%.2f, dx=%0.4f, t=%.1f)", cfl, dx, t(end));
            legend(labels, 'location', 'northoutside');
            
        case 2
            % upwind
            [u, t] = upwind(cfl, a, dt, t(end), N, u0, g);
            plot(x, u(end, :), 'linewidth', linewidth);
            hold on;
            labels{end+1} = sprintf("Up (cfl=%.2f, dx=%0.4f, t=%.1f)", cfl, dx, t(end));
            legend(labels, 'location', 'northoutside');
            
                
        case 3
            % lax-wendroff
            [u, t] = lw(cfl, a, dt, t(end), N, u0, g);
            plot(x, u(end, :), 'linewidth', linewidth);
            hold on;
            labels{end+1} = sprintf("LW (cfl=%.2f, dx=%0.4f, t=%.1f)", cfl, dx, t(end));
            legend(labels, 'location', 'northoutside');
            
            
        case 4
            plot(x, g(t(end) - x), 'linewidth', linewidth);
            hold on;
            labels{end+1} = sprintf("exact (dx=%0.4f, t=%.1f)", dx, t(end));
            legend(labels, 'location', 'northoutside');
            
            
        otherwise
            warning("You selected a non existing method! Ignoring it.");
            
    end % algo
    
    
end

function hitcallback_ex1(src,evnt)

if strcmp(evnt.Peer.Visible,'on')
    evnt.Peer.Visible = 'off';
else 
    evnt.Peer.Visible = 'on';
end

end
