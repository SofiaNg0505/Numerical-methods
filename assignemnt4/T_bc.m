function T_0t = T_bc(t, Tcool, Thot)
T_0t = zeros(1, length(t));
for timesteps = 1:length(t)
    if 0 <= t(timesteps) && t(timesteps) <= 0.1250 
        T_0t(timesteps) = Tcool + (Thot-Tcool) * sin(4*pi*t(timesteps));
    elseif 0.1250 <= t(timesteps) && t(timesteps) <= 1
        T_0t(timesteps) = Thot;
    elseif  t(timesteps) > 1 
        T_0t(timesteps) = Thot + Tcool * sin(5*pi * (t(timesteps) - 1));
    end
end
    
end