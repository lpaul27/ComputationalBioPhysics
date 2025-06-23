function [Cradius,x, y, vx, vy, vel_ang, x_time, y_time, theta_time, RadTracker, R, G, B, Pressure, exempt] = RadGrowth(Cradius0, Pressure, x, ...
    y, vel_ang, vx, vy, x_time, y_time, time, theta_time, RadTracker, R, G, B, exempt)
% Radial growth of a cell with each iteration

global dt NumCells critRad Ccyclet critical_pressure vels_med daughter_noise ...
    death_pressure death_rate

Cradius = Cradius0;
growth_rate = (pi * critRad^2) / (2* Ccyclet);

%% loops to check growth
for i = 1:NumCells
    % Loops over all cells
    % Check conditional of dormancy to grow
    if(Pressure(i,1) < critical_pressure && Cradius(i,1) <= critRad && exempt(i,1))
        Cradius(i,1) = Cradius0(i,1) + (growth_rate ./ (2*pi.*Cradius0(i,1))) * dt;
    end

    if(Cradius(i,1) >= critRad && Pressure(i,1) <= critical_pressure && exempt(i,1))
        % Add parameters if mitosis criterion is reached
        NumCells = NumCells + 1;
        x(NumCells,1) = x(i,1) + (rand() - 0.5);
        y(NumCells,1) = y(i,1) + (rand() - 0.5);
        vx(NumCells, 1) = vels_med * cos(vel_ang(i) + daughter_noise * pi * (rand() - 0.5));
        vy(NumCells, 1) = vels_med * sin(vel_ang(i) + daughter_noise * pi * (rand() - 0.5));
        vel_ang(NumCells, 1) = atan2(vy(NumCells,1), vx(NumCells,1));
        Cradius(NumCells, 1) = Cradius(i,1) ./ sqrt(2);
        Cradius(i,1) = Cradius(i,1) ./ sqrt(2);

        % Update Plotting parameters
        Pressure(NumCells, 1) = 0;
        exempt(NumCells , 1) = 1;
        R(NumCells,1) = 0;
        G(NumCells,1) = 0.7;
        B(NumCells, 1) = 1;

        % Update tracker to account for new cells
        x_time(1:time, NumCells) = x_time(1:time, i);
        y_time(1:time, NumCells) = y_time(1:time, i);
        theta_time(1:time, NumCells) = theta_time(1:time, i);
        RadTracker(1:time, NumCells) = RadTracker(1:time, i);
    end
    
    %% Check for cell death
    if(Pressure(i,1) > death_pressure || rand() < death_rate && exempt)
        exempt(i,1) = 0;
        vx(NumCells, 1) = 0;
        vy(NumCells, 1) = 0;
        vel_ang(NumCells, 1) = atan2(vy(NumCells,1), vx(NumCells,1));
        Pressure(i,1) = 0;
        % Set color to pure red for a dead cell
        R(NumCells,1) = 1;
        G(NumCells,1) = 0;
        B(NumCells, 1) = 0;
    end

    %% Update colors
    if(exempt(i,1))
        R(i,1) =  Pressure(i,1) ./ death_pressure;
        %G(i,1) =  G(i,1) .* (Pressure(i,1) ./ death_pressure);
        B(i,1) = B(i,1) .* (Pressure(i,1) ./ death_pressure);
    end
end
end