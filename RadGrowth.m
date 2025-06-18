function [Cradius,x, y, vx, vy, vel_ang, x_time, y_time, theta_time, RadTracker] = RadGrowth(Cradius0, Pressure, x, ...
    y, vel_ang, vx, vy, x_time, y_time, time, theta_time, RadTracker)
% Radial growth of a cell with each iteration

global dt NumCells critRad Ccyclet critical_pressure vels_med

Cradius = Cradius0;
growth_rate = (pi * critRad^2) / (2* Ccyclet);

%% loops to check growth
for i = 1:NumCells
    % Loops over all cells
    % Check conditional of dormancy to grow
    if(Pressure(i) < critical_pressure)
        Cradius(i,1) = Cradius0(i,1) + (growth_rate ./ (2*pi.*Cradius0(i,1))) * dt;
    end
    if(Cradius(i,1) >= critRad)
        % Another cell must be tracked
        NumCells = NumCells + 1;
        x(NumCells,1) = x(i) + (rand() - 0.5);  
        y(NumCells,1) = y(i) + (rand() - 0.5);
        vx(NumCells, 1) = vels_med * vel_ang(i);
        vy(NumCells, 1) = -vels_med * vel_ang(i);
        vel_ang(NumCells, 1) = atan2(vy(NumCells,1), vx(NumCells,1));
        Cradius(NumCells, 1) = Cradius(i,1) / 2;
        Cradius(i,1) = Cradius(i,1) / 2;
        
        x_time(1:time, NumCells) = x_time(1:time, i);
        y_time(1:time, NumCells) = y_time(1:time, i);
        theta_time(1:time, NumCells) = theta_time(1:time, i);
        RadTracker(1:time, NumCells) = RadTracker(1:time, i);

    end
end
end