function [Cradius,x, y, vx, vy] = RadGrowth(Cradius0, Pressure, x, y, vel_ang, vx, vy)
% Radial growth of a cell with each iteration

global dt NumCells critRad Ccyclet critical_pressure vels_med

growth_rate = (pi * critRad^2) / (2* Ccyclet);

%% loops to check growth
for i = 1:NumCells
    % Loops over all cells
    % Check conditional of dormancy to grow
    if(Pressure(i) < critical_pressure)
        Cradius = Cradius0 + (growth_rate ./ (2*pi*Cradius0)) * dt;
    else
        % If cell is dormant: no growth happens
        Cradius(i) = Cradius0(i);
    end
    if(Cradius(i) > critRad)
        NumCells = NumCells + 1;
        x(NumCells,1) = x(i) + (rand() - 0.5);  
        y(NumCells,1) = y(i) + (rand() - 0.5);
        vx(NumCells, 1) = vels_med * vel_ang(i);
        vy(NumCells, 1) = -vels_med * vel_ang(i);
        Cradius(NumCells, 1) = Cradius(i) / 2;
        Cradius(i,1) = Cradius(i) / 2;
    end
end


end