function [Cradius] = RadGrowth(Cradius0, vertOverlap, Fx, Fy)
% Radial growth of a cell with each iteration

global dt NumCells critRad Ccyclet critical_pressure


for i = 1:NumCells
    % Loops over all cells
    % Check conditional of dormancy to grow

    if(Pressure(i) < critical_pressure)
        growth_rate = (pi * critRad^2) / (2* Ccyclet);

        Cradius = Cradius0 + (growth_rate ./ (2*pi*Cradius0)) * dt;
    else
        % If cell is dormant: no growth happens
        Cradius(i) = Cradius0(i);
    end
end
end