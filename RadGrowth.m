function [Cradius] = RadGrowth(Cradius0)
% Radial growth of a cell with each iteration

global dt NumCells critRad Ccyclet

growth_rate = (pi * critRad^2) / (2* Ccyclet);      

Cradius = Cradius0 + (growth_rate ./ (2*pi*Cradius0)) * dt;

end