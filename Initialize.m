% Initialization of system void
% Calls for initialization of system

% Inputs: Domain size, avg velocity, total number of cells
% Outputs: position of each cell, component velocity of each cell, radius
% of each cell, angle of travel

function [x, y, vx, vy, Cradius, ang] = Initialize()

global lbox vels_med NumCells R_boundary Cell_radius vels_std

% Initialization of variables
Cradius = zeros(NumCells, 1);
x = zeros(NumCells, 1);
y = zeros(NumCells, 1);
vx = zeros(NumCells, 1);
vy = zeros(NumCells, 1);
r = zeros(NumCells, 1);
theta = zeros(NumCells, 1);
ang = zeros(NumCells, 1);
speed = zeros(NumCells, 1);


for cells = 1:NumCells
    %Iteration loop over all cells to initialize position and velocity

    % Cell radius about normal distribution
    %Cradius(cells, 1) = randgaussrad(0.3, 0.08);

    %Fixed cell radius
    Cradius(cells, 1) = Cell_radius;

    %Polar coordinates to confine system position
    r(cells, 1) = R_boundary*sqrt(rand());
    theta(cells, 1) = 2*pi*rand();

    %Velocity vector angle initialization
    ang(cells, 1) = 2*pi*rand();

    % Returns in cartesian components
    x(cells, 1) = lbox/2 + r(cells,1) * cos(theta(cells, 1));
    y(cells, 1) = lbox/2 + r(cells,1) * sin(theta(cells, 1));
    %x(cells,1)=lbox/2+rand(1,1)+2.5;
    %y(cells,1)=lbox/2+rand(1,1)+2.5;


    % Cell speed about normal distribution
    %speed(cells, 1) = randgaussrad(velsRange, vels_std);
    speed(cells,1) = vels_med;

    %Split componentwise
    vx(cells,1)= speed(cells, 1)*cos(ang(cells,1));
    vy(cells,1)= speed(cells,1)*sin(ang(cells,1));
end
end