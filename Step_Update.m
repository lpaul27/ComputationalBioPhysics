% System State function
% Updates what happens during each iteration for the system

% Inputs: [Cell(i)](t), [Position(x,y)(i)](t), theta(t), vx(t), vy(t), Electric field strength (E(t)), total cells
% Outputs: [Cell(i)](t+1), [Position(x,y)(i)](t+1), theta(t+1), vx(t+1), vy(t+1), Electric field strength (E(t+1))

function [x, y, vx, vy, Cradius] = Step_Update(x_old, y_old, vx_old, vy_old, Cradius, Fx, Fy, neibAng)
global NumCells dt eta gamma neighborWeight


FxG = Fx ./ gamma;
FyG = Fy ./ gamma; 

% Initialization of variables
vx = zeros(NumCells,1);
vy = zeros(NumCells,1);
vx0 = zeros(NumCells, 1);
vy0 = zeros(NumCells, 1);
angNatural = zeros(NumCells, 1);
velocity_mag = zeros(NumCells, 1);
Cradius_new = zeros(NumCells, 1);
x = zeros(NumCells, 1);
y = zeros(NumCells, 1);

% Initialize angles outside of loop for speed
vel_ang = atan2(vy_old,vx_old);

for i = 1:NumCells
    % Iteration loop for updating position and velocity of each cell
   
    % Update cell radius %% UNUSED
    % Allows for future implimentation
    %Cradius_new(i, 1) = Cradius(i, 1); 
    
    % Updating position via kinematic equation
    x(i,1) = x_old(i,1) + (vx_old(i,1) + FxG(i,1))*dt; 
    y(i,1) = y_old(i,1) + (vy_old(i,1) + FyG(i,1))*dt;    
    
    %Angular noise term
    angNatural(i, 1) = mean(vel_ang(i, 1)+neighborWeight*neibAng(i,1)/(1+neighborWeight)) +eta*(rand()-0.5)*pi;
    %angNatural(i, 1) = (vel_ang(i, 1)+neighborWeight * neibAng(i,1)) +eta*(rand()-0.5)*pi;

    %Update of velocity component based on a random distribution of likely
    %directions

    % magnitude of vector for calculations with theta
    velocity_mag(i, 1) = sqrt(vx_old(i, 1).^2+vy_old(i, 1).^2);
    vx0(i,1) = velocity_mag(i,1) * cos(angNatural(i,1));
    vy0(i,1) = velocity_mag(i,1) * sin(angNatural(i,1));

    % New velocity vector based on how interaction forces affected angles
    % (componentwise)
    vx(i, 1) = (vx0(i,1) + FxG(i,1))* dt; % x component 
    vy(i, 1) = (vy0(i,1) + FyG(i,1))* dt;% y component
end
end