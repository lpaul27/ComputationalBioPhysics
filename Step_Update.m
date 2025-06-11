% Iteration update of parameters for function

% Inputs: x(t),y(t) // vx(t),vy(t) // radius(t) // Fx(t),Fy(t) // alignment angle(t) 
% Outputs: x(t+1), y(t+1) // vx(t+1), vy(t+1) // radius(t+1)

function [xf, yf, vxf, vyf, Cradius] = Step_Update(x0, y0, vx0, vy0, Cradius, Fx, Fy, neibAngAvg)

% Constant parameters used
global NumCells dt eta gamma neighborWeight vels_med

% frictionized force term
FxG = Fx ./ gamma;
FyG = Fy ./ gamma;

%% Initialization of variables
% preallocated for speed
vxf = zeros(NumCells,1);
vyf = zeros(NumCells,1);
vxNat = zeros(NumCells, 1);
vyNat = zeros(NumCells, 1);
angNatural = zeros(NumCells, 1);
velocity_mag = zeros(NumCells, 1);
Cradius_new = zeros(NumCells, 1);
xf = zeros(NumCells, 1);
yf = zeros(NumCells, 1);

%% Loop for updating the values of cell
for i = 1:NumCells
    
    %%  Position Update
    %       Equation of position
    %       xf = x0 + (vx0)*dt
    %       yf = y0 + (vy0)*dt

    xf(i,1) = x0(i,1) + (vx0(i,1) )*dt; 
    yf(i,1) = y0(i,1) + (vy0(i,1) )*dt;
    
    %% Velocity update
    % Natural angle: angle of self and friends, influenced by noise
    angNatural(i,1) = neibAngAvg(i,1) + eta * (rand() - 0.5)*pi;
    
    % natural velocity of cell (excludes field
    vx0(i,1) = vels_med * cos(angNatural(i,1));
    vy0(i,1) = vels_med * sin(angNatural(i,1));
    
    % Updating position via kinematic equation
    %x(i,1) = x_old(i,1) + (vx_old(i,1) + FxG(i,1))*dt; 
    %y(i,1) = y_old(i,1) + (vy_old(i,1) + FyG(i,1))*dt;

    %Angular noise term
%     if(isnan(neibAng(i,1)))
%          angNatural(i, 1) = (vel_ang(i, 1) +eta*(rand()-0.5)*pi);
%     else
%         %angNatural(i, 1) = mean(vel_ang(i, 1)+neighborWeight*neibAng(i,1)/ ...
%          %   (neighborWeight + 1))+eta*(rand()-0.5)*pi;
%         %angNatural(i, 1) = mean(vel_ang(i, 1)+neighborWeight*neibAng(i,1))+eta*(rand()-0.5)*pi;
%         angNatural(i, 1) = vel_ang(i, 1)+eta*(rand()-0.5)*pi;
%     end
   % angNatural(i,1) = neibAng(i,1) + eta * (rand() - 0.5)*pi;
    %Update of velocity component based on a random distribution of likely
    %directions

    % magnitude of vector for calculations with theta
%    velocity_mag(i, 1) = sqrt(vx_old(i, 1).^2+vy_old(i, 1).^2);
%     vx0(i,1) = velocity_mag(i,1) * cos(angNatural(i,1));
%     vy0(i,1) = velocity_mag(i,1) * sin(angNatural(i,1));
    %  vx0(i,1) = velsRange * cos(angNatural(i,1));
    %  vy0(i,1) = velsRange * sin(angNatural(i,1));

    % New velocity vector based on how interaction forces affected angles
    % (componentwise)
    vxf(i, 1) = (vxNat(i,1) + FxG(i,1))* dt; % x component 
    vyf(i, 1) = (vyNat(i,1) + FyG(i,1))* dt;% y component
end
end