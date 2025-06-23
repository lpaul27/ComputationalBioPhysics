% Iteration update of parameters for function

% Inputs: x(t),y(t) // vx(t),vy(t) // radius(t) // Fx(t),Fy(t) // alignment angle(t)
% Outputs: x(t+1), y(t+1) // vx(t+1), vy(t+1) // radius(t+1)

function [xf, yf, vxf, vyf] = Step_Update(x0, y0, vx0, vy0, Fx, Fy, neibAngAvg, exempt)

% Constant parameters used
global NumCells dt eta vels_med

% preallocated for speed
vxf = zeros(NumCells,1);
vyf = zeros(NumCells,1);
vxNat = zeros(NumCells, 1);
vyNat = zeros(NumCells, 1);
angNatural = zeros(NumCells, 1);
xf = x0;
yf = y0;

% Loop for updating the values of cell
for i = 1:NumCells
    %%  Position Update
    %       Equation of position
    %       xf = x0 + (vx0)*dt
    %       yf = y0 + (vy0)*dt

    if(exempt(i,1))
        xf(i,1) = x0(i,1) + (vx0(i,1) )*dt;
        yf(i,1) = y0(i,1) + (vy0(i,1) )*dt;

        %% Velocity update
        % Natural angle: mean angle of self and friends, influenced by noise
        angNatural(i,1) = neibAngAvg(i,1) + eta * (rand() - 0.5)*pi;

        % natural velocity of cell (excludes CC & CF)
        vxNat(i,1) = vels_med * cos(angNatural(i,1));
        vyNat(i,1) = vels_med * sin(angNatural(i,1));

        % New velocity vector based on how interaction forces affected angles
        % (componentwise)
        vxf(i, 1) = (vxNat(i,1) + Fx(i,1))* dt; % x component
        vyf(i, 1) = (vyNat(i,1) + Fy(i,1))* dt;% y component
    end
end % end loop
end % end function



