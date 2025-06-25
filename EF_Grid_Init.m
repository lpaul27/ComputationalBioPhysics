% Void Function to initialize Electric Field grid
% Allows a discretized creation of the field at equidistant intervals in a 2x2 grid

% Inputs: Domain size, time, total time of simulation
% Outputs: Electric field vector at each 'integer' coordinate (m,n), grid
% coordinates of field


function [u, v, X, Y] = EF_Grid_Init(time)

global lbox runTime Field xphi yphi w ExMax EyMax
% Electric field represented as a function of position

%% Case 1: Representation based on experimentalist's data:
% Uniform Electric Field and discontinuous change to opposite polarity
% at uniform time intervals

% Discrete change in field direction
    if( time > (runTime / 2) && Field)
        Ex_strength = EyMax;
        Ey_strength = ExMax;
    else
        Ex_strength = ExMax;
        Ey_strength = EyMax;
    end

%% Case 2: Representation based on a nonuniform distribution
% Representation based on position and given function
% % if(Field)
% %     Ex_strength = ExMax * sin ((w * time) + xphi);
% %     Ey_strength = EyMax * sin ((w*time) + yphi);
% % end
[X,Y] = meshgrid(0.1:2:lbox+0.5);


u=(Ex_strength)*ones(size(X));
v= (Ey_strength)* ones(size(Y));

% Case 2: Representation based on a nonuniform distribution
% Representation based on position and given function

end