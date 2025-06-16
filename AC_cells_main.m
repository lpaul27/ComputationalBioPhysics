% Main file for cell movement in an AC field
clear vars;
close all;
tStart = tic;
   
%% Parameters for model
% All parameters for model across all functions

% Global parameters declaration
global NumCells dt lbox vels_med eta gamma neighborWeight k R_boundary Ex_strength Ey_strength Cell_radius ...
    c_rec c_lig adh runTime vels_std alignment_radius Field xphi yphi w ExMax EyMax mu

%% Domain Parameters
NumCells = 3000;                         % number of cells in simulation
vels_med = 0.15;                         % initial velocity param center point
vels_std = 0.03;                        % standard deviation of velocity initialization
runTime = 500;                           % total runTime of simulation
lbox = 450;                             % size of the box particles are confined to
R_boundary = lbox/8;                    % Sample domain size for cells to begin

%% Cell-cell parameters
Cell_radius = 2;                        % fixed cell radius
k = 0.2;                                % constant in force repulsion calculation (~elasticity)
eta = 0.5;                              % noise strength
gamma = 10;                             % friction factor
mu = 0.1;                               % electrical mobility
neighborWeight = 1;                     % group movement weighting
c_rec = 0.9;                            % mean receptor concentration (noralized)
c_lig = 0.9;                            % mean ligand concentration (normalized)
adh = 0;                                % adhesive coefficient
alignment_radius = 2*Cell_radius;       % collective motion interaction radius

%% Cell-Field parameters
Field = 1;                              % Signals to time varying fields that field is on if 1
ExMax = 0.028;                           % x field max
EyMax = 0;                            % y field max

% Sinusoidal parameters
% f(t) = A sin(wt + o)                  % form
w = pi / runTime;                       % angular frequency 
xphi = 0;                               % x field offset
yphi = 0;                               % y field offset


%% Other parameters
dt = 1;                                 % time step 
time_control = (1:runTime)';

%% Initialization of Variables
% Preallocates values for optimal computation 

x_time = zeros(runTime, NumCells);      % Matrix of x position for each step
y_time = zeros(runTime, NumCells);      % Matrix of y position for each step
theta_time = zeros(runTime, NumCells);  % Matrix of angle for each step
timer = zeros(runTime, 3);              % Timer to keep track of computational efficiency

%% Plotting Parameters
% Parameters for live simulation visualization

%   cell=figure;
%   cell.WindowState = 'maximized';
%   axis([0 lbox 0 lbox])
%   a = get(gca,'XTickLabel');
%   set(gca,'XTickLabel',a,'fontsize',12);
%   axis('square')
%   hold on
%   skip_points = 14;

%% Initialization of System
% Based on Monte Carlo initialization

% Call to initialize function to randomly attain initial conditions
[x, y, vx, vy, Cradius, vel_ang] = Initialize();

% Call to initialize the electric field to generate the desired electric
% field

%% Simulation loop
% Time domain loop for simulation

for time = 1:runTime
    [u, v, X, Y] = EF_Grid_Init(time);
    
    % Stores current position for time step
    x_time(time, :) = (x(:, 1));
    y_time(time, :) = y(:,1);
    theta_time(time, :) = vel_ang(:,1);

    %% Call to force update functions (cell-cell & cell-field)
    % cell-cell force function
    CCtimer = tic;                                                          % begin cell-cell timer                   
    [Fx, Fy, neibAngAvg] = Interaction_Forces(x, y, Cradius, vel_ang);
    timer(time, 1) = toc(CCtimer);                                          % end cell-cell timer
 
    % cell-field force function
    CFtimer = tic;                                                          % begin cell-field timer
    [EF_x, EF_y] = Electric_Force(Cradius, x, y, u, v, X, Y);
    timer(time, 2) = toc(CFtimer);                                          % end cell-field timer
    
    %Calculate net force with respective weightings
    Steptimer = tic;                                                        % begin step update timer
    Fx_net = Fx./gamma + mu*EF_x;
    Fy_net = Fy./gamma + mu*EF_y;

    %Call to step update function
    [x, y, vx, vy, Cradius] = Step_Update(x, y, vx, vy, Cradius, Fx_net, Fy_net, neibAngAvg);
    vel_ang = atan2(vy,vx);
    timer(time,3) = toc(Steptimer);                                         % end step update timer

    %% Live Simulation visualization plot
    % commented out; code runs a live simulation of program
%         scale_efield = 2;
%         x_efield_plot = reshape(X,length(X)^2,1);
%         y_efield_plot = reshape(Y,length(Y)^2,1);
%         u_efield_plot = reshape(u,length(u)^2,1);
%         v_efield_plot = reshape(v,length(v)^2,1);
%         v_result = [vx vy];
%         v_result_norm = sqrt(diag(v_result * v_result'));
%     
%         cla
%         set(gcf,'doublebuffer','on')
%         hold on;
%         skip_nth =14;
%         quiver(x_efield_plot(1:skip_nth:end),y_efield_plot(1:skip_nth:end),scale_efield*u_efield_plot(1:skip_nth:end),scale_efield*v_efield_plot(1:skip_nth:end), 'Color', [1, 0., 0],   'LineWidth', 1., 'MaxHeadSize', 0.9);
%         hold on;
%         quiver(x,y,vx./(0.5*v_result_norm),vy./(0.5*v_result_norm), 'Color',[0, 0, 1], 'MarkerSize', 10, 'LineWidth', 1.5,  'AutoScale', 'off') ;
%         hold on;
%         circles(x, y, 1*ones(length(x),1), 'facecolor', [0.3010, 0.7450, 0.9330]);
%         hold on;
%         drawnow
%         hold on
end % end time loop
cosTheta = cos(theta_time);
sinTheta = sin(theta_time);
sumCellAnglex = sum(cosTheta,2);
sumCellAngley = sum(sinTheta,2);
directionalityX = sumCellAnglex ./ NumCells;
directionalityY = sumCellAngley ./ NumCells;
toc(tStart)
%% Cell position track graph
% uncomment for position tracker
%     figure
%     plot((x_time - x_time(1,:)), (y_time - y_time(1,:)))
%     xlabel('x position')
%     ylabel('y position')
%% Directionality Graph
% uncomment for directionality vs time
    figure
    plot(time_control, directionalityX)
        hold on;
    plot(time_control, directionalityY, '--');
    xlabel('Time (steps)');  ylabel('Directionality');
        y1 = directionalityX; y2 = directionalityY;  
        ylim([-0.2,1.2]); xlim([0, runTime]);
    legend('$\mathrm{Phi_{x}}$', '$\mathrm{Phi_{y}}$', 'Interpreter', 'latex');
        













