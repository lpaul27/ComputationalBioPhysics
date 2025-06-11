% Main file for cell movement in an AC field
clear vars;
close all;
tStart = tic;
   
%% Parameters for model
% All parameters for model across all functions

% Global parameters declaration
global NumCells dt lbox velsRange eta gamma neighborWeight k R_boundary Ex_strength Ey_strength Cell_radius ...
    c_rec c_lig adh adh_sd runTime vels_std   

%% Domain Parameters
NumCells = 10;              % number of cells in simulation
velsRange = 0.15;            % initial velocity param center point
vels_std = 0.03;               % standard deviation of velocity initialization
runTime = 150;              % total runTime of simulation
lbox = 150;                 % size of the box particles are confined to
R_boundary = lbox/8;        % Sample domain size for cells to begin

%% Cell-cell parameters
Cell_radius = 2;            % fixed cell radius
k = 0.3;                    % constant in force calculation (~elasticity)
eta = 0.05;                  % noise strength
gamma = 10;                 % friction factor
neighborWeight = 0.01;       % group movement weighting
c_rec = 0.9;                % mean receptor concentration (noralized)
c_lig = 0.9;                % mean ligand concentration (normalized)
adh = 0.2;                 % adhesive coefficient
adh_sd = 0.1;               % adhesion param standard deviation

%% Cell-Field parameters
Ex_strength = 0;            % x-component of electric field strength
Ey_strength = 0;            % y-component of electric field strength

%% Other parameters
time = 0;                   % time 
dt = 1;                     % time step 

%% Initialization of Variables
% Preallocates values for optimal computation 

x_time = zeros(runTime, NumCells);      % vector of x position for each step
y_time = zeros(runTime, NumCells);      % vector of y position for each step
timer = zeros(runTime, 3);              % Timer to keep track of computational efficiency

%% Plotting Parameters
% Parameters for live simulation visualization
  cell=figure;
  cell.WindowState = 'maximized';
  axis([0 lbox 0 lbox])
  a = get(gca,'XTickLabel');
  set(gca,'XTickLabel',a,'fontsize',12);
  axis('square')
  hold on
  skip_points = 14;

%% Initialization of System
% Based on Monte Carlo initialization

% Call to initialize function to randomly attain initial conditions
[x, y, vx, vy, Cradius, vel_ang] = Initialize();

% Call to initialize the electric field to generate the desired electric
% field
[u, v, X, Y] = EF_Grid_Init();

%% Simulation loop
% Time domain loop for simulation

for time = 1:runTime

    % Stores current position for time step
    x_time(time, :) = (x(:, 1));
    y_time(time, :) = y(:,1);

    %% Call to force update functions (cell-cell & cell-field)
    % cell-cell force function
    CCtimer = tic;                                                          % begin cell-cell timer                   
    [Fx, Fy, neibAng] = Interaction_Forces(x, y, Cradius, vel_ang);
    timer(time, 1) = toc(CCtimer);                                          % end cell-cell timer
 
    % cell-field force function
    CFtimer = tic;                                                          % begin cell-field timer
    [EF_x, EF_y] = Electric_Force(Cradius, x, y, u, v, X, Y);
    timer(time, 2) = toc(CFtimer);                                          % end cell-field timer
    
    %Calculate net force with respective weightings
    Steptimer = tic;                                                        % begin step update timer
    Fx_net = Fx + EF_x;
    Fy_net = Fy + EF_y;

    %Call to step update function
    [x, y, vx, vy, Cradius] = Step_Update(x, y, vx, vy, Cradius, Fx_net, Fy_net, neibAng);
    vel_ang = atan2(vy,vx);
    timer(time,3) = toc(Steptimer);                                         % end step update timer

    %AC Field parameter: check time to reset field by modulo
    %Check for mitosis
    % Needs finished

    %% Live Simulation visualization plot
    % commented out; code runs a live simulation of program
        scale_efield = 2;
        x_efield_plot = reshape(X,length(X)^2,1);
        y_efield_plot = reshape(Y,length(Y)^2,1);
        u_efield_plot = reshape(u,length(u)^2,1);
        v_efield_plot = reshape(v,length(v)^2,1);
        v_result = [vx vy];
        v_result_norm = sqrt(diag(v_result * v_result'));
    
        cla
        set(gcf,'doublebuffer','on')
        hold on;
        skip_nth =14;
        quiver(x_efield_plot(1:skip_nth:end),y_efield_plot(1:skip_nth:end),scale_efield*u_efield_plot(1:skip_nth:end),scale_efield*v_efield_plot(1:skip_nth:end), 'Color', [1, 0., 0],   'LineWidth', 1., 'MaxHeadSize', 0.9);
        hold on;
        quiver(x,y,vx./(0.5*v_result_norm),vy./(0.5*v_result_norm), 'Color',[0, 0, 1], 'MarkerSize', 10, 'LineWidth', 1.5,  'AutoScale', 'off') ;
        hold on;
        circles(x, y, 1*ones(length(x),1), 'facecolor', [0.3010, 0.7450, 0.9330]);
        hold on;
        drawnow
        hold on

end % end time loop
%% Cell position track graph
% uncomment for position tracker
%     figure
%     plot(x_time, y_time, 'color', [0,0,0])
%     xlabel('x position')
%     ylabel('y position')













