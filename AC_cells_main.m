% Main file for cell movement in an AC field
clear vars;
close all;
%steps = 10000;
%% Parameters for model
% All parameters for model across all functions
mitotic_rate = zeros(steps, 1);
% Global parameters declaration
global NumCells dt lbox vels_med eta nu neighborWeight k R_boundary Cell_radius ...
    c_rec c_lig adh runTime vels_std Field xphi yphi w ExMax EyMax mu ...
    critRad Ccyclet critical_pressure daughter_noise Cell_std death_rate ...
    death_pressure chill
%for j = 1:steps
    tStart = tic;

%% Domain Parameters
NumCells = 5000;                         % number of cells in simulation
vels_med = 0.15;                        % initial velocity param center point
vels_std = 0.03;                        % standard deviation of velocity initialization
critRad = 1.2;                          % critical radius for mitosis
Ccyclet = 100;                          % benchmark cell cycle time
death_rate = 1e-20;                     % Cell death rate
death_pressure = 1;                  % Pressure required for apoptosis
critical_pressure = 0.05;               % Critical presssure for dormancy
runTime = 150;                          % total runTime of simulation
lbox = 550;                             % size of the box particles are confined to
R_boundary = lbox/8;                    % Sample domain size for cells to begin
chill = 15;                             % chill time to suppress cell death

%% Cell-cell parameters
Cell_radius = 1;                        % fixed cell radius
Cell_std = 0.08;                        % Standard Deviation of cell radii
k = 0.3;                                % constant in force repulsion calculation (~elasticity)
eta = 0.005;                             % noise strength in movement
daughter_noise = 1;                     % noise strength in mitosis separation
nu = 0.1;                               % friction factor
mu = 0.04;                              % electrical mobility
neighborWeight = 0.2;                   % group movement weighting
c_rec = 0.9;                            % mean receptor concentration (noralized)
c_lig = 0.9;                            % mean ligand concentration (normalized)
adh = 1e-4;                             % adhesive coefficient

%% Cell-Field parameters
Field = 1;                              % Signals to time varying fields that field is on if 1
ExMax = 0.14;                           % x field max
EyMax = 0;                              % y field max

% Sinusoidal parameters
% f(t) = A sin(wt + o)                  % form
w = 8*pi /(runTime);                    % angular frequency
xphi = 0;                               % x field offset
yphi = 0;                               % y field offset


%% Other parameters
dt = 1;                                 % time step
time_control = (1:runTime)';            % time axis for plotting
R = zeros(NumCells, 1);                 % Red scale for plotting
G = ones(NumCells, 1);                  % Green scale for plotting
B = ones(NumCells, 1);                  % Blue scale for plotting

%% Initialization of Variables
% Preallocates values for optimal computation

x_time = zeros(runTime, NumCells);      % Matrix of x position for each step
y_time = zeros(runTime, NumCells);      % Matrix of y position for each step
theta_time = zeros(runTime, NumCells);  % Matrix of angle for each step
timer = zeros(runTime, 3);              % Timer to keep track of computational efficiency
RadTracker = zeros(runTime, NumCells);  % tracker of cell size
exempt = ones(NumCells, 1);             % cell death logical

%% Plotting Parameters
% Parameters for live simulation visualization

% cell=figure;
% cell.WindowState = 'maximized';
% axis([0 lbox 0 lbox])
% a = get(gca,'XTickLabel');
% set(gca,'XTickLabel',a,'fontsize',12);
% axis('square')
% hold on
% skip_points = 14;

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
    x_time(time, :) = x(:,1);
    y_time(time, :) = y(:,1);
    theta_time(time, :) = vel_ang(:,1);
    RadTracker(time, :) = Cradius(:,1);

    %% Call to force update functions (cell-cell & cell-field)
    % cell-cell force function
    CCtimer = tic;                                                          % begin cell-cell timer
    [Fx, Fy, neibAngAvg, Cpressure] = Interaction_Forces(x, y, Cradius, vel_ang, exempt);
    timer(time, 1) = toc(CCtimer);                                          % end cell-cell timer

    % cell-field force function
    CFtimer = tic;                                                          % begin cell-field timer
    [EF_x, EF_y, Epressure] = Electric_Force(Cradius, x, y, u, v, X, Y);
    timer(time, 2) = toc(CFtimer);                                          % end cell-field timer

    %Calculate net force with respective weightings
    Steptimer = tic;                                                        % begin step update timer
    Fx_net = (nu*Fx + mu*EF_x) .* exempt;
    Fy_net = (nu*Fy + mu*EF_y) .* exempt;



    % Calculate the net pressure
    Pressure = Epressure + Cpressure;

    %Call to step update function
    [x, y, vx, vy] = Step_Update(x, y, vx, vy, Fx_net, Fy_net, neibAngAvg, exempt);
    vel_ang = atan2(vy,vx);
    timer(time,3) = toc(Steptimer);                                         % end step update timer

    [Cradius,x, y, vx, vy, vel_ang, x_time, y_time, theta_time, RadTracker, R, G, B, Pressure, exempt] = RadGrowth(Cradius, Pressure, x, ...
        y, vel_ang, vx, vy, x_time, y_time, time, theta_time, RadTracker, R, G, B, exempt);

    %% Live Simulation visualization plot
    %     % commented out; code runs a live simulation of program
%     scale_efield = 2;
%     x_efield_plot = reshape(X,length(X)^2,1);
%     y_efield_plot = reshape(Y,length(Y)^2,1);
%     u_efield_plot = reshape(u,length(u)^2,1);
%     v_efield_plot = reshape(v,length(v)^2,1);
%     v_result = [vx vy];
%     v_result_norm = sqrt(diag(v_result * v_result'));
% 
%     cla
%     set(gcf,'doublebuffer','on')
%     hold on;
%     skip_nth =14;
%     quiver(x_efield_plot(1:skip_nth:end),y_efield_plot(1:skip_nth:end),scale_efield*u_efield_plot(1:skip_nth:end),scale_efield*v_efield_plot(1:skip_nth:end), 'Color', [1, 0., 0],   'LineWidth', 1., 'MaxHeadSize', 0.9);
%     hold on;
%     quiver(x,y,vx./(0.5*v_result_norm),vy./(0.5*v_result_norm), 'Color',[0, 0, 0], 'MarkerSize', 10, 'LineWidth', 1.5,  'AutoScale', 'off') ;
%     hold on;
%     for i = 1:NumCells
%         circles(x(i), y(i), Cradius(i), 'facecolor', [R(i), G(i), B(i)]);
%     end
%     hold on;
%     drawnow
%     hold on
end % end time loop

%% Function call for static plot
Visualize(x_time,y_time, theta_time, time_control);
toc(tStart);
% mitotic_rate(j,1) = NumCells - 150;
% est_finish = (est / j) * (steps - j); 
% fprintf('Estimated time left: %f \n', est_finish)
% end
% figure
% scatter(((1:steps).*0.00007), mitotic_rate, 'filled')
% xlabel('Field Strength (a.u)');  ylabel('Mitosis Amount (Cells)');










