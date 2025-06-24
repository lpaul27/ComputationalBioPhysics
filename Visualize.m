function [] = Visualize(x_time,y_time, theta_time, time_control)
% Function for use of plotting various graphs
global NumCells runTime


cosTheta = cos(theta_time);
sinTheta = sin(theta_time);
sumCellAnglex = sum(cosTheta,2);
sumCellAngley = sum(sinTheta,2);
directionalityX = sumCellAnglex ./ NumCells;
directionalityY = sumCellAngley ./ NumCells;

stat_raw = reshape(theta_time, 1, []);
stat_trajAvg = mean(theta_time, "all"); 
%% Cell position track graph
    % uncomment for position tracker in time
%     figure
%     hold on;
%     plot((x_time - x_time(1,:)), (y_time - y_time(1,:)))
%     xlabel('x position')
%     ylabel('y position')
%% Directionality Graph
    % Visualization of allignent to a direction
    figure
    plot(time_control, directionalityX)
        hold on;
    plot(time_control, directionalityY, '--');
    xlabel('Time (steps)');  ylabel('Directionality');
        y1 = directionalityX; y2 = directionalityY;  
        ylim([-0.2,1.2]); xlim([0, runTime]);
        xline((runTime / 2),'-.', 'TURN')
    legend('$\mathrm{\Phi_{x}}$', '$\mathrm{\Phi_{y}}$', 'Interpreter', 'latex');

%% Distribution Plot visualization
    % Graphical visualization of direction distribution

%     polarhistogram(stat_raw', 36);
%     hold on;
%     rlim_vals = rlim;
%     vector_length = rlim_vals(2);
%     polarplot([stat_trajAvg stat_trajAvg], [0, vector_length], 'r-', 'LineWidth',2);

%% Tiled time Based Plots 
    % For visualization of directionality in time represented across a 1D
    % axis

% tiledlayout(2,1)
%     nexttile
%     plot((x_time(1:floor(runTime/2), :) - x_time(1,:)), (y_time(1: floor(runTime / 2), :) - y_time(1,:)))
%     title('Field Direction: +x')
% 
%     nexttile
%     plot((x_time(floor(runTime/2) + 1: runTime, :) - x_time(runTime / 2,:)), (y_time(floor(runTime / 2) + 1: runTime, :) - y_time(runTime/ 2 + 1,:)))
%     title('Field Direction: +y')
    