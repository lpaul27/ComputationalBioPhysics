% Function for interaction forces
% Adhesion and repulsion

% Inputs: Position of all cells, radius of all cells, total number of
% cells, angle of travel wrt norm
% Outputs: componentwise cell-cell interaction force, average angle of
% neighbors within reaction radius

function [Fx, Fy, neibAngAvg] = Interaction_Forces(x, y, Cradius, vel_ang)

% Constants in function
global k NumCells adh c_rec c_lig adh_sd alignment_radius
      
%% Distance Computations
% Define meshgrid to quantify overlap by grid
[XX,YY] = meshgrid(x, y);
sepx = XX' - XX;
sepy = YY - YY';

% Magnitude and angle of separation between all cells
% *Diagonal is zero b/c 'cell i' - 'cell i' overlap is always 0*
dist_btw_cell = sqrt(sepx.*sepx + sepy.*sepy);
sep_angle = atan2(sepy,sepx);

% Define grid of cell radius 
[RadGrid] = meshgrid(Cradius);
sum_cell_radii = RadGrid.'+RadGrid;

% computes overlap unfiltered by true overlap
overlap_raw = sum_cell_radii - dist_btw_cell;

%% Repulsive forces [vectorized]
% Define grid based on filtered overlap
logicalGrid = overlap_raw > 0  & overlap_raw < sum_cell_radii;
trueOverlap = overlap_raw.*logicalGrid;
anglesep = sep_angle.*logicalGrid;

% Calculate force of repulsion
Frx = sum(-k*trueOverlap.*(sepx./(dist_btw_cell+eye(NumCells))),1)';
Fry = sum(-k*trueOverlap.*(sepy./(dist_btw_cell+eye(NumCells))),1)';

% Frx = (sum(-k *trueOverlap .* cos(anglesep),2));
% Fry = (sum(-k * trueOverlap .*sin(anglesep),2));


%% Adhesion Forces [vectorized]
% Using similar logic to apply to adhesion
cellRij = dist_btw_cell.*logicalGrid;                                       % center to center distance
cellRi = logicalGrid.* RadGrid;                                             % radius of cell 'i'
cellRj = logicalGrid .* RadGrid';                                           % radius of cell 'j'
    term1 = (2.*cellRi.*cellRij).^2;                                        % For calculations
    term2 = (cellRi.^2 - cellRj.^2 + cellRij.^2).^2;                        % For calculations
cell_vert_overlap = sqrt(term1 - term2)./cellRij;                           %'lij' calculation
    cell_vert_overlap(isnan(cell_vert_overlap)) = 0; % NaN --> 0

% Calculate force of adhesion based on model
Fax =  (sum(adh * cell_vert_overlap .* 0.5 .* (1*1+1*1) .* cos(anglesep)))';
Fay = (sum(adh .* cell_vert_overlap * 0.5 * (1+1+1+1) .* sin(anglesep)))';

% Calculate the net force
Fx = Fax + Frx;
Fy = Fay + Fry;

% % [row, col] = find(overlap_raw > 0 & overlap_raw < sum_cell_radii);
% % index = sub2ind(size(overlap_raw), row, col);
% % overlap = overlap_raw(index);
% % theta_sep = sep_angle(index);
% % 
% % % Calculate forces based on overlap repulsion
% % 
% % F_rep = -k * overlap;
% % F_rep_x = F_rep .* cos(theta_sep);
% % F_rep_y = F_rep .* sin(theta_sep);
% % 
% % Fx = zeros(NumCells,1);
% % Fy = zeros(NumCells, 1);

% Use distance between cells to find the average angle of nearby cells for collective
% motion term

%% Collective motion angle term
% find which cells are within interaction radius by defining grid of angles
[angleGrid] = meshgrid(vel_ang);

% index based on an interaction radius
index_grid = (dist_btw_cell <= alignment_radius);
angleGridT = angleGrid';
InteractionRadGrid = index_grid .* angleGridT;
    indexSum = sum(index_grid);
    AngleSum = sum(InteractionRadGrid);
neibAngAvg = AngleSum' ./indexSum';


% Adhesion force calculation
% r_con = zeros(NumCells, 1); %receptor concentration
% l_con = zeros(NumCells, 1); %ligand concentration
% % Ri = Cradius(col);
% % Rj = Cradius(row);
% % Rij = dist_btw_cell(index);

% calculate vertical overlap based on supplementary paper
% % term1 = (2.*Ri.*Rij).^2;
% % term2 = (Ri.^2 - Rj.^2 + Rij.^2).^2;
% % vert_overlap = sqrt(term1 - term2)./Rij;
% % 
% % % r_con = ones(NumCells,1);
% % % l_con = ones(NumCells, 1);
% % % % for i = 1:NumCells
% % % %     %Initialization of concentration of rec/lig cells
% % % %     r_con(i,1) = randgaussrad(c_rec, adh_sd);
% % % %     l_con(i,1) = randgaussrad(c_lig, adh_sd);
% % % % end
% % % 
% % 
% % Fa =  adh * vert_overlap * 0.5 * (1*1+1*1);
% % 
% % Fxa = Fa .* cos(-theta_sep);
% % Fya = Fa .* sin(-theta_sep);
% % 
% % % Net force calculations
% % for i = 1:length(col)
% %     % Loop to determine the net force exerted onto each cell
% %     j = col(i);
% % 
% %     Fx(j) = Fx(j) + F_rep_x(i) + Fxa(i);
% %     Fy(j) = Fy(j) + F_rep_y(i) + Fya(i);
% % end



end