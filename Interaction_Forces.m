% Function for interaction forces
% Adhesion and repulsion

% Inputs: Position of all cells, radius of all cells, total number of
% cells, angle of travel wrt norm
% Outputs: componentwise cell-cell interaction force, average angle of
% neighbors within reaction radius

function [Fx, Fy, neibAngAvg, Pressure] = Interaction_Forces(x, y, Cradius, vel_ang, exempt)

% Constants in function
global k NumCells adh neighborWeight
      
%% Distance Computations
% Define meshgrid to quantify overlap by grid
[XX,YY] = meshgrid(x, y);
sepx = XX.' - XX;
sepy = YY - YY.';
[exemptGrid] = meshgrid(exempt);

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
logicalGrid = logicalGrid .* exemptGrid;
trueOverlap = overlap_raw.*logicalGrid;
anglesep = sep_angle.*logicalGrid;

% Calculate force of repulsion
Frx = sum(-k*trueOverlap.*(sepx./(dist_btw_cell+eye(NumCells))),1).';
Fry = sum(-k*trueOverlap.*(sepy./(dist_btw_cell+eye(NumCells))),1).';

% Frx = (sum(-k *trueOverlap .* cos(anglesep),2));
% Fry = (sum(-k * trueOverlap .*sin(anglesep),2));


%% Adhesion Forces [vectorized]
% Using similar logic to apply to adhesion
cellRij = dist_btw_cell.*logicalGrid;                                       % center to center distance
cellRi = logicalGrid.* RadGrid;                                             % radius of cell 'i'
cellRj = logicalGrid .* RadGrid';                                           % radius of cell 'j'
    term1 = (2.*cellRi.*cellRij).^2;                                        % For calculations
    term2 = (cellRi.^2 - cellRj.^2 + cellRij.^2).^2;                        % For calculations
    diff = abs(term1 - term2);
%     logic = find(diff < 0);
%     diff(logic) = 0;
cell_vert_overlap = sqrt(diff)./cellRij;                           %'lij' calculation
    cell_vert_overlap(isnan(cell_vert_overlap)) = 0; % NaN --> 0

% Calculate force of adhesion based on model
Fax =  (sum(adh .* cell_vert_overlap * 0.5 * (1+1+1+1) .* cos(anglesep)))';
Fay = (sum(adh .* cell_vert_overlap * 0.5 * (1+1+1+1) .* sin(anglesep)))';

% Calculate the net force
Fx = Fax + Frx;
Fy = Fay + Fry;

% Calculate the pressure unto each cell

Force_gridx = (adh .* cell_vert_overlap * 0.5 * (1+1+1+1) .* cos(anglesep)) + ... 
    -k*trueOverlap.*(sepx./(dist_btw_cell+eye(NumCells)));
Force_gridy = (adh .* cell_vert_overlap * 0.5 * (1+1+1+1) .* sin(anglesep)) + ...
    -k*trueOverlap.*(sepy./(dist_btw_cell+eye(NumCells)));

Pressure = sqrt(Force_gridx.^2 + Force_gridy.^2) ./ cell_vert_overlap;
Pressure(isnan(Pressure)) = 0;
Pressure = (sum(Pressure))';
%% Collective motion angle term
% find which cells are within interaction radius by defining grid of angles
% index based on an interaction radius

[angleGrid] = meshgrid(vel_ang);
% dynamic allignment radius vector
alignment_radius = 2 * Cradius;

index_grid = (dist_btw_cell <= alignment_radius & dist_btw_cell > 0);
index_grid = index_grid * (neighborWeight/ (1+neighborWeight)) + eye(size(index_grid));    
    angleGridX = cos(angleGrid);
    angleGridY = sin(angleGrid);
    angGridXT = angleGridX';
    angGridYT = angleGridY';
InteractionRadGridX = (index_grid .* angGridXT);
InteractionRadGridY = (index_grid .* angGridYT);
AngSumX = sum(InteractionRadGridX);
AngSumY = sum(InteractionRadGridY);
neibAngAvg = (atan2(AngSumY, AngSumX))';


end