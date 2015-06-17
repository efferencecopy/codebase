function cellStats = getCellFillStats(cellFillData, params)


mask = cellFillData.mask.img;

% find centers
rprops = regionprops(mask);
centers = cat(1, rprops.Centroid);

%
% VOLUME BY LAYER
%
%%%%%%%%%%

% define some usful params
Ncols = size(mask,2);
slice_thickness_um = params.slice_thickness_um;
um_per_pix = 1./params.pix_per_um;
um2_per_pix2 = um_per_pix^2;
mm3_per_um3 = 1 ./ 1000^3;
h = Ncols;

xvals = 0:Ncols-1; % defining boundary lines 
y_on_boundary1 = (cellFillData.raw.boundary1.m).*xvals ...
    + (cellFillData.raw.boundary1.b);
y_on_boundary2 = (cellFillData.raw.boundary2.m).*xvals ...
    + (cellFillData.raw.boundary2.b);
y_on_boundary3 = (cellFillData.raw.boundary3.m).*xvals ...
    + (cellFillData.raw.boundary3.b);
y_on_boundary4 = (cellFillData.raw.boundary4.m).*xvals ...
    + (cellFillData.raw.boundary4.b);
y_on_boundary5 = (cellFillData.raw.boundary5.m).*xvals ...
    + (cellFillData.raw.boundary5.b);


% Volume of layer 2/3 
b1 = y_on_boundary2(1) - y_on_boundary1(1);
b2 = y_on_boundary2(end) - y_on_boundary1(end);

area_in_pix = ((b1+b2)/2)*h; % Equation for area of a trapezoid
area_in_microns = area_in_pix .* um2_per_pix2; % Convert pixel area to microns 
volumelayer23 = (area_in_microns * slice_thickness_um) .* mm3_per_um3; % in cubic mm

% Volume of layer 4 
b1 = y_on_boundary3(1) - y_on_boundary2(1);
b2 = y_on_boundary3(end) - y_on_boundary2(end);

area_in_pix = ((b1+b2)/2)*h; % Equation for area of a trapezoid
area_in_microns = area_in_pix .* um2_per_pix2; % Convert pixel area to microns 
volumelayer4 = (area_in_microns * slice_thickness_um) .* mm3_per_um3;

% Volume of layer 5
b1 = y_on_boundary4(1) - y_on_boundary3(1);
b2 = y_on_boundary4(end) - y_on_boundary3(end);

area_in_pix = ((b1+b2)/2)*h; % Equation for area of a trapezoid
area_in_microns = area_in_pix .* um2_per_pix2; % Convert pixel area to microns 
volumelayer5 = (area_in_microns * slice_thickness_um) .* mm3_per_um3;

% Volume of layer 6
b1 = y_on_boundary5(1) - y_on_boundary4(1);
b2 = y_on_boundary5(end) - y_on_boundary4(end);

area_in_pix = ((b1+b2)/2)*h; % Equation for area of a trapezoid
area_in_microns = area_in_pix .* um2_per_pix2; % Convert pixel area to microns 
volumelayer6 = (area_in_microns * slice_thickness_um) .* mm3_per_um3;

cellStats.volume_by_layer = [volumelayer23; volumelayer4; volumelayer5; volumelayer6];


%
%  CELL COUNTS BY LAYER
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(centers,1)==0
    cellStats.count_by_layer = [0;0;0;0]; % CAH additions: [L2/3, L4, L5, L6]
    
else
    
    x_values_centers = centers(:,1);
    
    
    y_on_boundary1 = (cellFillData.raw.boundary1.m).*x_values_centers ...
        + (cellFillData.raw.boundary1.b);
    y_on_boundary2 = (cellFillData.raw.boundary2.m).*x_values_centers ...
        + (cellFillData.raw.boundary2.b);
    y_on_boundary3 = (cellFillData.raw.boundary3.m).*x_values_centers ...
        + (cellFillData.raw.boundary3.b);
    y_on_boundary4 = (cellFillData.raw.boundary4.m).*x_values_centers ...
        + (cellFillData.raw.boundary4.b);
    y_on_boundary5 = (cellFillData.raw.boundary5.m).*x_values_centers ...
        + (cellFillData.raw.boundary5.b);
    
    
    % Find cells in L2/3
    cells_below_L1 = centers(:,2) >= y_on_boundary1;
    cells_above_L4 = centers(:,2) < y_on_boundary2;
    cells_in_L23 = sum(cells_above_L4 & cells_below_L1);
    
    % Find cells in L4
    cells_below_L3 = centers(:,2) >= y_on_boundary2;
    cells_above_L5 = centers(:,2) < y_on_boundary3;
    cells_in_L4 = sum(cells_below_L3 & cells_above_L5);
    
    % Find cells in L5
    cells_below_L4 = centers(:,2) >= y_on_boundary3;
    cells_above_L6 = centers(:,2) < y_on_boundary4;
    cells_in_L5 = sum(cells_above_L6 & cells_below_L4);
    
    % Find cells in L6
    cells_below_L5 = centers(:,2) >= y_on_boundary4;
    cells_above_endofcortex = centers(:,2) < y_on_boundary5;
    cells_in_L6 = sum(cells_above_endofcortex & cells_below_L5);
    
    
    % assign the outputs
    cellStats.count_by_layer = [cells_in_L23; cells_in_L4; cells_in_L5; cells_in_L6];
    
end



%
%  CELL DEPTH FROM PIA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(centers,1)==0
    cellStats.cell_depth = [];
    
else
    
    
    
    % determine the Y-value for the L1-L2/3 boundary at the X-position
    % corresponding to the cell centers.
    cent_x = centers(:,1);
    cent_y = centers(:,2);
    y_on_boundary1 = (cellFillData.raw.boundary1.m).*cent_x + (cellFillData.raw.boundary1.b);
    cellStats.cell_depth = abs(y_on_boundary1 - cent_y);
    
end



%
%  CELL SIZE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(centers,1)==0
    cellStats.cell_size = [];
    
else
    
    % region props want to work on a 2D mask, so project the 3D mask onto 2
    % dimensions.
    tmp_size = regionprops(mask, 'area');
    cellStats.cell_size = cat(1, tmp_size(:).Area);
    
end







