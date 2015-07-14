function cellStats = getCellFillStats(cellFillData, params)

% import the raw data
mask = cellFillData.mask.img;

%
% calculate the centers of each region in the mask. Occasionally a single
% region is split into non-contiguous volumes. this is rare, and I should
% fix these cases.
%
% also calculate the volume, and diameter of each region
%
%%%%%%%%%%%%%%%%%%%%%%%%
uniqueRegions = unique(mask(:));
uniqueRegions = uniqueRegions(uniqueRegions>0);
Nregions = numel(uniqueRegions);
centers = nan(Nregions, 3);
cellStats.cell_volume = nan(Nregions, 1);
cellStats.cell_diam = nan(Nregions, 1);
cellStats.cell_on_edge = false(Nregions, 1);
cellStats.minor_axis = nan(Nregions, 1);
cellStats.major_axis = nan(Nregions, 1);
cellStats.orientation = nan(Nregions, 1);
for i_reg = 1:Nregions
    
    tmp_mask = mask == uniqueRegions(i_reg);
    rprops_3D = regionprops(tmp_mask, 'Area', 'Centroid', 'PixelIdxList');
    
    if ~all(cat(1, rprops_3D(:).Area) == 0) % ignore regions that have no pixels
        
        % deal with cases where there are multiple regions with the same
        % label number
        if numel(rprops_3D) > 1
            tmp_vols = cat(1, rprops_3D(:).Area);
            [~, idx] = max(tmp_vols);
            
            % modify the tmp_mask so that all the non-cellfill regions are
            % gone
            tmp_mask = false(size(tmp_mask));
            tmp_mask(rprops_3D(idx).PixelIdxList) = true;
            rprops_3D = rprops_3D(idx);
            
            fprintf('found %d subregions:', numel(tmp_vols))
            eval('transpose(tmp_vols)')
        end
        
        % calculated the 2D stats
        rprops_2D = regionprops(max(tmp_mask, [], 3), 'EquivDiameter', 'MajorAxisLength', 'MinorAxisLength', 'Orientation');
        
                
        % now assign the centers, volume, diam
        centers(i_reg,:) = rprops_3D.Centroid;
        cellStats.cell_volume(i_reg) = rprops_3D.Area;
        cellStats.cell_diam(i_reg) = rprops_2D.EquivDiameter;
        cellStats.minor_axis(i_reg) = rprops_2D.MinorAxisLength;
        cellStats.major_axis(i_reg) = rprops_2D.MajorAxisLength;
        cellStats.orientation(i_reg) = rprops_2D.Orientation;
        
        
        % flag cases where the region touches any of the edges of the
        % counting window
        critval = 1;
        [row, col, ~] = ind2sub(size(tmp_mask), find(tmp_mask==1));
        on_L_edge = sum(col == 1) >= critval;
        on_R_edge = sum(col == size(tmp_mask, 2)) >= critval;
        on_T_edge = sum(row == 1) >= critval;
        on_B_edge = sum(row == size(tmp_mask, 1)) >= critval;
        cellStats.cell_on_edge(i_reg) = any([on_L_edge, on_R_edge, on_T_edge, on_B_edge]);        
        
    else
        disp('found region with no volume')
    end
end



%
% LAMINAR VOLUME
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
    cellStats.layerAssignments = nan(Nregions,1);
else
    
    cellStats.layerAssignments = nan(size(centers,1),1);
    
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
    
    
    % Find cells in L1
    cells_above_L1 = centers(:,2) < y_on_boundary1;
    cellStats.layerAssignments(cells_above_L1) = 1;
    
    % Find cells in L2/3
    cells_below_L1 = centers(:,2) >= y_on_boundary1;
    cells_above_L4 = centers(:,2) < y_on_boundary2;
    cells_in_L23 = sum(cells_above_L4 & cells_below_L1);
    cellStats.layerAssignments(cells_above_L4 & cells_below_L1) = 23;
    
    % Find cells in L4
    cells_below_L3 = centers(:,2) >= y_on_boundary2;
    cells_above_L5 = centers(:,2) < y_on_boundary3;
    cells_in_L4 = sum(cells_below_L3 & cells_above_L5);
    cellStats.layerAssignments(cells_below_L3 & cells_above_L5) = 4;
    
    % Find cells in L5
    cells_below_L4 = centers(:,2) >= y_on_boundary3;
    cells_above_L6 = centers(:,2) < y_on_boundary4;
    cells_in_L5 = sum(cells_above_L6 & cells_below_L4);
    cellStats.layerAssignments(cells_above_L6 & cells_below_L4) = 5;
    
    % Find cells in L6
    cells_below_L5 = centers(:,2) >= y_on_boundary4;
    cells_above_endofcortex = centers(:,2) < y_on_boundary5;
    cells_in_L6 = sum(cells_above_endofcortex & cells_below_L5);
    cellStats.layerAssignments(cells_above_endofcortex & cells_below_L5) = 6;
    
    % Find cells that are below the cortex
    cells_below_cortex = centers(:,2) > y_on_boundary5;
    cellStats.layerAssignments(cells_below_cortex) = inf;
    
    
    % assign the outputs
    cellStats.count_by_layer = [cells_in_L23; cells_in_L4; cells_in_L5; cells_in_L6];
    
end



%
%  CELL DEPTH FROM PIA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(centers,1)==0
    cellStats.cell_depth = nan(Nregions,1);
    
else
    
    
    
    % determine the Y-value for the L1-L2/3 boundary at the X-position
    % corresponding to the cell centers.
    cent_x = centers(:,1);
    cent_y = centers(:,2);
    y_on_boundary1 = (cellFillData.raw.boundary1.m).*cent_x + (cellFillData.raw.boundary1.b);
    cellStats.cell_depth = abs(y_on_boundary1 - cent_y);
    
end


% some error checking
assert(~any(isnan([cellStats.volume_by_layer; cellStats.count_by_layer])), 'ERROR: NaN output of cellcount analysis')
assert(~any(isnan(cellStats.cell_depth)), 'ERROR: NaN output for cell depth');
assert(~any(isnan(cellStats.cell_volume)), 'ERROR: NaN output for cell volume');
assert(~any(isnan(cellStats.cell_diam)), 'ERROR: NaN output for cell diam');
assert(~any(isnan(cellStats.layerAssignments)), 'ERROR: NaN output for layer assignment');
N = [size(cellStats.cell_depth);...
    size(cellStats.cell_volume);...
    size(cellStats.cell_diam);...
    size(cellStats.layerAssignments)];
assert(size(unique(N, 'rows'), 1) == 1, 'ERROR: inconsistent numbers of outputs');



