%% make a histogram of cell depths for SOM+ and PV+ cells



fin

% read in the spreadsheet that defines a bunch of things for the population
% analysis

cd([GL_DOCUPATH, filesep, 'Other_Workbooks'])
[~,txt] =xlsread('Amy_counting_data.xlsx',1);
cellType = txt(2:end, 1);
fileName = txt(2:end, 2);
mouseName = txt(2:end, 3);
brainArea = txt(2:end,4);
startingPath = '\\crash.dhe.duke.edu\data\home\amy\counting\';


% initalize the output variables
depth = {};
xpos = {};
for i_file = 1:size(fileName,1)
    
    fprintf('File <%d> of <%d>\n', i_file, size(fileName,1));
    
    
    % locate the .mat file for the cell fill mask that contains the laminar
    % boundaries already marked.
    tmp_fname = fileName{i_file}; %specify a file
    tmp_area = brainArea{i_file};
    tmp_mouse = mouseName{i_file};
    optpath = [startingPath, tmp_mouse, filesep, 'confocal', filesep, tmp_area];
    
    out = findfile(tmp_fname, optpath, '.mat'); % recursively look for it
    load(out); % load the data
    mask = cellFillData.mask.img;
    
    % find centers
    rprops = regionprops(mask, 'centroid');
    centers = cat(1, rprops.Centroid);
    if isempty(centers)
        fprintf('file %s from area %s has no cells\n', tmp_fname, tmp_area);
        continue
    end
    
    % define the upper most boundary that devides layer 1 from layer 2/3
    if ~isfield(cellFillData.raw, 'boundary1')
        error('Could not find Layer 1 boundary')
    end
    
    % now figure out the vertical distance to the L1 boundary
    m = cellFillData.raw.boundary1.m;
    b = cellFillData.raw.boundary1.b;
    predY = m .* centers(:,1) + b;
    dy = centers(:,2) - predY;
    
    % assign to a cell array for now
    depth{i_file} = dy;
    xpos{i_file} = centers(:,1);
    
end


intype = {'SOM', 'PV'};
for i_type = 1:numel(intype)

   figure
   % just plot a histogram of all the data (across areas)
   maxDepth = max(cat(1, depth{:}));
   binWidth = 25;
   nBins = ceil(maxDepth./binWidth);
   bins = 0 : binWidth : binWidth*nBins;
   
   l_type =  cellfun(@(x,y) regexpi(x,y), cellType, repmat(intype(i_type), size(cellType)), 'uniformoutput', false);
   l_type = cellfun(@(x) ~isempty(x), l_type);
   
   tmp_depths = cat(1, depth{l_type});
   counts = histc(tmp_depths, bins);
   
   bar(bins, counts, 'histc')
   
   
    
end






