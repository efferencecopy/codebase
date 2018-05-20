%% import the data

fin

% read in the spreadsheet that defines a bunch of things for the population
% analysis
%
% sheet 1 = PVcre cells
% sheet 2 = GIN cells
% sheet 3 = SOMcre cells
%

CELLTYPE = 'PVcre';

xlspath = [GL_DOCUPATH 'Other_workbooks', filesep, 'Interneuron_density_analysis.xlsx'];
[~, txt, raw] =xlsread(xlspath, CELLTYPE);
header = txt(1,:);
hidx = struct();
for i_head = 1:numel(header)
    if numel(header{i_head}) > 0
        head_name = regexprep(header{i_head}, ' ', '_');
        head_name = regexprep(head_name, '/', '');
        hidx.(head_name) = i_head;
    end
end
mouseName = txt(2:end, hidx.Mouse);
fileName = txt(2:end, hidx.Filename);
brainArea_fpath = txt(2:end, hidx.HVA_directory);
brainArea_analysis = txt(2:end, hidx.HVA_analysis);
Nfiles = size(fileName, 1);
ap_dist = cat(1, raw{2:end, hidx.ap_distance_to_cc});
ap_dist(Nfiles+1:end) = []; % cut the junk that comes along for the ride
ml_dist = cat(1, raw{2:end, hidx.ml_distance});
ml_dist(Nfiles+1:end) = [];

startingPath = [GL_DATPATH(1:end-5), 'SOM_PV_Density', filesep];

% initalize the output variables
results = [];

for i_file = 1:Nfiles
    
    fprintf('Analyzing file %d of %d \n', i_file, Nfiles)
    
    % locate the .mat file for the cell fill mask that contains the laminar
    % boundaries already marked.
    tmp_fname = fileName{i_file}; %specify a file
    tmp_area = brainArea_fpath{i_file};
    tmp_mouse = mouseName{i_file};
    optpath = [startingPath, CELLTYPE, filesep, tmp_mouse, filesep, 'confocal', filesep, tmp_area];
    
    out = findfile(tmp_fname, optpath, '.mat'); % recursively look for it
    load(out); % load the data
    
    % define some params by getting the info for this stack
    plate_slice = regexpi(tmp_fname, 'p[\d]s[\d]', 'match');
    d = dir(optpath);
    istif = cellfun(@(x) ~isempty(x), regexp({d(:).name}, 'tif'));
    ismatch = cellfun(@(x) ~isempty(x), regexp({d(:).name}, plate_slice));
    assert(sum(istif & ismatch)==1, 'ERROR: did not find match')
    info_idx = find(istif & ismatch);
    out = findfile(d(info_idx).name, optpath, '.tif');
    info = imfinfo(out);
    resolution = info(1).XResolution;
    params.pix_per_um = resolution; % this comes from the LSM files, but may be inconsistent from file to file...
    thickness = getSliceThickness(out, info);
    if strcmpi(tmp_mouse, 'AK090314A') && strcmpi(plate_slice, 'p6s1')
        params.slice_thickness_um = 35; % this is just hard coded
    elseif strcmpi(tmp_mouse, 'CH_150612_A') && strcmpi(plate_slice, 'p5s6')
        params.slice_thickness_um = 35; % this is just hard coded
    else
        params.slice_thickness_um = 70; % this is just hard coded
    end
    
    % get the cell stats
    cellStats = getCellFillStats(cellFillData, params); % run analysis code
    
    % Add output variables to arrays
    results(i_file).cellsByLayer = cellStats.count_by_layer;
    results(i_file).volumeByLayer = cellStats.volume_by_layer;
    results(i_file).cellDepths = cellStats.cell_depth;
    results(i_file).cellVolume = cellStats.cell_volume;
    results(i_file).cellDiam = cellStats.cell_diam;
    results(i_file).minorAxis = cellStats.minor_axis;
    results(i_file).majorAxis = cellStats.major_axis;
    results(i_file).orientation = cellStats.orientation;
    results(i_file).layerAssignments = cellStats.layerAssignments;
    results(i_file).cell_on_edge= cellStats.cell_on_edge;
    results(i_file).thickness = thickness;

end

fprintf('  ** Done loading data **\n')

%% PLOTTING ROUTINES:  CELL DENSITY

close all; clc

Nmice = size(unique(mouseName),1);
Nareas = size(unique(brainArea_analysis), 1);
mice = unique(mouseName);
areas = unique(brainArea_analysis);

% Pre-assemble the population structure. Assume that all mice were tested in
% every brain area, and that there are 4 layers per brain area. This is not
% true in every case, so pad with NaNs.

popdat = [];
allHVAs = {'PM', 'AL', 'AM', 'RL', 'ERC'};
for i_area = 1:numel(allHVAs)
    popdat.(allHVAs{i_area}).totalVolume = nan(4, Nmice); %The array is Nlayers x Nmice
    popdat.(allHVAs{i_area}).cellCount = nan(4, Nmice);
end

for i_mouse = 1:Nmice
   
    for i_area = 1:Nareas
        
        % grab the raw data;
        idx = strcmpi(mouseName, mice{i_mouse}) & strcmpi(brainArea_analysis, areas{i_area});
        if sum(idx)==0
         
            continue % this brain area was not tested in this mouse.
            
        elseif sum(idx) >= 4
            
            % A hack to only include 4 of the slices from EB_150427_A
            if strcmpi(mice{i_mouse}, 'EB_150427_A')
                tmp_inds = find(idx);
                if strcmpi(areas(i_area), 'pm')
                    idx(1:tmp_inds(end)-4) = false;
                elseif strcmpi(areas(i_area), 'al')
                    idx(tmp_inds(1):tmp_inds(1)+3) = false;
                    idx(tmp_inds(1)+4+4:tmp_inds(end)) = false;
                end
            end
            
%             % A hack to only include 4 of the slices from CH_150612_A
%             if strcmpi(mice{i_mouse}, 'CH_150612_A')
%                 tmp_inds = find(idx);
%                 if strcmpi(areas(i_area), 'pm')
%                     idx(tmp_inds(2)) = false;
%                     idx(tmp_inds(7):tmp_inds(end)) = false;
%                 end
%             end
                    
            
            
            % make sure that there are no duplicate files
            assert(size(unique(fileName(idx)),1)== sum(idx), 'ERROR: duplicate files found')
            
            % combine data across the 4 brain slices
            tmp_volume = cat(2, results(idx).volumeByLayer);
            tmp_counts = cat(2, results(idx).cellsByLayer);
            
            % sum across the 4 slices to arrive at a representative value for
            % each mouse
            tmp_volume = sum(tmp_volume, 2);
            tmp_counts = sum(tmp_counts, 2);
            
        else
            error('Incorrect number of matches')
        end
        
        % the population data. The array is Nlayers x Nmice
        hva = areas{i_area};
        popdat.(hva).totalVolume(:, i_mouse) = tmp_volume;
        popdat.(hva).cellCount(:, i_mouse) = tmp_counts;
        
        
    end    
    
end



% simple plot of volume, counts, density for each area and layer
areas = {'PM', 'AL', 'AM', 'RL', 'ERC'};
figure
set(gcf, 'position', [184    35   764   746])
for i_area = 1:numel(areas);
    
    tmp_volume = popdat.(areas{i_area}).totalVolume;
    tmp_counts = popdat.(areas{i_area}).cellCount;
    
    % density = {nAreas}[nLayers, nMice]
    volume{i_area} = tmp_volume;
    counts{i_area} = tmp_counts;
    density{i_area} = tmp_counts ./ tmp_volume;
    
    subplot(6, 3, (i_area-1)*3 + 1)
    plot(tmp_counts, '.-')
    xlabel('layer')
    ylabel('counts')
    title(areas{i_area})
    box off
    axis tight
    
    subplot(6, 3, (i_area-1)*3 + 2)
    plot(tmp_volume, '.-')
    xlabel('layer')
    ylabel('volume')
    title(areas{i_area})
    box off
    axis tight
    
    subplot(6, 3, (i_area-1)*3 + 3)
    plot(density{i_area}, '.-')
    xlabel('layer')
    ylabel('density')
    title(areas{i_area})
    box off
    axis tight
    
    
end


% Plot of density for each layer, comparing across areas
figure
set(gcf, 'position', [223    31   573   754]);
layerLabel = {'2/3', '4', '5', '6'};
for i_layer = 1:4
    
    tmp_density = [];
    tmp_volume = [];
    tmp_counts = [];
    for i_area = 1:numel(areas)
            tmp_density(i_area, :) = density{i_area}(i_layer,:); % [Nareas, Nmice]
            tmp_volume(i_area,:) = volume{i_area}(i_layer,:); % [Nareas, Nmice]
            tmp_counts(i_area,:) = counts{i_area}(i_layer,:); % [Nareas, Nmice]
    end
    
    subplot(4,3, (i_layer-1)*3 +1)
    hold on,
    plot(tmp_counts, '.--', 'linewidth', 0.25)
    plot(nanmean(tmp_counts, 2), '.k-', 'linewidth', 3)
    box off; axis tight
    ylabel('counts')
    xlabel('Brain Area')
    title(sprintf('Layer %s', layerLabel{i_layer}));
    
    subplot(4,3, (i_layer-1)*3 +2)
    hold on,
    plot(tmp_volume, '.--', 'linewidth', 0.25)
    plot(nanmean(tmp_volume, 2), '.k-', 'linewidth', 3)
    box off; axis tight
    ylabel('volume')
    xlabel('Brain Area')
    title(sprintf('Layer %s', layerLabel{i_layer}));
    
    subplot(4,3, (i_layer-1)*3 +3)
    hold on,
    plot(tmp_density, '.--', 'linewidth', 0.25)
    plot(nanmean(tmp_density, 2), '.k-', 'linewidth', 3)
    box off; axis tight
    ylabel('density')
    xlabel('Brain Area')
    title(sprintf('Layer %s', layerLabel{i_layer}));
    
end


% Plot of density, integrated across specific layers, and comparing across
% brain areas
layers = [1]; % 1= L2/3, 2=L4, 3=L5, 4=L6
normArea = 2; % 1=PM, 2=AL
densityAcrossLayers = [];
for i_area = 1:numel(areas);
    
    tmp_volume = volume{i_area}(layers,:); % [Nareas, Nmice]
    tmp_counts = counts{i_area}(layers,:); % [Nareas, Nmice]
    
    tmp_density = sum(tmp_counts, 1) ./ sum(tmp_volume, 1);
    densityAcrossLayers(i_area,:) = tmp_density; % notice that the dim is [Nareas x Nmice]
end

figure, hold on,
plot(densityAcrossLayers, 'o--')
nareas = numel(areas);
xbar = nanmean(densityAcrossLayers,2);
sem = nanstd(densityAcrossLayers,[], 2) ./ sqrt(sum(~isnan(densityAcrossLayers), 2));
errorbar(1:nareas, xbar, sem, 'k', 'linewidth', 3)
set(gca, 'xtick', 1:nareas, 'xTickLabel', areas)
legend(unique(mouseName))

% now plotting percent change relative to a single area
figure, hold on,
diff_vals = bsxfun(@minus, densityAcrossLayers, densityAcrossLayers(normArea,:));
prcnt_change = bsxfun(@rdivide, diff_vals, densityAcrossLayers(normArea,:));
plot(prcnt_change, 'o--')
nareas = numel(areas);
xbar = nanmean(prcnt_change, 2);
sem = nanstd(prcnt_change,[], 2) ./ sqrt(sum(~isnan(prcnt_change), 2));
errorbar(1:nareas, xbar, sem, 'k', 'linewidth', 3)
set(gca, 'xtick', 1:nareas, 'xTickLabel', areas)
ylabel(sprintf('percent change (relative to %s)', areas{normArea}))
legend(unique(mouseName))

%% PLOTTING ROUTINES: CELL SIZE



% concatenate the data across mice and slices
areas = {'PM','AL', 'RL', 'ERC'};
for i_area = 1:numel(areas)
    
    % grab the raw data;
    idx = strcmpi(brainArea_analysis, areas{i_area});
    popdat.(areas{i_area}).cellDepth = cat(1, results(idx).cellDepths);
    popdat.(areas{i_area}).cellDiam = cat(1, results(idx).cellDiam);
    popdat.(areas{i_area}).minorAxis = cat(1, results(idx).minorAxis);
    popdat.(areas{i_area}).majorAxis = cat(1, results(idx).majorAxis);
    popdat.(areas{i_area}).ori = cat(1, results(idx).orientation);
    popdat.(areas{i_area}).layerAssignments = cat(1, results(idx).layerAssignments);
    
    % which cells were touching an edge
    on_edge = cat(1, results(idx).cell_on_edge);
    
    % remove the cells that were on an edge:
    popdat.(areas{i_area}).cellDepth(on_edge) = [];
    popdat.(areas{i_area}).cellDiam(on_edge) = [];
    popdat.(areas{i_area}).layerAssignments(on_edge) = [];
    popdat.(areas{i_area}).minorAxis(on_edge) = [];
    popdat.(areas{i_area}).majorAxis(on_edge) = [];
    popdat.(areas{i_area}).ori(on_edge) = [];
end

% plot a histogram of cell sizes across all areas
bigdataset_diam = [];
bigdataset_ori = [];
bigdataset_MM = [];
bigdataset_layers = [];
for i_area = 1:numel(areas)
    bigdataset_diam = cat(1, bigdataset_diam, popdat.(areas{i_area}).cellDiam);
    bigdataset_ori = cat(1, bigdataset_ori, popdat.(areas{i_area}).ori);
    bigdataset_layers = cat(1, bigdataset_layers, popdat.(areas{i_area}).layerAssignments);
    
    tmp = popdat.(areas{i_area}).majorAxis ./ popdat.(areas{i_area}).minorAxis;
    bigdataset_MM = cat(1, bigdataset_MM, tmp);
end

figure
set(gcf, 'position', [53    13   799   769])
layerTypes = [23, 4, 5, 6];
for i_lyr = 1:4
    
    idx = bigdataset_layers == layerTypes(i_lyr);
    
    % cell sizes
    subplot(4,3, (i_lyr-1)*3 + 1)
    bins = linspace(min(bigdataset_diam), max(bigdataset_diam), 80);
    tmp = bigdataset_diam(idx);
    hist(tmp, bins); hold on,
    ylims = get(gca, 'ylim');
    plot([mean(tmp), mean(tmp)], ylims, 'b')
    xlim([bins(1), bins(end)])
    title(sprintf('Layer %d', layerTypes(i_lyr)))
    ylabel('counts')
    xlabel('cell diameter')
    
    % polar plot of orientation
    subplot(4,3, (i_lyr-1)*3 + 2)
    ori_rad = (bigdataset_ori(idx) ./ 180) .* pi;
    rose(ori_rad, 40)
    title('ori from horiz')
    axis equal
    
    
    % major/minor ratio
    subplot(4,3, (i_lyr-1)*3 + 3)
    bins = logspace(log10(min(bigdataset_MM)), log10(max(bigdataset_MM)), 100);
    N = histc(bigdataset_MM(idx), bins);
    bar(bins, N);
    gm = geomean(bigdataset_MM(idx))
    set(gca, 'xscale', 'log');
    set(gca, 'xtick', [1 2 4], 'xticklabel', [1 2 4])
%     child = get(gca, 'children');
%     set(child(1), 'visible', 'off');
%     set(child(2), 'edgealpha', 0.1);
    hold on,
    plot([gm, gm], [0, max(N)], 'b')
    xlim([1, bins(end)])
    title(sprintf('Layer %d', layerTypes(i_lyr)))
    ylabel('counts')
    xlabel('Major : Minor ratio')
    
    
end

%% WITHIN MOUSE COMPARISON OF CELL SIZES ACROSS LAMINAE AND AREAS



Nmice = size(unique(mouseName),1);
Nareas = size(unique(brainArea_analysis), 1);
mice = unique(mouseName);
areas = unique(brainArea_analysis);

% Pre-assemble the population structure. Assume that all mice were tested in
% every brain area, and that there are 4 layers per brain area. This is not
% true in every case, so pad with NaNs.

allHVAs = {'PM','AL', 'ERC', 'RL'};
for i_area = 1:numel(allHVAs)
    % means and SEM
    popdat.(allHVAs{i_area}).diam_xbar = nan(4, Nmice); %The array is Nlayers x Nmice
    popdat.(allHVAs{i_area}).diam_sem = nan(4, Nmice); %The array is Nlayers x Nmice
end

for i_mouse = 1:Nmice
   
    for i_area = 1:Nareas
        
        % grab the raw data;
        idx = strcmpi(mouseName, mice{i_mouse}) & strcmpi(brainArea_analysis, areas{i_area});
        if sum(idx)==0
            
            continue % this brain area was not tested in this mouse.
            
        end
        
        % make sure that there are no duplicate files
        assert(size(unique(fileName(idx)),1)== sum(idx), 'ERROR: duplicate files found')
        
        % combine data across the 4 brain slices
        tmp_diam = cat(1, results(idx).cellDiam);
        tmp_layers = cat(1, results(idx).layerAssignments);
        tmp_on_edge = cat(1, results(idx).cell_on_edge);
        
        % any nan values?
        assert(~any(isnan([tmp_diam; tmp_layers])),...
            'ERROR: found a NaN value')
        
        layerTypes = [23, 4, 5, 6];
        for i_lyr = 1:4;
            
            idx = tmp_layers == layerTypes(i_lyr);
            
            % remove the cells that were on the edges
            idx = idx & ~tmp_on_edge;
            
            % compile the raw data
            xbar = mean(tmp_diam(idx));
            sem = stderr(tmp_diam(idx));
            
            % compile the mean and SEM
            popdat.(areas{i_area}).diam_xbar(i_lyr, i_mouse) = xbar;
            popdat.(areas{i_area}).diam_sem(i_lyr, i_mouse) = sem;
            
        end
        
    end    
    
end

% plot xbar across areas for each layer

figure
set(gcf, 'position', [235    50   244   723])
allHVAs = {'PM','AL', 'ERC', 'RL'};
for i_lyr = 1:4
    
    % make some tmp data [Nareas x Nmice] where each element is the xbar
    % for layer = i_lyr
    tmp_dat = nan(Nareas, Nmice);
    for i_area = 1:numel(allHVAs)
        tmp_dat(i_area, :) = popdat.(allHVAs{i_area}).diam_xbar(i_lyr, :);
    end
    
    subplot(4,1,i_lyr)
    plot([1:numel(allHVAs)], tmp_dat, '-')
    title(sprintf('Layer: %d', layerTypes(i_lyr)))
    xlabel('Brain Area')
    ylabel('cell diameter')
    set(gca, 'xtick', 1:numel(allHVAs), 'xticklabel', allHVAs)
end




%% PLOTTING ROUTINES: CELL DEPTHS


% plot a histograms of cell depths
figure
set(gcf, 'position', [108    88   863   663])
for i_area = 1:numel(areas)
    
    subplot(numel(areas),1,i_area)
    tmp = popdat.(areas{i_area}).cellDepth;
    bins = linspace(0, 1000, 100);
    hist(tmp, bins)
    set(get(gca, 'children'), 'edgealpha', 0.1)
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor',[0 .5 .5])
    ylims = get(gca, 'ylim');
    hold on,
    xbar = nanmean(tmp);
    plot([xbar, xbar], ylims, 'b', 'linewidth', 3)
    xlim([0, 1000])
    title(areas{i_area})
    ylabel('counts')
    xlabel('cell depth from L1')
    box off
end



% scatter plot of cell size vs. cell depth
figure
set(gcf, 'position', [109   384   887   368])
bigdata_depth = [];
bigdata_diam = [];
for i_area = 1:numel(areas)
    
    subplot(1,numel(areas),i_area)
    tmp_depth = popdat.(areas{i_area}).cellDepth;
    tmp_size = popdat.(areas{i_area}).cellDiam;
    plot(tmp_depth, tmp_size, '.')
    ylim([0,25])
    xlim([0,1000])
    xlabel('cell depth')
    ylabel('cell diameter (um)')
    [r, p] = corr(tmp_depth, tmp_size, 'type', 'spearman');
    title(sprintf('%s = %.2f', areas{i_area}, p))
    
    bigdata_depth = cat(1, bigdata_depth, tmp_depth);
    bigdata_diam = cat(1, bigdata_diam, tmp_size);
end

figure
plot(bigdata_depth, bigdata_diam, '.')
ylim([0,25])
xlim([0,1000])
xlabel('cell depth')
ylabel('cell diameter (um)')
[r, p] = corr(bigdata_depth, bigdata_diam, 'type', 'spearman');
title(sprintf('All Areas Combined = %.2f', p))



%% TOPOGRAPHY OF CELL DENSITY 

close all

% integrate the cell counts and volume across all layers. The output is a
% single value for each data file
layers = [1:3];
tmp_counts = cat(2, results(:).cellsByLayer);
tmp_counts = sum(tmp_counts(layers, :), 1);
tmp_volume = cat(2, results(:).volumeByLayer);
tmp_volume = sum(tmp_volume(layers, :), 1);
tmp_density = tmp_counts ./ tmp_volume;


% iterate over mice, and brain areas. color code each counting site with
% clr_raw, and then color code the mean location with a square (with
% clr_avg).
Nmice = size(unique(mouseName),1);
Nareas = size(unique(brainArea_analysis), 1);
mice = unique(mouseName);
areas = unique(brainArea_analysis);
PLOTRAW = false;

figure; hold on,
set(gcf, 'position', [296   142   575   620])
for i_mouse = 1:Nmice
    for i_area = 1:Nareas
        
        idx = strcmpi(mice{i_mouse}, mouseName) & strcmpi(areas{i_area}, brainArea_analysis);
        
        X = -ml_dist(idx);
        Y = -ap_dist(idx);
        
        tmp_density_cat = sum(tmp_counts(idx)) ./ sum(tmp_volume(idx)); % represents all slices for a brain area
        
        [clr_raw, clr_avg] = hvaPlotColor(areas{i_area});
        p = plot3(mean(X), mean(Y), tmp_density_cat, 's', 'color', clr_avg, 'markerfacecolor', clr_avg);
        bdfxn = @(x,y,z) title(z);
        set(p, 'buttonDownFcn', {bdfxn, mice{i_mouse}})
        
        
        if PLOTRAW
            plot3(X, Y, tmp_density(idx), 'o', 'color', clr_raw)
        end
    end
end
xlim([-1800 310])
ylim([-1800 310])
xlabel('M/L distance from Callosum')
ylabel('A/P distance from Callosum')
zlabel('Density')
%set(gca, 'view', [-48, 32])


%% TOPOGRAPHY OF COUNTING REGIONS ACROSS ALL MICE


PLOTRAW = false;
NORMDIST = true;
PLOTINJ = false;

% GRAB ALL THE DATA FROM ALL MICE
[all_ap_dist, all_ml_dist, all_norm_ap, all_norm_ml, all_inj_ap] = deal([]);
[all_mouseName, all_brainArea_analysis] = deal({});
for i_type = 1:3
    xlspath = [GL_DOCUPATH 'Other_workbooks', filesep, 'Interneuron_density_analysis.xlsx'];
    [~, txt, raw] =xlsread(xlspath, i_type);
    error('need to generate a new header index dict')
    
    all_mouseName = cat(1, all_mouseName, txt(2:end, hidx.Mouse));
    all_brainArea_analysis = cat(1, all_brainArea_analysis, txt(2:end, hidx.HVA_analysis));
    Nfiles = size(txt, 1)-1;
    
    tmp_ap = cat(1, raw{2:end, hidx.ap_distance_to_cc});
    tmp_ap(Nfiles+1:end) = []; % cut the junk that comes along for the ride
    all_ap_dist = cat(1, all_ap_dist, tmp_ap);
    
    tmp_ml = cat(1, raw{2:end, hidx.ml_distance_to_cc});
    tmp_ml(Nfiles+1:end) = [];
    all_ml_dist = cat(1, all_ml_dist, tmp_ml);
    
    tmp_norm_fact = cat(1, raw{2:end, hidx.AP_norm_fact});
    tmp_norm_fact(Nfiles+1:end) = []; 
    all_norm_ap = cat(1, all_norm_ap, tmp_norm_fact);
    
    tmp_norm_fact = cat(1, raw{2:end, hidx.ML_norm_fact});
    tmp_norm_fact(Nfiles+1:end) = [];
    all_norm_ml = cat(1, all_norm_ml, tmp_norm_fact);
    
    tmp_inj_ap = cat(1, raw{2:end, hidx.V1_inj_um_from_cc});
    tmp_inj_ap(Nfiles+1:end) = [];
    all_inj_ap = cat(1, all_inj_ap, tmp_inj_ap);
end


% convert the ML distances and scale factors to um (from pix). For the
% retiga camera on the nikon scope 2x, 1000 um = 235 pix
umPerPix = 1000/235;
all_ml_dist = all_ml_dist .* umPerPix;
all_norm_ml = all_norm_ml .* umPerPix;


% iterate over mice, and brain areas. color code each counting site with
% clr_raw, and then color code the mean location with a square (with
% clr_avg).
Nmice = size(unique(all_mouseName),1);
Nareas = size(unique(all_brainArea_analysis), 1);
mice = unique(all_mouseName);
areas = unique(all_brainArea_analysis);

figure; hold on,
set(gcf, 'position', [296   142   575   620])
for i_mouse = 1:Nmice
    for i_area = 1:Nareas
        
        idx = strcmpi(mice{i_mouse}, all_mouseName) & strcmpi(areas{i_area}, all_brainArea_analysis);
        
        if sum(idx) == 0
            continue
        end
        
        X = -all_ml_dist(idx);
        Y = -all_ap_dist(idx);
        inj_y = -unique(all_inj_ap(idx));
        inj_x = -2500; % in um.
        
        if NORMDIST
            norm_ml = unique(all_norm_ml(idx));
            norm_ap = unique(all_norm_ap(idx));
            
            % including this term effecitvely means the norm_ap does
            % nothing, and I'm scaling both dimensions by the norm_ml
            % factor
            norm_fact_equator = norm_ap ./ norm_ml;
            
            X = X./norm_ml;
            Y = Y./norm_ap*norm_fact_equator;
            inj_x = inj_x./norm_ml;
            inj_y = inj_y./norm_ap*norm_fact_equator;
        end
        
        if PLOTINJ
            if any(strcmpi(areas{i_area}, {'AL', 'PM'}))
                try
                    plot([inj_x, mean(X)], [inj_y, mean(Y)], '-', 'color', [0.3 0.3 0.3])
                catch
                    disp('need to add inj site AP data')
                end
            end
        end
        
                
        if PLOTRAW
            plot(X, Y, 'o', 'color', clr_raw)
        end
        
        % do this last so that the average point are on top of everything
        % else
        [clr_raw, clr_avg] = hvaPlotColor(areas{i_area});
        p = plot(mean(X), mean(Y), 's', 'color', clr_avg, 'markerfacecolor', clr_avg);
        bdfxn = @(x,y,z) title(z);
        set(p, 'buttonDownFcn', {bdfxn, mice{i_mouse}})
        
        

        
        
    end
end
axis equal
xlabel('M/L distance from Callosum')
ylabel('A/P distance from Callosum')









        