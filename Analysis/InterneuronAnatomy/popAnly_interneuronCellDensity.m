%% import the data

fin

% read in the spreadsheet that defines a bunch of things for the population
% analysis
%
% sheet 1 = PV-Cre cells
% sheet 2 = GIN cells
% sheet 3 = SOM-Cre cells
%

CELLTYPE = 2;

cd ~/Crash/Data/SOM_PV_Density
[~,txt] =xlsread('counting_data_sheet.xlsx', CELLTYPE);
mouseName = txt(2:end, 1);
fileName = txt(2:end, 2);
brainArea = txt(2:end,3);
startingPath = '~/Crash/Data/SOM_PV_Density/';


% initalize the output variables
results = [];
for i_file = 1:size(fileName,1)
    
    % locate the .mat file for the cell fill mask that contains the laminar
    % boundaries already marked.
    tmp_fname = fileName{i_file}; %specify a file
    tmp_area = brainArea{i_file};
    tmp_mouse = mouseName{i_file};
    optpath = [startingPath, tmp_mouse, filesep, 'confocal', filesep, tmp_area];
    
    out = findfile(tmp_fname, optpath, '.mat'); % recursively look for it
    load(out); % load the data
    
    % define some params by getting the info for this stack
    udrscr = regexp(tmp_fname,'_');
    info = imfinfo(findfile(tmp_fname(udrscr(2)+1:end), optpath, '.tif'));
    resolution = info(1).XResolution;
    params.pix_per_um = resolution; % this comes from the LSM files, but may be inconsistent from file to file...
    params.slice_thickness_um = 70; % this is just hard coded
    
    cellStats = getCellFillStats(cellFillData, params); % run analysis code
    assert(~any(isnan([cellStats.volume_by_layer; cellStats.count_by_layer])), 'ERROR: NaN output of cellcount analysis')
    
    % Add output variables to arrays
    results(i_file).cellsByLayer = cellStats.count_by_layer;
    results(i_file).volumeByLayer = cellStats.volume_by_layer;
    results(i_file).cellDepths = cellStats.cell_depth;
    results(i_file).cellSize = cellStats.cell_size;
    
end

%% PLOTTING ROUTINES:  CELL DENSITY



Nmice = size(unique(mouseName),1);
Nareas = size(unique(brainArea), 1);
mice = unique(mouseName);
areas = unique(brainArea);

% Pre-assemble the population structure. Assume that all mice were tested in
% every brain area, and that there are 4 layers per brain area. This is not
% true in every case, so pad with NaNs.
popdat = [];
allHVAs = {'PM','AL', 'ERC', 'RL'};
for i_area = 1:numel(allHVAs)
    popdat.(allHVAs{i_area}).totalVolume = nan(4, Nmice); %The array is Nlayers x Nmice
    popdat.(allHVAs{i_area}).cellCount = nan(4, Nmice);
end

for i_mouse = 1:Nmice
   
    for i_area = 1:Nareas
        
        % grab the raw data;
        idx = strcmpi(mouseName, mice{i_mouse}) & strcmpi(brainArea, areas{i_area});
        if sum(idx)==0
         
            continue % this brain area was not tested in this mouse.
            
        elseif sum(idx) == 4
            
            % make sure that there are no duplicate files
            assert(size(unique(fileName(idx)),1)==4, 'ERROR: duplicate files found')
            
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
areas = {'PM','AL', 'RL', 'ERC'};
figure
set(gcf, 'position', [184    35   764   746])
for i_area = 1:numel(areas);
    
    tmp_volume = popdat.(areas{i_area}).totalVolume;
    tmp_counts = popdat.(areas{i_area}).cellCount;
    
    % density = {nAreas}[nLayers, nMice]
    volume{i_area} = tmp_volume;
    counts{i_area} = tmp_counts;
    density{i_area} = tmp_counts ./ tmp_volume;
    
    subplot(4, 3, (i_area-1)*3 + 1)
    plot(tmp_counts, '.-')
    xlabel('layer')
    ylabel('counts')
    title(areas{i_area})
    box off
    axis tight
    
    subplot(4, 3, (i_area-1)*3 + 2)
    plot(tmp_volume, '.-')
    xlabel('layer')
    ylabel('volume')
    title(areas{i_area})
    box off
    axis tight
    
    subplot(4, 3, (i_area-1)*3 + 3)
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
layers = [1:3]; % 1= L2/3, 2=L4, 3=L5, 4=L6
densityAcrossLayers = [];
for i_area = 1:numel(areas);
    tmp_density = density{i_area};
    tmp_density = sum(tmp_density(layers, :), 1); % there could be some NaNs, but I'm preserving them (and not using nansum) b/c I want to get an accurate N for the SEM...
    densityAcrossLayers(i_area,:) = tmp_density; % notice that the dim is [Nareas x Nmice]
end

figure, hold on,
plot(densityAcrossLayers, '.--')
nareas = numel(areas);
xbar = nanmean(densityAcrossLayers,2);
sem = nanstd(densityAcrossLayers,[], 2) ./ sqrt(sum(~isnan(densityAcrossLayers), 2));
errorbar(1:nareas, xbar, sem, 'k', 'linewidth', 3)
set(gca, 'xtick', 1:nareas, 'xTickLabel', areas)

% now plotting percent change relative to PM
figure, hold on,
prcnt_change_from_PM = bsxfun(@rdivide, densityAcrossLayers, densityAcrossLayers(1,:));
plot(prcnt_change_from_PM, '.--')
nareas = numel(areas);
xbar = nanmean(prcnt_change_from_PM, 2);
sem = nanstd(prcnt_change_from_PM,[], 2) ./ sqrt(sum(~isnan(prcnt_change_from_PM), 2));
errorbar(1:nareas, xbar, sem, 'k', 'linewidth', 3)
set(gca, 'xtick', 1:nareas, 'xTickLabel', areas)


%% PLOTTING ROUTINES: CELL DEPTH AND SIZE



% concatenate the data across mice and slices
areas = {'PM','AL', 'RL', 'ERC'};
for i_area = 1:numel(areas)
    
    % grab the raw data;
    idx = strcmpi(brainArea, areas{i_area});
    popdat.(areas{i_area}).cellDepth = cat(1, results(idx).cellDepths);
    popdat.(areas{i_area}).cellSize = cat(1, results(idx).cellSize);
    
end

% plot a histogram of cell sizes
figure
set(gcf, 'position', [104   164   871   601])
bigdataset = [];
for i_area = 1:numel(areas)
    
    subplot(numel(areas),1,i_area)
    tmp = popdat.(areas{i_area}).cellSize;
    bigdataset = cat(1, bigdataset, tmp);
    bins = linspace(0, 2500, 200);
    hist(tmp, bins)
    set(get(gca, 'children'), 'edgealpha', 0.1)
    ylims = get(gca, 'ylim');
    hold on,
    xbar = mean(tmp);
    plot([xbar, xbar], ylims, 'b', 'linewidth', 3)
    xlim([0, 2500])
    title(areas{i_area})
    ylabel('counts')
    xlabel('cell volume')
end


% plot cell sizes across all areas
figure
bins = linspace(0, 2500, 100);
hist(bigdataset, bins);
set(get(gca, 'children'), 'edgealpha', 0.1)
ylims = get(gca, 'ylim');
xlim([0, 2500])
title('All Areas Combined')
ylabel('counts')
xlabel('cell volume')


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
set(gcf, 'position', [108    88   863   663])
for i_area = 1:numel(areas)
    
    subplot(1,numel(areas),i_area)
    tmp_depth = popdat.(areas{i_area}).cellDepth;
    tmp_size = popdat.(areas{i_area}).cellSize;
    plot(tmp_depth, tmp_size, '.')
    ylim([0,2500])
    xlim([0,1000])
end





        