% load in all the data

fin

xlspath = [GL_DOCUPATH 'Other_workbooks', filesep, 'axon_anatomy_population_analysis.xlsx'];
[~, ~, raw] =xlsread(xlspath, 'All Images');
header = raw(1,:);
raw(1,:) = []; % remove the header row
Nfiles = size(raw, 1);

mouseName = raw(:, strcmpi(header, 'mouse'));
fileName = raw(:, strcmpi(header, 'filename'));
brainArea = raw(:, strcmpi(header, 'brain area'));
projectFile = raw(:, strcmpi(header, 'folder'));
cre_dependent = raw(:, strcmpi(header, 'cre dependent'));


% figure out which color channel to analyze
virus = raw(:, strcmpi(header, 'virus'));
color_channel = nan(Nfiles, 1);
for i_fid = 1:Nfiles
    if regexpi(virus{i_fid}, 'gfp')
        color_channel(i_fid) = 2;
    elseif regexpi(virus{i_fid}, 'tomato')
        color_channel(i_fid) = 1;
    else 
        error('Did not identify color channel')
    end
end

% figure out the genotype folder name
genotype = raw(:, strcmpi(header, 'genotype'));
genotype_folder = repmat({'null'}, Nfiles, 1);
layer_line = repmat({'null'}, Nfiles, 1);
for i_fid = 1:Nfiles
    if regexpi(genotype{i_fid}, 'EMX')
        genotype_folder{i_fid} = 'GIN';
        layer_line{i_fid} = 'emx';
    elseif regexpi(genotype{i_fid}, 'PV-Cre')
        genotype_folder{i_fid} = 'PV-Cre';
        layer_line{i_fid} = 'none';
    elseif regexpi(genotype{i_fid}, 'SOM-Cre')
        genotype_folder{i_fid} = 'SOM-Cre';
        layer_line{i_fid} = 'none';
    elseif regexpi(genotype{i_fid}, 'Calb1'), 
        genotype_folder{i_fid} = 'Calb1-cre';
        layer_line{i_fid} = 'calb1';
    elseif regexpi(genotype{i_fid}, 'Tlx3')
        genotype_folder{i_fid} = 'TLX-cre';
        layer_line{i_fid} = 'tlx3';
    else
        error('Did not identify the genotype correctly')
    end
    fprintf('%s, %s\n', layer_line{i_fid}, genotype{i_fid}) 
end



dat = struct(); % a cell array of structures
for i_fid = 1:Nfiles
    
    fprintf('file %d of %d: %s\n', i_fid, Nfiles, fileName{i_fid})
    
    % make sure the slice is from one of the defined arease
    slice_hva = lower(brainArea{i_fid});
    if ~any(strcmp(slice_hva, {'al', 'lm', 'am', 'pm', 'erc'}))
        continue
    end
    if  ~isfield(dat, mouseName{i_fid})
        hvastruct = struct('images', {{}}, 'profiles', {{}}, 'layers', struct());
        dat.(mouseName{i_fid}) =  struct('al', hvastruct,...
                                         'pm', hvastruct,...
                                         'lm', hvastruct,...
                                         'am', hvastruct,...
                                         'erc', hvastruct);
    end
    
    
    % save some meta data
    dat.(mouseName{i_fid}).cre_dependent = cre_dependent{i_fid};
    
    % derive the file name for the .tiff file from the layered_cellfill_...
    % name.
    fname_tiff = fileName{i_fid};
    prefixes = {'layered_', 'cellFillData_'};
    for i_pfx = 1:numel(prefixes)
        prefix_stop_idx = regexpi(fname_tiff, prefixes{i_pfx}, 'end');
        if prefix_stop_idx
            fname_tiff = fname_tiff(prefix_stop_idx+1:end);
        end
    end
    
    % fix a few file names
    if strcmpi(fname_tiff(end-2:end), '_ch')
        fname_tiff = fname_tiff(1:end-3);
    end
      
    % find and load the file
    startingPath = [GL_DATPATH(1:end-5), projectFile{i_fid}, filesep, genotype_folder{i_fid}, filesep];
    optpath = [startingPath, mouseName{i_fid}, filesep, 'confocal', filesep, slice_hva];
    tiff_path = findfile(fname_tiff, optpath, '.tif'); % recursively look for it
    tiff_info = imfinfo(tiff_path);
    n_frames = numel(tiff_info);
    n_pix_x = tiff_info(1).Width;
    n_pix_y = tiff_info(1).Height;
    
    % unpack the image and analyze only the axon channel
    tiff_raw = nan(n_pix_y, n_pix_x, n_frames);
    for i_frame = 1:n_frames
        frame_raw = imread(tiff_path, 'tif', i_frame);
        tiff_raw(:,:,i_frame) = frame_raw(:,:,color_channel(i_fid));
    end
    
    % get the intensity profile
    intensity_profile = permute(mean(mean(tiff_raw, 1), 2), [3,1,2]);
    threshold_intensity = max(intensity_profile) ./ 2;
    
    % store the slices above threshold
    l_above_threshold = intensity_profile >= threshold_intensity;
    tiff_raw = double(tiff_raw(:, :, l_above_threshold));
    idx_img = numel(dat.(mouseName{i_fid}).(slice_hva).images) + 1;
    dat.(mouseName{i_fid}).(slice_hva).images{idx_img} = tiff_raw;
    dat.(mouseName{i_fid}).(slice_hva).profiles{idx_img} = sum(sum(tiff_raw, 2), 3);
    
    % store some metatdata about layer boundaries
    meta_path = findfile(fileName{i_fid}, optpath, '.mat'); % recursively look for it
    meta_data = load(meta_path); % load the data
    meta_data = meta_data.cellFillData.raw;
    
    % define the mean position of the layer boundaries
    midpoint = round(meta_data.info.Width ./ 2);
    L1_L23 = round((meta_data.boundary1.m .* midpoint) + (meta_data.boundary1.b));
    L23_L4 = round((meta_data.boundary2.m .* midpoint) + (meta_data.boundary2.b));
    L4_L5 =  round((meta_data.boundary3.m .* midpoint) + (meta_data.boundary3.b));
    L5_L6 =  round((meta_data.boundary4.m .* midpoint) + (meta_data.boundary4.b));
    bottom = round((meta_data.boundary5.m .* midpoint) + (meta_data.boundary5.b));
    
    dat.(mouseName{i_fid}).(slice_hva).layers.L1_L23{idx_img} = L1_L23;
    dat.(mouseName{i_fid}).(slice_hva).layers.L23_L4{idx_img} = L23_L4;
    dat.(mouseName{i_fid}).(slice_hva).layers.L4_L5{idx_img} = L4_L5;
    dat.(mouseName{i_fid}).(slice_hva).layers.L5_L6{idx_img} = L5_L6;
    dat.(mouseName{i_fid}).(slice_hva).layers.bottom{idx_img} = bottom;
    
end


%% Just plot the data, one subplot per slice



all_hvas = {'lm', 'al', 'pm', 'am'};
mouse_names = fieldnames(dat);
for mouse = mouse_names'
    mouse = mouse{1};
    
    nrows = 0;
    for i_hva = 1:numel(all_hvas)
        nrows = max(nrows, numel(dat.(mouse).(all_hvas{i_hva}).images));
    end
    
    hf = figure;
    hf.Name = mouse;
    for i_hva = 1:numel(all_hvas)
        n_imgs = numel(dat.(mouse).(all_hvas{i_hva}).images);
        for i_img = 1:n_imgs
            plt_idx = sub2ind([nrows, numel(all_hvas)], i_img, i_hva);
            ha = subplot(numel(all_hvas), nrows,plt_idx);
            img_raw = dat.(mouse).(all_hvas{i_hva}).images{i_img};
            hi = image(max(img_raw, [], 3));
            colormap((gray(256)))
            ha.XTick = [];
            ha.YTick = [];
            hold on,
            l1_yy = dat.(mouseName{i_fid}).(slice_hva).layers.L1_L23{i_img};
            plot([1, size(img_raw,2)], [l1_yy, l1_yy], 'g--')
            if i_img == 1
                ha.YLabel.String = all_hvas{i_hva};
            end
        end
    end    
end


%% Plot the profile of fluorescence intensity (one for each slice)

close all

mouse_names = fieldnames(dat);
for mouse = mouse_names'
    mouse = mouse{1};
    hf = figure;
    hf.Name = mouse;
    all_hvas = {'al', 'pm', 'am', 'lm'};
    max_y = 0;
    max_x = 0;
    min_x = inf;
    ha = [];
    for i_hva = 1:numel(all_hvas)
        n_imgs = numel(dat.(mouse).(all_hvas{i_hva}).profiles);
        ha(i_hva) = subplot(numel(all_hvas),1,i_hva); hold on,
        ylabel(all_hvas{i_hva});
        for i_img = 1:n_imgs
            
            % raw profile
            profile_raw = dat.(mouse).(all_hvas{i_hva}).profiles{i_img};
            
            % layer boundaries
            l1_idx = dat.(mouse).(all_hvas{i_hva}).layers.L1_L23{i_img};
            l6_idx = dat.(mouse).(all_hvas{i_hva}).layers.L5_L6{i_img};
            bottom = min(dat.(mouse).(all_hvas{i_hva}).layers.bottom{i_img}, numel(profile_raw));
            offset = mean(profile_raw(l6_idx:bottom));
            
            % adjust profile for differences in bkgnd
            %profile_raw = profile_raw - offset;
            
            % plot, but define x=0 as the L1/L23 boundary
            xx = [1:numel(profile_raw)] - l1_idx;
            hp = plot(xx, profile_raw);
            
            % store the max axis vals
            max_y = max(max_y, max(profile_raw));
            max_x = max(max_x, xx(end));
            min_x = min(min_x, xx(1));
            
        end
    end
    set(ha, 'YLim', [0, max_y], 'XLim', [min_x, max_x])
    
end



%% Plot the average profile of fluorescence intensity

close all
clc

NORMALIZE = false;
DELETE_PIA = false;
all_hvas = {'lm', 'al', 'pm', 'am'};
hva_clrs = {'r', 'm', 'b', 'c'};

mouse_names = fieldnames(dat);
for mouse = mouse_names'
    mouse = mouse{1};
    hf = figure; hold on,
    legend_text = {};
    for i_hva = 1:numel(all_hvas)
        if isempty(dat.(mouse).(all_hvas{i_hva}).images)
            continue % no data to plot
        end
        legend_text{end+1} = all_hvas{i_hva};
        xx_profile = {};
        yy_profile = {};
        common_right_edge = inf;
        common_left_edge = -inf;
        n_imgs = numel(dat.(mouse).(all_hvas{i_hva}).images);
        for i_img = 1:n_imgs
            % raw profile
            yy_profile{i_img} = dat.(mouse).(all_hvas{i_hva}).profiles{i_img};
            
            % layer boundaries
            l1_idx = dat.(mouse).(all_hvas{i_hva}).layers.L1_L23{i_img};
            l6_idx = dat.(mouse).(all_hvas{i_hva}).layers.L5_L6{i_img};
            bottom = min(dat.(mouse).(all_hvas{i_hva}).layers.bottom{i_img}, numel(yy_profile{i_img}));
            offset = mean(yy_profile{i_img}(l6_idx:bottom));
            
            % adjust profile for differences in bkgnd
            %yy_profile{i_img} = yy_profile{i_img} - offset;
            
            % plot, but define x=0 as the L1/L23 boundary
            xx_profile{i_img} = [1:numel(yy_profile{i_img})] - l1_idx;
            
            common_left_edge = max(common_left_edge, xx_profile{i_img}(1));
            common_right_edge = min(common_right_edge, xx_profile{i_img}(end));
            
        end
        
        
        new_yy = cellfun(@(x,y) x(find(y==common_left_edge) : find(y==common_right_edge)), yy_profile, xx_profile, 'uniformoutput', false);
        new_yy = cat(2, new_yy{:});
        avg_yy = mean(new_yy, 2);
        new_xx = [common_left_edge : common_right_edge];
        if NORMALIZE
            avg_yy = avg_yy ./ max(avg_yy(new_xx>0));
        end
        if DELETE_PIA
            avg_yy(avg_yy < -0.20) = NaN;
        end
        plot(new_xx, avg_yy, '-', 'color', hva_clrs{i_hva}, 'linewidth', 2)
    end
    hf.Name = mouse;
    ylabel(all_hvas{i_hva});
    legend(legend_text, 'location', 'best')
end






