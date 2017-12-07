% a quick, hasty analysis of Axon ramifications in the HVAs
% using Amy's data from GIN x EMX mice.
%
% flex.ChR2 injected into V1

% loop through the GIN x EMX mice
% pull out "cropped" images from AL and PM, extract the red channel
% for now, just overlay the data and look at it
fin

CELLTYPE = 'PVcre';
COLOR_CHANNEL = 2;
startingPath = [GL_DATPATH(1:end-5), 'SOM_PV_Density', filesep];

xlspath = [GL_DOCUPATH 'Other_workbooks', filesep, 'Interneuron_density_analysis.xlsx'];
[~, ~, raw] =xlsread(xlspath, CELLTYPE);
header = raw(1,:);
raw(1,:) = []; % remove the header row

mouseName = raw(:, strcmpi(header, 'mouse'));
fileName = raw(:, strcmpi(header, 'filename'));
brainArea = raw(:, strcmpi(header, 'brain area'));
Nfiles = size(fileName, 1);

dat = struct(); % a cell array of structures
for i_fid = 1:Nfiles
    
    % make sure the slice is from AL or PM, store the brain area
    slice_hva = lower(brainArea{i_fid});
    if ~any(strcmp(slice_hva, {'al', 'pm', 'erc'}))
        continue
    end
    
    % define a file name
    fname_layered = fileName{i_fid};
    prefix = 'layered_cellFillData_';
    assert(strncmp(fname_layered, prefix, numel(prefix)), 'ERROR: incorrect prefix')
    fname_tiff = fname_layered(numel(prefix)+1:end);
    if strcmpi(fname_tiff(end-2:end), '_ch')
        fname_tiff = fname_tiff(1:end-3);
    end
    if strcmpi(fname_tiff(end-3:end), '.mat')
        fname_tiff = fname_tiff(1:end-4);
    end
        
    
    % find and load the file
    optpath = [startingPath, mouseName{i_fid}, filesep, 'confocal', filesep, slice_hva];
    tiff_path = findfile(fname_tiff, optpath, '.tif'); % recursively look for it
    tiff_info = imfinfo(tiff_path);
    n_frames = numel(tiff_info);
    n_pix_x = tiff_info(1).Width;
    n_pix_y = tiff_info(1).Height;
    
    tiff_raw = nan(n_pix_y, n_pix_x, n_frames);
    for i_frame = 1:n_frames
        frame_raw = imread(tiff_path, 'tif', i_frame);
        tiff_raw(:,:,i_frame) = frame_raw(:,:,COLOR_CHANNEL);
    end
    
    % get the intensity profile
    intensity_profile = permute(mean(mean(tiff_raw, 1), 2), [3,1,2]);
    threshold_intensity = max(intensity_profile) ./ 2;
    
    % store the slices above threshold
    if  ~isfield(dat, mouseName{i_fid})
        hvastruct = struct('images', {{}}, 'profiles', {{}}, 'layers', struct());
        dat.(mouseName{i_fid}) =  struct('al', hvastruct,...
                                         'pm', hvastruct,...
                                         'erc', hvastruct);
    end
    l_above_threshold = intensity_profile >= threshold_intensity;
    tiff_raw = double(tiff_raw(:, :, l_above_threshold));
    idx = numel(dat.(mouseName{i_fid}).(slice_hva).images) + 1;
    dat.(mouseName{i_fid}).(slice_hva).images{idx} = tiff_raw;
    dat.(mouseName{i_fid}).(slice_hva).profiles{idx} = sum(sum(tiff_raw, 2), 3);
    
    % store some metatdata about layer boundaries
    meta_path = findfile(fname_layered, optpath, '.mat'); % recursively look for it
    meta_data = load(meta_path); % load the data
    meta_data = meta_data.cellFillData.raw;
    
    % define the mean position of the layer boundaries
    midpoint = round(meta_data.info.Width ./ 2);
    L1_L23 = round((meta_data.boundary1.m .* midpoint) + (meta_data.boundary1.b));
    L23_L4 = round((meta_data.boundary2.m .* midpoint) + (meta_data.boundary2.b));
    L4_L5 =  round((meta_data.boundary3.m .* midpoint) + (meta_data.boundary3.b));
    L5_L6 =  round((meta_data.boundary4.m .* midpoint) + (meta_data.boundary4.b));
    bottom = round((meta_data.boundary5.m .* midpoint) + (meta_data.boundary5.b));
    
    dat.(mouseName{i_fid}).(slice_hva).layers.L1_L23 = L1_L23;
    dat.(mouseName{i_fid}).(slice_hva).layers.L23_L4 = L23_L4;
    dat.(mouseName{i_fid}).(slice_hva).layers.L4_L5 = L4_L5;
    dat.(mouseName{i_fid}).(slice_hva).layers.L5_L6 = L5_L6;
    dat.(mouseName{i_fid}).(slice_hva).layers.bottom = bottom;
    
end


%% Just plot the data, one subplot per slice



all_hvas = {'al', 'pm', 'erc'};
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
    all_hvas = {'al', 'pm', 'erc'};
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
            l1_idx = dat.(mouse).(all_hvas{i_hva}).layers.L1_L23;
            l6_idx = dat.(mouse).(all_hvas{i_hva}).layers.L5_L6;
            bottom = min(dat.(mouse).(all_hvas{i_hva}).layers.bottom, numel(profile_raw));
            offset = mean(profile_raw(l6_idx:bottom));
            
            % adjust profile for differences in bkgnd
            profile_raw = profile_raw - offset;
            
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

NORMALIZE = true;
DELETE_PIA = false;
all_hvas = {'al', 'pm'};
hva_clrs = {'r', 'b', 'k'};

mouse_names = fieldnames(dat);
for mouse = mouse_names'
    mouse = mouse{1};
    hf = figure; hold on,
    for i_hva = 1:numel(all_hvas)
        xx_profile = {};
        yy_profile = {};
        common_right_edge = inf;
        common_left_edge = -inf;
        n_imgs = numel(dat.(mouse).(all_hvas{i_hva}).images);
        for i_img = 1:n_imgs
            % raw profile
            yy_profile{i_img} = dat.(mouse).(all_hvas{i_hva}).profiles{i_img};
            
            % layer boundaries
            l1_idx = dat.(mouse).(all_hvas{i_hva}).layers.L1_L23;
            l6_idx = dat.(mouse).(all_hvas{i_hva}).layers.L5_L6;
            bottom = min(dat.(mouse).(all_hvas{i_hva}).layers.bottom, numel(yy_profile{i_img}));
            offset = mean(yy_profile{i_img}(l6_idx:bottom));
            
            % adjust profile for differences in bkgnd
            yy_profile{i_img} = yy_profile{i_img} - offset;
            
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
end






