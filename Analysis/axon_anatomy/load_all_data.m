function dat = load_all_data()

global GL_DOCUPATH
global GL_DATPATH

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
        genotype_folder{i_fid} = 'PVCre';
        layer_line{i_fid} = 'none';
    elseif regexpi(genotype{i_fid}, 'SOM-Cre')
        genotype_folder{i_fid} = 'SOMCre';
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
        hvastruct = struct('images', {{}}, 'profiles', {{}}, 'n_zplanes', {{}}, 'z_profile', {{}}, 'layers', struct());
        dat.(mouseName{i_fid}) =  struct('al', hvastruct,...
                                         'pm', hvastruct,...
                                         'lm', hvastruct,...
                                         'am', hvastruct,...
                                         'rl', hvastruct,...
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
    n_pix = size(tiff_raw,2) .* size(tiff_raw,3);
    dat.(mouseName{i_fid}).(slice_hva).profiles{idx_img} = sum(sum(tiff_raw, 2), 3) ./ n_pix;
    dat.(mouseName{i_fid}).(slice_hva).n_zplanes{idx_img} = sum(l_above_threshold);
    dat.(mouseName{i_fid}).(slice_hva).z_profile{idx_img} = intensity_profile(l_above_threshold);
    
    % store some metatdata about layer boundaries
    meta_path = findfile(fileName{i_fid}, optpath, '.mat'); % recursively look for it
    meta_data = load(meta_path); % load the data
    meta_data = meta_data.cellFillData.raw;
    
    % define the mean position of the layer boundaries
    xx = 1:meta_data.info.Width;
    L1_L23 = round((meta_data.boundary1.m .* xx) + (meta_data.boundary1.b));
    L23_L4 = round((meta_data.boundary2.m .* xx) + (meta_data.boundary2.b));
    L4_L5 =  round((meta_data.boundary3.m .* xx) + (meta_data.boundary3.b));
    L5_L6 =  round((meta_data.boundary4.m .* xx) + (meta_data.boundary4.b));
    bottom = round((meta_data.boundary5.m .* xx) + (meta_data.boundary5.b));
    
    dat.(mouseName{i_fid}).(slice_hva).layers.L1_L23{idx_img} = L1_L23;
    dat.(mouseName{i_fid}).(slice_hva).layers.L23_L4{idx_img} = L23_L4;
    dat.(mouseName{i_fid}).(slice_hva).layers.L4_L5{idx_img} = L4_L5;
    dat.(mouseName{i_fid}).(slice_hva).layers.L5_L6{idx_img} = L5_L6;
    dat.(mouseName{i_fid}).(slice_hva).layers.bottom{idx_img} = bottom;
    
end
