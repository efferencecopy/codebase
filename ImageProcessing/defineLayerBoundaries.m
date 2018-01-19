function defineLayerBoundaries


% import the data
[raw, mask, fpath] = io_importData;

% initialize the gui window display the raw image and the mask
gui_initialize(raw, mask, fpath);


end



function [raw, mask, fpath] = io_importData

    % find the path to the mask file (created in the opticaldissector.m program)
    [fnames, fpath] = uigetfile({['*'], ['all files']}, 'Select the mask file', 'multiselect', 'on');
    if ~iscell(fnames)
        tmpname = fnames;
        fnames = {};
        fnames{1} = tmpname;
    end

    
    % figure out which one was the .mat and the .tif
    numel(fnames)
    assert(any([1,2] == numel(fnames)), 'ERROR: please select a .tif and (optionally) a .mat file too')
    idx_raw = cellfun(@(x,y) ~isempty(regexpi(x,y)), fnames, repmat({'\.tif'}, size(fnames)));
    idx_mask = ~idx_raw;
    
    % load the raw image
    tiffstack = [fpath, fnames{idx_raw}];
    info = imfinfo(tiffstack, 'tif');
    raw.info = info(1);
    raw.info.Nframes = numel(info);
    [~, raw.info.shortName] = fileparts(raw.info.Filename);
    
    % unpack each frame of the raw image. raw.img will now have
    % dimensionality [x, y, Nframes, Ncolors]
    raw.img = nan(raw.info.Height, raw.info.Width, raw.info.Nframes, 3);
    for fr = 1:raw.info.Nframes
        raw.img(:,:,fr,:) = imread(tiffstack, 'tif', 'index', fr);
    end
    
    % load the mask
    if any(idx_mask)
        load([fpath, fnames{idx_mask}]); % loads a structure called "cellFillData"
        mask = cellFillData.mask;
    else
        nrows = size(raw.img, 1);
        ncols = size(raw.img, 2);
        nplanes = size(raw.img, 3);
        mask.img = zeros([nrows, ncols, nplanes]);
        mask.Ncells = 0;
    end
    
end

function io_exportData(varargin)
    
    % grab the user data
    cellFillData = get(gcf, 'userdata');
    
    % delete the actual raw image so that the file size isn't huge
    cellFillData.raw.img = [];
    
    % delete the handles (which cause matlab to complain when the .mat file
    % is opened for post-hoc anslysis)
    cellFillData.h = [];
    
    % cd to the appropriate directory
    dir_start = pwd;
    cd(cellFillData.fpath)
    
    % create a new .mat file name
    fname_out = ['layered_cellFillData_', cellFillData.raw.info.shortName, '.mat'];
    
    % export the data with the new file name
    uisave('cellFillData', fname_out)
    
    % cd back to the original directory, and reset the value of the
    % uibutton to zero
    cd(dir_start)
    set(varargin{1}, 'value', 0)
    
end


function gui_initialize(raw, mask, fpath)
    
    % for the main GUI window
    h.fig.main = figure;
    set(h.fig.main, 'position', [807 211 337 462]);
    
    % define a bunch of buttons to define the laminar boundaries
    h.gui.boundary1 = uicontrol('style', 'togglebutton',...
                         'units', 'normalized',...
                         'string', 'L1 to L2/3',...
                         'Position', [0.1, 0.9, 0.40, 0.05],...
                         'BackgroundColor', [0.7 0.7 0.7],...
                         'Callback', {@mask_defineBoundary, 'boundary1'});
    
    h.gui.boundary2 = uicontrol('style', 'togglebutton',...
                         'units', 'normalized',...
                         'string', 'L2/3 to L4',...
                         'Position', [0.1, 0.8, 0.40, 0.05],...
                         'BackgroundColor', [0.7 0.5 0.5],...
                         'Callback', {@mask_defineBoundary, 'boundary2'});
                     
    h.gui.boundary3 = uicontrol('style', 'togglebutton',...
                         'units', 'normalized',...
                         'string', 'L4 to L5',...
                         'Position', [0.1, 0.7, 0.40, 0.05],...                         
                         'BackgroundColor', [0.5 0.7 0.5],...
                         'Callback', {@mask_defineBoundary, 'boundary3'});
                     
    h.gui.boundary4 = uicontrol('style', 'togglebutton',...
                         'units', 'normalized',...
                         'string', 'L5 to L6',...
                         'Position', [0.1, 0.6, 0.40, 0.05],...                         
                         'BackgroundColor', [0.7 0.7 0.5],...
                         'Callback', {@mask_defineBoundary, 'boundary4'});
                     
    h.gui.boundary5 = uicontrol('style', 'togglebutton',...
                         'units', 'normalized',...
                         'string', 'Bottom of Cortex',...
                         'Position', [0.1, 0.5, 0.40, 0.05],...                         
                         'BackgroundColor', [0.5 0.5 0.7],...
                         'Callback', {@mask_defineBoundary, 'boundary5'});
    
   % define buttons to change the color chanel viewed on the raw data and
   % which focal plane is visible.
   nFocalPlanes = size(raw.img, 3);
   defaultFocalPlane = round(nFocalPlanes/2);
   h.gui.slider = uicontrol('style', 'slider',...
                         'units', 'normalized',...
                         'position', [0.2, 0.24, 0.6, 0.03],...
                         'Value', defaultFocalPlane,...
                         'max', nFocalPlanes,...
                         'min', 1,...
                         'SliderStep', [1./nFocalPlanes, 1./nFocalPlanes.*5]);
   
   % make a listener for the slider that will enable smooth control
   h.gui.sliderListen = addlistener(h.gui.slider, 'Value', 'PostSet',@(x,y) gui_updateRawImg(nan,nan,h.fig.main));
                    
   h.gui.colorChannel = uicontrol('style', 'listbox',...
                          'units', 'normalized',...
                          'position', [0.3, 0.33, 0.4, 0.13],...
                          'string', {'Red', 'Green', 'Blue'},...
                          'value', 2,...
                          'Callback', {@gui_updateRawImg, h.fig.main});
                      
                      
   % a button for exporting the data
   h.gui.export = uicontrol('style', 'togglebutton',...
                         'units', 'normalized',...
                         'string', 'Export Data',...
                         'Position', [0.25, 0.1, 0.50, 0.08],...
                         'Callback', {@io_exportData});
   
   
   % display the raw image and the mask
   h.fig.mask = figure;
   proj = max(mask.img,[],3);
   proj(proj>0) = 0.85;
   h.mask.proj = imshow(proj);
   hold on
   
   % add the lines for different boundaries
   h.mask.boundLine1 = plot(nan, nan, 'w--', 'linewidth', 3);
   h.mask.boundLine2 = plot(nan, nan, 'r--', 'linewidth', 3);
   h.mask.boundLine3 = plot(nan, nan, 'g--', 'linewidth', 3);
   h.mask.boundLine4 = plot(nan, nan, 'y--', 'linewidth', 3);
   h.mask.boundLine5 = plot(nan, nan, 'b--', 'linewidth', 3);
   set(gcf, 'position', [407 47 360 681])
   
   
   h.fig.raw = figure;
   h.raw.img = imshow(raw.img(:,:,defaultFocalPlane,2)./255);
   hold on,
   h.raw.boundLine1 = plot(nan, nan, 'w--', 'linewidth', 3);
   h.raw.boundLine2 = plot(nan, nan, 'r--', 'linewidth', 3);
   h.raw.boundLine3 = plot(nan, nan, 'g--', 'linewidth', 3);
   h.raw.boundLine4 = plot(nan, nan, 'y--', 'linewidth', 3);
   h.raw.boundLine5 = plot(nan, nan, 'b--', 'linewidth', 3);
   set(gcf, 'position', [25 42 354 681])
   
   % set the user data field in the main gui window
   udat.raw = raw;
   udat.mask = mask;
   udat.h = h;
   udat.fpath = fpath;
   set(h.fig.main, 'userdata', udat)
   
end

function gui_updateRawImg(varargin)

    % grab the udata
    udat = get(varargin{3}, 'userdata');
    
    % define the focal plane and color channel
    colorChannel = get(udat.h.gui.colorChannel, 'value');
    focalPlane = round(get(udat.h.gui.slider, 'value'));
    
    % update the raw image
    figure(udat.h.fig.raw)
    set(udat.h.raw.img, 'CData', udat.raw.img(:,:,focalPlane, colorChannel)./255);
    
end

function mask_defineBoundary(varargin)
    
    % grab the user data
    udat = get(gcf, 'userdata');
    
    % which boundary is this supposed to be?
    boundary = varargin{3};
    
    % grab two points using ginput. make sure to force the raw image to be
    % the current axes
    figure(udat.h.fig.raw);
    [x,y] = ginput(2);
    
    % define a slope and y intercept based off the 2 points
    m = diff(y)./diff(x);
    b = y(1) - (m.*x(1));
    
    % push these values into the appropriate field of the udat
    udat.raw.(boundary).m = m;
    udat.raw.(boundary).b = b;
    
    % reset the uibutton
    set(varargin{1}, 'value', 0);
    
    % reset the user data
    set(udat.h.fig.main, 'userdata', udat);
    
    % refresh the laminar boundaries on the mask and raw images
    mask_drawLaminarBoundaries(udat);

end


function mask_drawLaminarBoundaries(udat)
    
    
    boundType = {'boundary1', 'boundary2', 'boundary3', 'boundary4', 'boundary5'};
    boundHand = {'boundLine1', 'boundLine2', 'boundLine3', 'boundLine4', 'boundLine5'};
    for a = 1:numel(boundType)
        
       % make sure the slope and intercept are defined
       if ~isfield(udat.raw, boundType{a})
           continue
       end
       
       % define a domain, and calculate the yvals
       x = [1, size(udat.raw.img, 2)];
       m = udat.raw.(boundType{a}).m;
       b = udat.raw.(boundType{a}).b;
       y = m.*x + b;
       
       % update the line in the raw image
       figure(udat.h.fig.raw)
       set(udat.h.raw.(boundHand{a}), 'xdata', x, 'ydata', y)
       
       % update the line on the mask
       figure(udat.h.fig.mask)
       set(udat.h.mask.(boundHand{a}), 'xdata', x, 'ydata', y)
       
    end
    
    
    
    
    
    
end



