function defineLayerBoundaries

% game plan
%
% Make gui with some buttons
%   * change color chanel
%   * export data
%   * scroll through focal planes

% import the data
[raw, mask] = io_importData;

% initialize the gui window display the raw image and the mask
gui_initialize(raw, mask);


end



function [raw, mask] = io_importData

    % find the path to the mask file (created in the opticaldissector.m program)
    [fnames, fpath] = uigetfile({['*'], ['all files']}, 'Select the mask file', 'multiselect', 'on');

    
    % figure out which one was the .mat and the .tif
    assert(numel(fnames) == 2, 'ERROR: please select a .mat AND a .tiff file')
    idx_raw = cellfun(@(x,y) ~isempty(regexpi(x,y)), fnames, repmat({'.tif'}, size(fnames)));
    idx_mask = ~idx_raw;
    
    % load the mask
    load([fpath, fnames{idx_mask}]); % loads a structure called "cellFillData"
    mask = cellFillData.mask;
    
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
    

end

function gui_initialize(raw, mask)
    
    % for the main GUI window
    h.fig.main = figure;
    set(h.fig.main, 'position', [1000 151 337 462]);
    
    % define a bunch of buttons to define the laminar boundaries
    h.gui.boundary1 = uicontrol('style', 'togglebutton',...
                         'units', 'normalized',...
                         'string', 'L1 to L2/3',...
                         'Position', [0.1, 0.9, 0.40, 0.05],...
                         'BackgroundColor', [0.7 0.7 0.7],...
                         'Callback', {@mask_defineBoundary});
    
    h.gui.boundary2 = uicontrol('style', 'togglebutton',...
                         'units', 'normalized',...
                         'string', 'L2/3 to L4',...
                         'Position', [0.1, 0.78, 0.40, 0.05],...
                         'BackgroundColor', [0.7 0.5 0.5],...
                         'Callback', {@mask_defineBoundary});
                     
    h.gui.boundary3 = uicontrol('style', 'togglebutton',...
                         'units', 'normalized',...
                         'string', 'L4 to L5',...
                         'Position', [0.1, 0.66, 0.40, 0.05],...                         
                         'BackgroundColor', [0.5 0.7 0.5],...
                         'Callback', {@mask_defineBoundary});
                     
    h.gui.boundary4 = uicontrol('style', 'togglebutton',...
                         'units', 'normalized',...
                         'string', 'Bottom of Cortex',...
                         'Position', [0.1, 0.54, 0.40, 0.05],...                         
                         'BackgroundColor', [0.5 0.5 0.7],...
                         'Callback', {@mask_defineBoundary});
    
   % define buttons to change the color chanel viewed on the raw data and
   % which focal plane is visible.
   nFocalPlanes = size(raw.img, 3);
   defaultFocalPlane = round(nFocalPlanes/2);
   h.gui.slider = uicontrol('style', 'slider',...
                         'units', 'normalized',...
                         'position', [0.2, 0.24, 0.6, 0.03],...
                         'Callback', {@gui_updateRawImg},...
                         'Value', defaultFocalPlane,...
                         'max', nFocalPlanes,...
                         'min', 1,...
                         'SliderStep', [1./nFocalPlanes, 1./nFocalPlanes.*5]);
                    
   h.gui.colorChannel = uicontrol('style', 'listbox',...
                          'units', 'normalized',...
                          'position', [0.3, 0.33, 0.4, 0.13],...
                          'string', {'Red', 'Green', 'Blue'},...
                          'value', 2,...
                          'Callback', {@gui_updateRawImg});
                      
                      
   % a button for exporting the data
   h.gui.boundary4 = uicontrol('style', 'togglebutton',...
                         'units', 'normalized',...
                         'string', 'Export Data',...
                         'Position', [0.25, 0.1, 0.50, 0.08],...
                         'Callback', {@io_exportData});
   
   
   % display the raw image and the mask
   h.fig.mask = figure;
   set(gcf, 'position', [546 19 289 660])
   proj = max(mask.img,[],3);
   h.mask.proj = imshow(proj./max(proj(:)));
   hold on
   
   % add the lines for different boundaries
   h.mask.boundLine1 = plot(nan, nan, 'w--', 'linewidth', 3);
   h.mask.boundLine2 = plot(nan, nan, 'r--', 'linewidth', 3);
   h.mask.boundLine3 = plot(nan, nan, 'g--', 'linewidth', 3);
   h.mask.boundLine4 = plot(nan, nan, 'b--', 'linewidth', 3);
   
   
   h.fig.raw = figure;
   set(gcf, 'position', [174    16   289   660])
   h.raw.img = imshow(raw.img(:,:,defaultFocalPlane,2)./255);
   hold on,
   h.raw.boundLine1 = plot(nan, nan, 'w--', 'linewidth', 3);
   h.raw.boundLine2 = plot(nan, nan, 'r--', 'linewidth', 3);
   h.raw.boundLine3 = plot(nan, nan, 'g--', 'linewidth', 3);
   h.raw.boundLine4 = plot(nan, nan, 'b--', 'linewidth', 3);
   
   
   % set the user data field in the main gui window
   udat.raw = raw;
   udat.mask = mask;
   udat.h = h;
   set(h.fig.main, 'userdata', udat)
   
end

function gui_updateRawImg(varargin)

    % grab the udata
    udat = get(gcf, 'userdata');
    
    % define the focal plane and color channel
    colorChannel = get(udat.h.gui.colorChannel, 'value');
    focalPlane = round(get(udat.h.gui.slider, 'value'));
    
    % update the raw image
    figure(udat.h.fig.raw)
    set(udat.h.raw.img, 'CData', udat.raw.img(:,:,focalPlane, colorChannel)./255);
    
end






