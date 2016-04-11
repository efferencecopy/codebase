function opticaldissector(tiffstack, colorchannel)

    %  open the tiff stack
    %
    % if 'tiffstack' is a path, open the stack. if not, open a ui control to
    % navigate to the data file
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~exist('tiffstack', 'var') || isempty(tiffstack)
        [filename, pathname] = uigetfile({'*.tif;*.tiff'}, 'Pick a tiff stack');
        tiffstack = [pathname, filesep, filename];
    end
    if ~exist('colorchannel', 'var')
        colorchannel = [];
    end
    raw = io_unpacktiff(tiffstack, colorchannel); % import the data.



    %
    % setup the gui and update the mask
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%
    gui_initialize(raw)
    mask_update()


    % clear the command window b/c the gui is going to print out a bunch of
    % text
    clc

end % main function



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       SUPPORTING  sub-FUNCTIONS           
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function raw = io_unpacktiff(tiffstack, defaultcolor)
    
    % grab and store info for the tiff stack
    info = imfinfo(tiffstack, 'tif');
    raw.info = info(1);
    raw.info.Nframes = numel(info);
    
    % make a short hand name that will be used for autosave and exporting
    % the data
    [~, raw.info.shortName] = fileparts(raw.info.Filename);
    
    
    % this gui only works with monochrome images, if the tiff has more than
    % one channel, force the user to select just one channel
    if ~isempty(defaultcolor)
        colorchannel = defaultcolor;
    else
        str = {'Red', 'Green', 'Blue'};
        colorchannel = listdlg('PromptString','Select a color channel:',...
                                'SelectionMode','single',...
                                'ListString',str,...
                                'ListSize', [160 100]);
    end
    
    
    % unpack each frame of the image
    raw.img = nan(raw.info.Height, raw.info.Width, raw.info.Nframes);
    for fr = 1:raw.info.Nframes
        tmp = imread(tiffstack, 'tif', 'index', fr);
        raw.img(:,:,fr) = tmp(:,:,colorchannel);
    end
    
end

function io_exportData(varargin)
    
    % grab the contents of the user data
    cellFillData = get(gcf, 'userdata');
    cellFillData.h = [];
    
    % save the data interactively
    uisave('cellFillData', ['cellFillData_', cellFillData.raw.info.shortName, '.mat'])
    
end

function io_importSavedMask(varargin)
    
    % grab the user data
    udat = get(gcf, 'userdata');
    set(udat.h.importMask, 'value', 0);
    
    % allow the user to specify a pre-existing mask file, which will be a
    % .mat file
    [filename, pathname] = uigetfile({'*.mat'}, 'Pick a mask file');
    if ~filename
        return
    end
    
    % loads a structure called "cellFillData" which should contain the
    % saved mask;
    load([pathname, filesep, filename]); 
    
    % grab the saved mask and do a quick check
    savedMask = cellFillData.mask.img;
    savedNCells = cellFillData.mask.Ncells;
    sz_saved = size(savedMask);
    sz_current = size(udat.mask.img);
    assert(all(sz_saved == sz_current), 'ERROR: Saved mask is the wrong size');
    
    % overwrite the current mask data
    udat.mask.img = savedMask;
    udat.mask.Ncells = savedNCells;
    
    % set the udat and replot the main gui window
    set(udat.h.main, 'userdata', udat);
    gui_updateImages(udat.h.main);
    
    % display some text
    fprintf('Imported pre-existing mask with %d cells\n', udat.mask.Ncells);
    
    
end



function gui_initialize(raw)
    
    % the top level figure
    h.main = figure;
    set(h.main, 'name', sprintf('Optical Dissector'),...
                'units', 'normalized',...
                'position', [0.3628 0.0394 0.5998 0.8900],...
                'windowkeypressfcn', {@gui_keypress})
            
    % useful constants    
    lims = [raw.info.MinSampleValue(1) raw.info.MaxSampleValue(1)];
    raw.clrmap = colormap('gray');
    
    % add an axis for the main focal plane
    gl.Zplane = round(raw.info.Nframes/2); % start in the middle
    h.ax_curFrame = axes('position', [0, 0.02, .4, .96]);
    h.img_curFrame = imshow(raw.img(:,:,gl.Zplane), lims,...
                            'parent', h.ax_curFrame,...
                            'colormap', raw.clrmap);
    set(h.img_curFrame, 'ButtonDownFcn', {@gui_mainImgClick});
    
                        
    % add an axis for a slice through the vertical dimension
    h.ax_sliceVert = axes('position', [.32, 0.02, .2, .96]);
    gl.Xplane = round(raw.info.Width/2); % default starting slice is in the middle
    gl.sliceExapandScalar = ceil(50./raw.info.Nframes);
    vertSlice = squeeze(raw.img(:, gl.Xplane, :));
    vertSlice = repmat(vertSlice, gl.sliceExapandScalar, 1);
    vertSlice = reshape(vertSlice, raw.info.Height, []);
    h.img_sliceVert = imshow(vertSlice, lims,...
                            'parent', h.ax_sliceVert,...
                            'colormap', raw.clrmap);
    set(h.img_sliceVert, 'ButtonDownFcn', {@gui_vertSliceImgClick});

    
    % add an axis for a slice through the horizontal dimension
    h.ax_sliceHoriz = axes('position', [0.47, 0.50, .3, .31]);
    gl.Yplane = round(raw.info.Height/2); % default starting slice is in the middle
    horizSlice = squeeze(raw.img(gl.Yplane,:, :));
    horizSlice = repmat(horizSlice, gl.sliceExapandScalar, 1);
    horizSlice = reshape(horizSlice, raw.info.Width, [])';
    h.img_sliceHoriz = imshow(horizSlice, lims,...
                            'parent', h.ax_sliceHoriz,...
                            'colormap', raw.clrmap);
    set(h.img_sliceHoriz, 'ButtonDownFcn', {@gui_horzSliceImgClick});
    
    
    % add a button to mark a filled region as multiple cells
    h.multipleCells = uicontrol('style', 'togglebutton',...
                                'units', 'normalized',...
                                'string', 'Multiple Cells',...
                                'Position', [0.47, 0.265, 0.2 0.05],...
                                'Callback', {@mask_markAsMultiple});
    
    
    % add a button to save the data and dump the result onto the command
    % window
    h.export = uicontrol('style', 'togglebutton',...
                         'units', 'normalized',...
                         'string', 'Export Data',...
                         'Position', [0.72, 0.265, 0.20, 0.05],...
                         'Callback', {@io_exportData});
                     
    % add a button to import a previously saved mask file
    h.importMask = uicontrol('style', 'togglebutton',...
                         'units', 'normalized',...
                         'string', 'Import Mask',...
                         'Position', [0.47, 0.20, 0.20, 0.05],...
                         'Callback', {@io_importSavedMask});
                     
    % add a button to quickly check your work
    h.quickCheck = uicontrol('style', 'togglebutton',...
                         'units', 'normalized',...
                         'string', 'Quick Check',...
                         'Position', [0.72, 0.20, 0.20, 0.05],...
                         'Callback', {@mask_quickCheck}); 
                     
    % add a little text box for an online cell counter
    h.cellCounter = uicontrol('style', 'text',...
                              'units', 'normalized',...
                              'string', 'Cell Count: 0',...
                              'fontsize', 14,...
                              'position', [0.58, 0.12, 0.2, 0.05],...
                              'BackgroundColor', get(h.main, 'color'));
    
                        
    % package the relevant info into the user data field
    udat.gl = gl;
    udat.h = h;
    udat.raw = raw;
    set(h.main, 'userdata', udat);
    
end

function gui_updateImages(h_main)
    
    % grab stuff from the user data field
    if exist('h_main', 'var')
        udat = get(h_main, 'userdata');
    else
        udat = get(gcf, 'userdata');
    end
    Xidx = udat.gl.Xplane;
    Yidx = udat.gl.Yplane;
    Zidx = udat.gl.Zplane;
    maxX = udat.raw.info.Width;
    maxY = udat.raw.info.Height;
    maxdac = udat.raw.info.MaxSampleValue(1);
    
    % make an RGB version of the image, but scale certian pixels according
    % to the cell mask
    ch = {'R', 'G', 'B'};
    clr = [1 0.7 0.7];
    maskvect = udat.mask.img > 0;
    for a = 1:3
        tmp = udat.raw.img./maxdac;
        tmp = tmp(:);
        tmp(maskvect(:)) = tmp(maskvect(:)).* clr(a);
        fov.(ch{a}) = reshape(tmp, maxY, maxX, []);
    end
        
    
    % new main image focal plane
    tmp = cat(3, fov.R(:,:, Zidx), fov.G(:,:,Zidx), fov.B(:,:,Zidx));
    set(udat.h.img_curFrame, 'CData', tmp)
     
    % new vertical slice
    for a = 1:3
        tmp = squeeze(fov.(ch{a})(:, Xidx, :));
        tmp = repmat(tmp, udat.gl.sliceExapandScalar, 1);
        tmp = reshape(tmp, maxY, []);
        vert(:,:,a) = tmp;
    end
    set(udat.h.img_sliceVert, 'CData', vert)
    
    % new horizontal slice    
    for a = 1:3
        tmp = squeeze(fov.(ch{a})(Yidx,:, :));
        tmp = repmat(tmp, udat.gl.sliceExapandScalar, 1);
        tmp = reshape(tmp, maxX, [])';
        horiz(:,:,a) = tmp;
    end
    set(udat.h.img_sliceHoriz, 'CData', horiz)

    
    % update the online counter
    msg = sprintf('Cell Count:  %d', udat.mask.Ncells);
    set(udat.h.cellCounter, 'string', msg);
    
end

function gui_getNewFocalPlane(h_slice, ~) % currently unused.
    
    % grab stuff from the user data field
    udat = get(gcf, 'userdata');
    
    
    % we don't know which slice the user clicked on, but it shouldn't
    % matter b/c dimensionality of each slice is the same in the z
    % direction.
    clickmode = get(gcf, 'selectiontype');
    if strcmp(clickmode, 'alt') % a right click
        pt = get(get(h_slice, 'parent'), 'currentpoint');
        
        newFocalPlane = pt(1, 2);
        newFocalPlane = newFocalPlane ./ udat.gl.sliceExapandScalar; % adj for the expansion
        udat.gl.Zplane = round(newFocalPlane);
        
        % set the userdata
        set(gcf, 'userdata', udat);
        
        % update the images
        gui_updateImages;
    end
    
end



function gui_vertSliceImgClick(varargin)

    % grab the udat
    udat = get(gcf, 'userdata');

    % first, update the vertical and horizontal slices by grabing the X,Y pos
    % of the mouse click
    pt = get(udat.h.ax_sliceVert, 'currentpoint');
    udat.gl.Zplane = round(pt(1,1) ./ udat.gl.sliceExapandScalar);
    udat.gl.Yplane = round(pt(1,2));
    set(gcf, 'userdata', udat);
    gui_updateImages()

    % When the user uses the left mouse button, proceed to update the mask,
    % if they pressed the right button, do nothing. (right clicks are just
    % for navigation)
    switch get(gcf, 'SelectionType')
        case 'alt'
            return
        case 'normal'
            mask_update;
    end

end

function gui_horzSliceImgClick(varargin)

    % grab the udat
    udat = get(gcf, 'userdata');

    % first, update the vertical and horizontal slices by grabing the X,Y pos
    % of the mouse click
    pt = get(udat.h.ax_sliceHoriz, 'currentpoint');
    udat.gl.Xplane = round(pt(1,1));
    udat.gl.Zplane = round(pt(1,2) ./ udat.gl.sliceExapandScalar);
    set(gcf, 'userdata', udat);
    gui_updateImages()

    % When the user uses the left mouse button, proceed to update the mask,
    % if they pressed the right button, do nothing. (right clicks are just
    % for navigation)
    switch get(gcf, 'SelectionType')
        case 'alt'
            return
        case 'normal'
            mask_update;
    end

end

function gui_mainImgClick(varargin)

    % grab the udat
    udat = get(gcf, 'userdata');

    % first, update the vertical and horizontal slices by grabing the X,Y pos
    % of the mouse click
    pt = get(udat.h.ax_curFrame, 'currentpoint');
    udat.gl.Xplane = round(pt(1,1));
    udat.gl.Yplane = round(pt(1,2));
    set(gcf, 'userdata', udat);
    gui_updateImages()
    
    % When the user uses the left mouse button, proceed to update the mask,
    % if they pressed the right button, do nothing. (right clicks are just
    % for navigation)
    switch get(gcf, 'SelectionType')
        case 'alt'
            return
        case 'normal'
            mask_update;
    end

end

function gui_keypress(~, key)

    % grab the user data
    udat = get(gcf, 'userdata');

    % route the key press to the appropriate place
    switch lower(key.Key)
        case {'z', 'downarrow'}
            udat.gl.Zplane = min([udat.gl.Zplane+1, udat.raw.info.Nframes]);
            set(gcf, 'userdata', udat)
            gui_updateImages
        case {'a', 'uparrow'}
            udat.gl.Zplane = max([udat.gl.Zplane-1, 1]);
            set(gcf, 'userdata', udat)
            gui_updateImages
        case 'leftarrow'
            udat.gl.Xplane = max([udat.gl.Xplane-1, 1]);
            set(gcf, 'userdata', udat)
            gui_updateImages
        case 'rightarrow'
            udat.gl.Xplane = min([udat.gl.Xplane+1, udat.raw.info.Width]);
            set(gcf, 'userdata', udat)
            gui_updateImages
        case 'c'
            udat.gl.Yplane = min([udat.gl.Yplane+1, udat.raw.info.Height]);
            set(gcf, 'userdata', udat)
            gui_updateImages
        case 'd'
            udat.gl.Yplane = max([udat.gl.Yplane-1, 1]);
            set(gcf, 'userdata', udat)
            gui_updateImages
        case 'space'
            why; % egg :)
        case 'backspace'
            set(gcf, 'userdata', udat);
            mask_removeLastCell

    end
end


function mask_update()

% pull out the udat
udat = get(gcf, 'userdata');

% Make the mask if it's not defined yet (e.g., the first time you left
% click in the main window.)
if ~isfield(udat, 'mask')
    udat.mask.img = zeros(size(udat.raw.img));
    udat.mask.Ncells = 0;
    
    % Set the user data and return
    set(gcf, 'userdata', udat)
    return
end

% create a temporary mask that fits around the point selected
xmax = udat.raw.info.Width;
ymax = udat.raw.info.Height;
zmax = udat.raw.info.Nframes;
[x,y,z] = ndgrid(1:ymax, 1:xmax, 1:zmax);
x = x-udat.gl.Yplane;
y = y-udat.gl.Xplane;
z = z-udat.gl.Zplane;
ellipsoid = sqrt( (x.^2)/4 + (y.^2)/4 + z.^2 );

cellsize = 3; % default: 4
surroundsize = 7; % default: 13
cellmask = ellipsoid < cellsize;
surroundmask =  (ellipsoid < surroundsize) & ~((ellipsoid < cellsize+1.5) | udat.mask.img); % making a shell around the cell mask

% are there any existing (other cells) in the window just defined? If so,
% return without adding any new cells.
existingCells = udat.mask.img(cellmask) > 0;
if any(existingCells)
    disp('This cell has already been marked.')
    return % do nothing
end


% what's the mean of the surrounding pixels?
cellmask = cellmask(:);
surroundmask = surroundmask(:);
tmpraw = udat.raw.img(:);
mean_inside = mean(tmpraw(cellmask));
mean_surround = mean(tmpraw(surroundmask));
SNR = mean_inside/mean_surround;

% only accept SNR > 1.3
if SNR < 1.5
    fprintf('SNR was too low: %.3f\n', SNR)
    return
end

% fill up a region in the cell mask
FILLTYPE = 'region_fill';
switch FILLTYPE
    case 'simple'
        
        % update the cell counter.
        udat.mask.Ncells = udat.mask.Ncells + 1;
        
        % highlight the cell filled region
        udat.mask.img(reshape(cellmask, ymax, xmax, zmax)) = udat.mask.Ncells;
        
    case 'region_fill'
        tmpmask = mask_regionFill(udat, cellmask, mean_inside);
        tmpmask = tmpmask==1;
        tmpmask = tmpmask & ~(udat.mask.img); % no overlap with existing cells.
        
        if sum(tmpmask(:))>0
            udat.mask.Ncells = udat.mask.Ncells + 1;
            udat.mask.img(tmpmask==1) = udat.mask.Ncells;
            fprintf('Adding cell #%d\n', udat.mask.Ncells);
        else
            fprintf('Cell found but not marked, try clicking on a different part of the cell\n')
        end
        
end

% Set the user data and update the images
set(gcf, 'userdata', udat)
gui_updateImages


end

function mask_quickCheck(varargin)
    
    % grab the udat
    udat = get(gcf, 'userdata');
    set(udat.h.quickCheck, 'value', 0);
    
    % some useful info
    maxX = udat.raw.info.Width;
    maxY = udat.raw.info.Height;
    maxdac = udat.raw.info.MaxSampleValue(1);
    
    
    % project all the cells down the z dimension. Do the same with the mask
    proj_raw = max(udat.raw.img, [], 3);
    proj_mask = max(udat.mask.img, [], 3);
    
    % make an RGB version of the image, but scale certian pixels according
    % to the cell mask
    clr = [1 0.7 0.7];
    maskvect = proj_mask > 0;
    rgbimg = nan(maxY, maxX, 3);
    for a = 1:3
        tmp = proj_raw./maxdac;
        tmp = tmp(:);
        tmp(maskvect(:)) = tmp(maskvect(:)).* clr(a);
        rgbimg(:,:,a) = reshape(tmp, maxY, maxX);
    end
    
    % Display the result
    h_check = figure;
    imshow(rgbimg)
    set(h_check, 'position', [16 5 376 797], 'toolBar', 'none')
    
    
end

function mask_removeLastCell()

    % grab the udat
    udat = get(gcf, 'userdata');

    % remove all instances of the last cell, and decrement the cell counter
    lastCell = udat.mask.Ncells;
    if lastCell > 0;
        tmp = udat.mask.img;
        sz3D = size(tmp);
        tmp = tmp(:);
        idx = tmp == lastCell;
        tmp(idx) = 0;
        udat.mask.img = reshape(tmp, sz3D);
        fprintf('Removing cell #%d\n', lastCell);
        udat.mask.Ncells = udat.mask.Ncells - 1;
    end

    % set the userdata and update
    set(gcf, 'userdata', udat);
    gui_updateImages
    
end

function mask_markAsMultiple(varargin)
    
    % grab the user data
    udat = get(gcf, 'userdata');
    
    % assume that the cell(s) in question were the last ones filled, so
    % they should have the highest number in the mask. Make a new 3D mask
    % that's a binary mask with ones where the cell(s) of interest is located
    if ~isfield(udat.mask, 'Ncells') || udat.mask.Ncells == 0
        fprintf('No cells have been filled yet \n')
        set(udat.h.multipleCells, 'value', 0) 
        return
    end
    orig_mask = udat.mask.img == udat.mask.Ncells;
    
    % hack off all the things you don't care about
    linIdx = find(orig_mask(:) > 0);
    [x,y,z] = ind2sub(size(orig_mask), linIdx);
    xmin = max([1, min(x)-10]);
    xmax = min([size(orig_mask,1), max(x)+10]);
    ymin = max([1, min(y)-10]);
    ymax = min([size(orig_mask,2), max(y)+10]);
    zmin = max([1, min(z)-4]);
    zmax = min([size(orig_mask,3), max(z)+4]);
    tmp_mask = orig_mask(xmin:xmax, ymin:ymax, zmin:zmax);
    
    % alert the user that this could take some time
    fprintf('*** Entering automated segmentation routine *** \n')
    fprintf('*** This could take up to 1 minute to complete *** \n')
    
    
    % present the before image
    h_check = figure;
    set(h_check, 'position', [220 160 784 636])
    subplot(1,2,1)
    [x,y,z] = meshgrid(1:size(tmp_mask,2), 1:size(tmp_mask,1), 1:size(tmp_mask,3));
    isosurface(x,y,z,tmp_mask,0.9), 
    axis equal, title('Before')
    xlabel x, ylabel y, zlabel z
    view(3), camlight, lighting gouraud
    
    % work on the watershed analysis
    D = bwdist(~tmp_mask); % compute the distance transform...
    se = strel('ball', 4, 4);
    D = imerode(D, se);
    se = strel('ball', 7, 7);
    D = imdilate(D, se);
    D = -D;    
    D(~tmp_mask) = -Inf;
    L = watershed(D);
    
    % scan for stuff that's too small to be a cell
    objects = unique(L(:));
    objects = objects(objects>1);
    l_valid = false(numel(objects),1);
    cellSize = 125;
    for a = 1:numel(objects)
        sz = sum(L(:) == objects(a));
        l_valid(a) = (sz > cellSize);
    end
    objects = objects(l_valid);
    
    
    subplot(1,2,2)
    clrval = linspace(0.7, 0.99, numel(objects));
    for a = 1:numel(objects)
        isosurface(x,y,z,L==objects(a),clrval(a))
    end
    axis equal
    title('Segmented objects')
    xlabel x, ylabel y, zlabel z
    view(3), camlight, lighting gouraud
    drawnow
    
    % user input
    h_accept = uicontrol('style', 'togglebutton',...
            'units', 'normalized',...
            'string', 'Accept',...
            'Position', [0.5, 0.05, 0.13 0.05],...
            'Callback', {@set_returnval});
    
    h_reject = uicontrol('style', 'togglebutton',...
            'units', 'normalized',...
            'string', 'Reject',...
            'Position', [0.65, 0.05, 0.13 0.05],...
            'Callback', {@set_returnval});
        
    h_manual = uicontrol('style', 'togglebutton',...
            'units', 'normalized',...
            'string', 'Manual Entry',...
            'Position', [0.80, 0.05, 0.13 0.05],...
            'Callback', {@set_returnval});
        
    drawnow
    rotate3d on;
    
    % helper function
    function set_returnval(varargin)
        set(h_check, 'name', 'resolved');
    end
    
    % wait for the user to push one of the two buttons. This is important
    % b/c we don't want the user to make any more dicisions in the main gui
    % window until this issue is resolved.
    waitfor(h_check, 'name', 'resolved');
    
    % now figure out which button they pressed.
    assert(sum([get(h_accept, 'value'), get(h_reject, 'value'), get(h_manual, 'value')]) == 1,...
            'ERROR: Please enter only one selection... unsure how this happened... Cells were not added or deleted');
        
    if get(h_accept, 'value')
        decision = 'accept';
    elseif get(h_reject, 'value')
        decision = 'reject';
    elseif get(h_manual, 'value')
        decision = 'manual';
    else
        error('Unknown selection')
    end


    % now deal with the decision
    switch decision
        case 'accept'
            
            % provide some feedback to the user
            fprintf('Accepted: %d new cell(s) added\n', numel(objects))
            
            % delete the representation in the old mask, but don't disturb
            % the other pixels. Turning the tmp_mask back into a 'double'
            % is important because it starts out as a logical, and logical
            % arrays can only have 1 and 0, but I want each cloud of cells
            % to have different integers.
            tmp_mask = udat.mask.img(xmin:xmax, ymin:ymax, zmin:zmax);
            idx = tmp_mask == udat.mask.Ncells;
            tmp_mask(idx) = 0;
            
            % reduce the cell count back down by 1
            udat.mask.Ncells = udat.mask.Ncells-1;
            
            % update the tmp_mask (which is a subset of the orig_mask)
            for a = 1:numel(objects)
                udat.mask.Ncells = udat.mask.Ncells + 1;
                tmp_mask(L==objects(a)) = udat.mask.Ncells;
                fprintf('  Adding cell #%d\n', udat.mask.Ncells)
            end
            
            % put the tmp_mask back into the original mask
            udat.mask.img(xmin:xmax, ymin:ymax, zmin:zmax) = tmp_mask;
            
            
        case 'manual'
            
            % promp the user to specify how many cells are present
            prompt={'How many cells are present?'};
            name='Manual Entry For Multiple Cells';
            numlines=1;
            answer = inputdlg(prompt,name,numlines);
            answer = str2double(answer);
            fprintf('Manual Entry: %d new cell(s) added\n', answer)
            
            % use the original cell volume, but fill it up with multiple
            % integers (one for each "cell"). This will mean that multiple
            % cells with fill the exact same valume. This is a hack, but
            % allows the program to count multiple cells even when the
            % region filling algorithm doesn't work well AND the automated
            % routine to split cells doesn't work either.
            tmp_mask = udat.mask.img(xmin:xmax, ymin:ymax, zmin:zmax);
            sz_mask = size(tmp_mask);
            tmp_mask = tmp_mask(:); % for linear indexing
            idx = find(tmp_mask == udat.mask.Ncells);
            [~,~,z] = ind2sub(sz_mask, idx);
            
            % reset the cell counter
            udat.mask.Ncells = udat.mask.Ncells-1;
            
            % re-fill the original cell volume with the appropriate number
            % of cells
            if answer == 0
                tmp_mask(idx) = 0;
                fprintf('Manual Entry: previously selected cell has been deleted\n')
                
            else
                
                % try to maintain the footprint in the xy dims and divide
                % the volume through the z dim
                assert(numel(unique(z)) > answer, 'ERROR: you asked for more cells than can be accomodated in this volume')
                
                zPlanes = unique(z);
                zPlanes = round(linspace(zPlanes(1), zPlanes(end), answer+1));
                assert( numel(unique(zPlanes))-1 == answer, 'ERROR: sub-dividing through the z-dim has failed')
                
                
                % update cell counter, provide feedback, and adjust the
                % entries in the cell volume
                for a = 1:answer
                    udat.mask.Ncells = udat.mask.Ncells+1;
                    fprintf('  Adding cell #%d\n', udat.mask.Ncells);
                    
                    l_pix = (z >= zPlanes(a)) & (z <= zPlanes(a+1));
                    tmp_mask(idx(l_pix)) = udat.mask.Ncells;
                end
                
            end
            
            % put the tmp_mask back into the original mask
            tmp_mask = reshape(tmp_mask, sz_mask);
            udat.mask.img(xmin:xmax, ymin:ymax, zmin:zmax) = tmp_mask;
            
            
            
        case 'reject'
            
            % provide some feedback to the user
            fprintf('Rejected: cell #%d deleted from mask\n', udat.mask.Ncells)
            
            % for now, just delete that region from the mask
            idx = udat.mask.img == udat.mask.Ncells;
            udat.mask.img(idx) = 0;
            udat.mask.Ncells = udat.mask.Ncells-1;
           
    end
    
    
   % close the h_check figure, and reset the 'multiple cells' button on the
   % main gui window
   close(h_check) 
   set(udat.h.multipleCells, 'value', 0) 
    
   % apply the changes
   set(udat.h.main, 'userdata', udat);
   
   % update the gui
   gui_updateImages(udat.h.main);
   
   

end

function tmpmask = mask_regionFill(udat, cellmask, mean_inside)
    
    % find the starting location. Define this point as the point of maximal
    % brightness inside the cellmask
    tmpimg = udat.raw.img(:);
    tmpmask = cellmask(:); 
    val = max(tmpimg(tmpmask));
    l_eqmax = (tmpimg == val) & tmpmask;
    [x,y,z] = ind2sub(size(udat.raw.img), find(l_eqmax==1, 1, 'first'));

    
    % check this nhood. if mean(nhood)>critval, iterate over pixels in
    % nhood. iterate recursively over all nhoods. Define some things that
    % will be in the scope of the nested-subfunction:
    maxX = udat.raw.info.Width;
    maxY = udat.raw.info.Height;
    maxZ = udat.raw.info.Nframes;
    tmpmask = nan(size(udat.raw.img));
    
    % check the neighborhood surrounding the brightest pixel
    addToQueue = check_nhood(x,y,z, mean_inside);
    queue = addToQueue;
    
    % while loop to grow the region
    tic;
    TIMEOUT = false;
    while (size(queue,1) > 0) && ~TIMEOUT ;
               
       % try to do things via cellfun (empirically faster than not using
       % cellfun)
       x = mat2cell(queue(:,1), ones(size(queue,1),1));
       y = mat2cell(queue(:,2), ones(size(queue,1),1));
       z = mat2cell(queue(:,3), ones(size(queue,1),1));
       crit = repmat({mean_inside}, size(queue,1), 1);
       new_nhood = cellfun(@check_nhood, x, y, z, crit, 'uniformoutput', false);
       queue = vertcat(new_nhood{:});
       queue = unique(queue, 'rows');
             
       
       % time out after 10 seconds.
       TIMEOUT = toc>10;
       if TIMEOUT; disp('TIMED OUT'); end
        
    end
    
    %
    % nested sub function
    %
    %%%%%%%%%%%%%%%%%%%%
    function new_nhood = check_nhood(x,y,z, critval)
        
        % if this place has been tested, move on,
        if ~isnan(tmpmask(x,y,z));
            new_nhood = [];
            return
        end
        
        % otherwise, make a neighborhood and test it relative to the crit val.
        nhood = fullfact([3 3 3])-2;
        nhood_idx = bsxfun(@minus, [x,y,z], nhood);
        
        % make sure none of the neighboor hood indicies are out of bounds.
        toolow = nhood_idx < 1;
        toohigh = bsxfun(@minus, [maxY, maxX, maxZ], nhood_idx) < 0;
        l_oob = sum([toolow, toohigh], 2) > 0;
        nhood_idx(l_oob,:) = [];
        
        % test the neighborhood against a criterion
        ind = sub2ind(size(udat.raw.img), nhood_idx(:,1), nhood_idx(:,2), nhood_idx(:,3));
        aboveCrit = mean(udat.raw.img(ind)) > critval;
        
        % add to the queue if the center of the Nhood exceeds the crit val.
        if ~aboveCrit
            tmpmask(x,y,z) = 0;
            new_nhood = [];
        else
            tmpmask(x,y,z) = 1;
            new_nhood = nhood_idx;
            
            % eliminate points that have already been tested
            ind = sub2ind(size(udat.raw.img), nhood_idx(:,1), nhood_idx(:,2), nhood_idx(:,3));
            l_notTested = isnan(tmpmask(ind));
            new_nhood = new_nhood(l_notTested,:);
            
        end
        
        
    end
    
    
    
end


%
% TO DO
%
% 1) Autosave
%
% 3) standard size of window for slices
%
% 6) targeted delete
%
% 7) indicators on marigins to aid in navigation.




