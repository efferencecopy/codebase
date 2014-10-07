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
    
    % save the data interactively
    uisave('cellFillData', ['cellFillData_', cellFillData.raw.info.shortName, '.mat'])
    
end



function gui_initialize(raw)
    
    % the top level figure
    h.main = figure;
    set(h.main, 'name', sprintf('Optical Dissector'),...
                'units', 'normalized',...
                'position', [0.2, 0, 0.6, 1],...
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
    h.ax_sliceVert = axes('position', [.27, 0.02, .2, .96]);
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
    h.ax_sliceHoriz = axes('position', [0.45, 0.50, .3, .31]);
    gl.Yplane = round(raw.info.Height/2); % default starting slice is in the middle
    horizSlice = squeeze(raw.img(gl.Yplane,:, :));
    horizSlice = repmat(horizSlice, gl.sliceExapandScalar, 1);
    horizSlice = reshape(horizSlice, raw.info.Width, [])';
    h.img_sliceHoriz = imshow(horizSlice, lims,...
                            'parent', h.ax_sliceHoriz,...
                            'colormap', raw.clrmap);
    set(h.img_sliceHoriz, 'ButtonDownFcn', {@gui_horzSliceImgClick});
    
    % add a button to save the data and dump the result onto the command
    % window
    h.export = uicontrol('style', 'togglebutton',...
                         'units', 'normalized',...
                         'string', 'Export Data',...
                         'Position', [0.72, 0.265, 0.20, 0.05],...
                         'Callback', {@io_exportData});
                        
    % package the relevant info into the user data field
    udat.gl = gl;
    udat.h = h;
    udat.raw = raw;
    set(h.main, 'userdata', udat);
    
end

function gui_updateImages()
    
    % grab stuff from the user data field
    udat = get(gcf, 'userdata');
    Xidx = udat.gl.Xplane;
    Yidx = udat.gl.Yplane;
    Zidx = udat.gl.Zplane;
    maxX = udat.raw.info.Width;
    maxY = udat.raw.info.Height;
    maxdac = 255;
    
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

cellsize = 4;
surroundsize = 13;
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

% only accept SNR > 1.4
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
% 4) incrase size of main window
%
% 5) select cells in "slice" windows
%
% 6) targeted delete
%
% 7) indicators on marigins to aid in navigation.
%
% 8) slice selection with arrow keys?




