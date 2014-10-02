function opticaldissector(tiffstack, colorchannel)

%  open the tiff stack
%
% if 'tiffstack' is a path, open the stack. if not, open a ui control to
% navigate to the data file
%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('tiffstack', 'var') || isempty(tiffstack)
    [filename, pathname] = uigetfile('*.tif*', 'Pick a tiff stack');
    tiffstack = [pathname, filesep, filename];
end
if ~exist('colorchannel', 'var')
    colorchannel = [];
end
raw = io_unpacktiff(tiffstack, colorchannel); % import the data.



%
% setup the gui
%
%%%%%%%%%%%%%%%%%%%%%%%%%
gui_initialize(raw)





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


function gui_initialize(raw)
    
    % the top level figure
    h.main = figure;
    set(h.main, 'name', sprintf('Optical Dissector'),...
                'units', 'normalized',...
                'position', [0.2, 0, 0.6, 1])
            
    % useful constants    
    lims = [raw.info.MinSampleValue(1) raw.info.MaxSampleValue(1)];
    raw.clrmap = colormap('gray');%pmkmp(256, 'CubicL');
    
    % add an axis for the main focal plane
    gl.Zplane = 1;
    h.ax_curFrame = axes('position', [0.10, 0.15, .31, .8]);
    h.img_curFrame = imshow(raw.img(:,:,gl.Zplane), lims,...
                            'parent', h.ax_curFrame,...
                            'colormap', raw.clrmap);
    set(h.img_curFrame, 'ButtonDownFcn', {@gui_mainImgClick});
                        
    % add an axis for a slice through the vertical dimension
    h.ax_sliceVert = axes('position', [0.35, 0.70, .6, .31]);
    gl.Xplane = round(raw.info.Width/2); % default starting slice is in the middle
    gl.sliceExapandScalar = round(50./raw.info.Nframes);
    vertSlice = squeeze(raw.img(:, gl.Xplane, :));
    vertSlice = repmat(vertSlice, gl.sliceExapandScalar, 1);
    vertSlice = reshape(vertSlice, raw.info.Height, [])';
    h.img_sliceVert = imshow(vertSlice, lims,...
                            'parent', h.ax_sliceVert,...
                            'colormap', raw.clrmap);
    set(h.img_sliceVert, 'ButtonDownFcn', {@gui_getNewFocalPlane});
    
    % add an axis for a slice through the horizontal dimension
    h.ax_sliceHoriz = axes('position', [0.35, 0.60, .3, .31]);
    gl.Yplane = round(raw.info.Height/2); % default starting slice is in the middle
    horizSlice = squeeze(raw.img(gl.Yplane,:, :));
    horizSlice = repmat(horizSlice, gl.sliceExapandScalar, 1);
    horizSlice = reshape(horizSlice, raw.info.Width, [])';
    h.img_sliceHoriz = imshow(horizSlice, lims,...
                            'parent', h.ax_sliceHoriz,...
                            'colormap', raw.clrmap);
    set(h.img_sliceHoriz, 'ButtonDownFcn', {@gui_getNewFocalPlane});
    
                        
    % package the relevant info into the user data field
    udat.gl = gl;
    udat.h = h;
    udat.raw = raw;
    set(h.main, 'userdata', udat);
    
end

function gui_updateImages()
    
    % grab stuff from the user data field
    udat = get(gcf, 'userdata');
    
    % new main image focal plane
    set(udat.h.img_curFrame, 'CData', udat.raw.img(:,:, udat.gl.Zplane))
     
    % new vertical slice
    vertSlice = squeeze(udat.raw.img(:, udat.gl.Xplane, :));
    vertSlice = repmat(vertSlice, udat.gl.sliceExapandScalar, 1);
    vertSlice = reshape(vertSlice, udat.raw.info.Height, [])';
    set(udat.h.img_sliceVert, 'CData', vertSlice)
    
    % new horizontal slice
    horizSlice = squeeze(udat.raw.img(udat.gl.Yplane,:, :));
    horizSlice = repmat(horizSlice, udat.gl.sliceExapandScalar, 1);
    horizSlice = reshape(horizSlice, udat.raw.info.Width, [])';
    set(udat.h.img_sliceHoriz, 'CData', horizSlice)

end

function gui_getNewFocalPlane(h_slice, ~)
    
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


function mask_update()

% pull out the udat
udat = get(gcf, 'userdata');

% Make the mask if it's not defined yet (e.g., the first time you left
% click in the main window.)
if ~isfield(udat, 'mask')
    udat.mask.img = zeros(size(udat.raw.img));
    udat.mask.Ncells = 0;
end

% create a temporary mask that fits around the point selected
xmax = udat.raw.info.Width;
ymax = udat.raw.info.Height;
zmax = udat.raw.info.Nframes;
[x,y,z] = ndgrid(1:ymax, 1:xmax, 1:zmax);
x = x-udat.gl.Yplane;
y = y-udat.gl.Xplane;
z = z-udat.gl.Zplane;
ellipsoid = sqrt( (x.^2)/9 + (y.^2)/9 + z.^2 );

cellmask = ellipsoid < 2;
surroundmask =  (ellipsoid < 5) & ~((ellipsoid < 2.5) | udat.mask.img); % making a shell around the cell mask

% what's the mean of the surrounding pixels?
cellmask = cellmask(:);
surroundmask = surroundmask(:);
tmpraw = udat.raw.img(:);
mean_inside = mean(tmpraw(cellmask));
mean_surround = mean(tmpraw(surroundmask));
SNR = mean_inside/mean_surround;

% only accept SNR > 1.4
if SNR < 1.4
    fprintf('SNR was too low: %.3f\n', SNR)
    return
end

% add some some info to the cell mask and update the cell counter.


end






