function raw_img = loadTIFF(fpath, info, frameNums)
%
% EXAMPLE:  raw_img = loadTIFF(fpath, info, [frameNums])
%
% fpath -> a fully qualified path to the data file
% info  -> an "info" structure from loadTIFF_info
% frameNums  -> optional, a 2 element vector describing the first and last frame
%
% raw_img -> the output image, as a uint of the same bitdepth as the original image
%
% C.Hass 3/2016


% if the total number of frames wasn't specified, or if the user didn't
% supply an info structure, then load one in. This is slow, so it's
% preferable to use loadTIFF_info and then to get the Nframes some other
% way
if ~exist('info', 'var') || ~isfield(info, 'Nframes')
    tmpinfo = imfinfo(fpath);
    info.height_pix = tmpinfo(1).Height;
    info.width_pix = tmpinfo(1).Width;
    info.bitdepth = tmpinfo(1).BitDepth;
    info.Nframes = numel(tmpinfo);
end

% which frames should be extracted. Default is all of them...
if exist('frameNums', 'var')
    startFrame = frameNums(1);
    endFrame = frameNums(2);
    Nframes = endFrame - startFrame + 1;
else 
    startFrame = 1;
    Nframes = info.Nframes;
    endFrame = Nframes;
end

% what's the bit depth
switch info.bitdepth;
    case 8
        bitdepth = 'uint8';
    case 16
        bitdepth = 'uint16';
    case 32
        bitdepth = 'uint32';
    otherwise
        error('unknown bitdepth')
end     




% open the file
tifflink.FileID = tifflib('open',fpath,'r');
tifflink.FileName = fpath;
tifflink.Mode = 'r';

% stuff to make sure the tifflib code runs smoothly
assert(getTag(tifflink, 'Photometric') == 1, 'ERROR: min is not set to equal black')
assert(getTag(tifflink, 'PlanarConfiguration') == Tiff.PlanarConfiguration.Chunky, 'ERROR: planar config must be set to Chunky')
assert(getTag(tifflink, 'SamplesPerPixel') == 1, 'ERROR samps per pix must equal 1')
assert(~tifflib('isTiled',tifflink.FileID), 'ERROR: tiled images not yet supported')
imageDepth = getTag(tifflink, 'ImageDepth');
isTile = tifflib('isTiled',tifflink.FileID);
if (imageDepth > 1) && ~isTile
    warning(message('MATLAB:imagesci:Tiff:strippedImageDepthRead', imageDepth));
end

h = getTag(tifflink, 'ImageLength');
rps = getTag(tifflink, 'RowsPerStrip');
rps = min(rps,h);

% preallocate space
raw_img = zeros(info.height_pix, info.width_pix, Nframes, bitdepth);


% pull out the frames
for i_frame = startFrame:endFrame
    
    tifflib('setDirectory',tifflink.FileID,i_frame-1);
    
    % Go through each strip of data.
    for r = 1:rps:h
        row_inds = r:min(h,r+rps-1);
        stripNum = tifflib('computeStrip',tifflink.FileID,r-1);
        raw_img(row_inds,:, i_frame) = tifflib('readEncodedStrip',tifflink.FileID,stripNum-1);
    end
end

% close the file
tifflib('close',tifflink.FileID);
tifflink.FileID = uint64(0);


end


function tagValue = getTag(obj, tagId)
    assert(ischar(tagId))
    tagValue = tifflib('getField',obj.FileID,Tiff.TagID.(tagId));
end

