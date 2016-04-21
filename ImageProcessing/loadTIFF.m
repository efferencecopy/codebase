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


warnID = 'MATLAB:imagesci:tiffmexutils:libtiffWarning';
warning('off', warnID);

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

% preallocate space
raw_img = zeros(info.height_pix, info.width_pix, Nframes, bitdepth);

tifflink = Tiff(fpath, 'r');
for i_frame = startFrame:endFrame
   tifflink.setDirectory(i_frame);
   raw_img(:,:,i_frame)=tifflink.read;
end
tifflink.close();


