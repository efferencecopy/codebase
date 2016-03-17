function raw_img = loadTIFF(fpath, info, frameNums)
%
% EXAMPLE:  raw_img = loadTIFF(fpath, info, [frameNums])
%
% fpath -> a fully qualified path to the data file
% info  -> an "info" structure from imfinfo.m
% frameNums  -> optional, a 2 element vector describing the first and last frame
%
% raw_img -> the output image, as a uint of the same bitdepth as the original image
%
% C.Hass 3/2016


warnID = 'MATLAB:imagesci:tiffmexutils:libtiffWarning';
warning('off', warnID);

% pull out the dimensions of the image
width_pix = info(1).Width;
height_pix = info(1).Height;
Nframes = length(info);

% which frames should be extracted. Default is all of them...
if exist('frameNums', 'var')
    startFrame = frameNums(1);
    endFrame = frameNums(2);
    Nframes = endFrame - startFrame + 1;
else 
    startFrame = 1;
    endFrame = Nframes;
end

% what's the bit depth
switch info(1).BitDepth;
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
raw_img = zeros(height_pix, width_pix, Nframes, bitdepth);

tifflink = Tiff(fpath, 'r');
for i_frame = startFrame:endFrame
   tifflink.setDirectory(i_frame);
   raw_img(:,:,i_frame)=tifflink.read;
end
tifflink.close();



