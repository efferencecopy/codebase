function [dFoF, Fo] = dfof_from_tiffstack(img_raw, frameRate, window_size_sec)
%
% EXAMPLE: dFoF = dfof_from_tiffstack(img_raw, frameRate, window_size_sec)

% add eps to all the values of the raw image. This avoids 0./0 errors,
% which happen for black pixels.
img_raw = img_raw + eps;

% setup t box car filter coeffs
window_size_samps = round(frameRate .* window_size_sec);
B_box = (1/window_size_samps).*ones(1,window_size_samps);

% front-pad the image stack so that the filter kernal comes to steady state
% before it hits the data
assert(size(img_raw,3)>=5, 'ERROR: not enough frames for baseline')
pad = repmat(mean(img_raw(:,:,1:5),3), [1,1,window_size_samps+5]);
tmp_img = cat(3, pad, img_raw);

% run the filter
Fo = filter(B_box, 1, tmp_img, [], 3);

% remove the front padding
Fo(:,:,1:window_size_samps+5) = [];

% now just do the math
dFoF =  (img_raw - Fo) ./ Fo;

