function [dFoF, Fo] = dfof_from_tiffstack(img_raw, NsampsPerTrial)
%
% EXAMPLE: dFoF = dfof_from_tiffstack(img_raw, frameRate, window_size_sec)

% add eps to all the values of the raw image. This avoids 0./0 errors,
% which happen for black pixels.
zeroFrames = sum(sum(img_raw,1),2) == 0;
fprintf('  found %d blank frames.', sum(zeroFrames));
fprintf('  found %d blank pixels.  ', sum(img_raw(:)==0));
img_raw = img_raw+eps;

% setup the box car filter coeffs
window_size_samps = NsampsPerTrial .* 3; % Estimate Fo as the average F over the past three trials
B_box_Fo = (1/window_size_samps).*ones(1,window_size_samps);



% make sure that the background period doesn't have any frames that are all
% blank, which will cause the dfof to be wonky. this could hapen if the
% user turns on the light after hitting the go button. This fix assumes
% that the problem is isolated to the first few seconds, if not, then there
% will be additional problems later.
bkgnd = img_raw(:,:,1:window_size_samps); %#ok<*BDSCI>
normal_stderr = std(img_raw(:,:,window_size_samps+1:2*window_size_samps), [], 3);
normal_raw = mean(img_raw(:,:,window_size_samps+1:2*window_size_samps), 3);
critval = 3.* normal_stderr;
pix_oob = abs(bsxfun(@minus, bkgnd, normal_raw)) > repmat(critval, [1,1,window_size_samps]);
Percent_pix_oob = sum(sum(pix_oob,1),2) ./ numel(img_raw(:,:,1));
l_oob_frames = Percent_pix_oob > 0.50; % if half the pix are black, assume the frame is a problem

inds = find(l_oob_frames);
for i_oob = 1:sum(l_oob_frames)
    img_raw(:,:,inds(i_oob)) = normal_raw;
end


% front-pad the image stack so that the filter kernal comes to steady state
% before it hits the data
bkgnd = mean(img_raw(:,:,1:window_size_samps),3);
pad_length = 5;
pad = repmat(bkgnd, [1,1,window_size_samps+pad_length]);
tmp_img = cat(3, pad, img_raw);

% back-pad the imagestack so that the filter can be non-causal (by shifting
% the output)
half_window = ceil( window_size_samps ./2 );
endvals = mean(img_raw(:,:,end-(half_window-1):end),3);
pad = repmat(endvals, [1, 1, half_window]);
tmp_img = cat(3, tmp_img, pad);

% run the filter for Fo
Fo = filter(B_box_Fo, 1, tmp_img, [], 3);
startidx = pad_length + window_size_samps + half_window + 1;
Fo = Fo(:,:,startidx:end);

% now just do the math
dFoF =  (img_raw - Fo) ./ Fo;

% crazy idea: subtract the mean across all pixels to eliminate image wide
% noise that has nothing to do with IAF or hemodynamics. Make sure that the
% thing you subtract off has the same sigma as each pixel time series.
xbar = mean(mean(dFoF,1),2);
sigma = std(dFoF,[],3);
scaleFactor = sigma ./ std(xbar(:)); % a scale factor to equate sigma on a pix by pix basis
rsub = bsxfun(@times, xbar, scaleFactor);
assert(all(all((sigma - std(rsub,[],3))<1e-10)), 'ERROR: sigmas are not the same')
dFoF = dFoF - rsub;


% pix = 150;
% tmp_raw = permute(img_raw(pix:pix+10, pix:pix+10, :), [3,1,2]);
% tmp_raw = reshape(tmp_raw, size(tmp_raw, 1), []);
% tmp_raw = mean(tmp_raw, 2);
% 
% tmp_fo = permute(Fo(pix:pix+10, pix:pix+10, :), [3,1,2]);
% tmp_fo = reshape(tmp_fo, size(tmp_fo, 1), []);
% tmp_fo = mean(tmp_fo, 2);
% 
% figure
% hold on,
% plot(tmp_raw, 'k')
% plot(tmp_fo, 'r', 'linewidth', 2)
