function [dFoF, Fo] = dfof_from_tiffstack(img_raw, NsampsON, NsampsOFF)
%
% EXAMPLE: dFoF = dfof_from_tiffstack(img_raw, frameRate, window_size_sec)


% throw an error if there are frames where all the pixels are black. This
% is potentially a problem for the method I'm using to calculate the Fo.
l_zero = sum(sum(img_raw, 1), 2) == 0;
assert(sum(l_zero) == 0, 'ERROR: found frames that were blank') % consider using iafExtractRuns_offperiod


% add eps to all the values of the raw image. This avoids 0./0 errors,
% which happen for black pixels.
img_raw = img_raw+eps;


% setup the box car filter coeffs
NsampsPerTrial = NsampsON+NsampsOFF;
window_size_samps = NsampsPerTrial .* 3; % Estimate Fo as the average F over three trials (1.5 trials pre and post)
B_box_Fo = (1/window_size_samps).*ones(1,window_size_samps);


% front-pad the image stack so that the filter kernal comes to steady state
% before it hits the data
bkgnd = mean(img_raw(:,:,1:window_size_samps),3);
pad_length = 5;
front_pad = repmat(bkgnd, [1,1,window_size_samps+pad_length]);

% back-pad the imagestack so that the filter can be non-causal (by shifting
% the output)
half_window = ceil( window_size_samps ./2 );
endvals = mean(img_raw(:,:,end-(half_window-1):end),3);
back_pad = repmat(endvals, [1, 1, half_window]);
tmp_img = cat(3, front_pad, img_raw, back_pad);


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



%
%  FIGURES FOR DEBUGGING
%
%

pix = 150;
tmp_raw = permute(img_raw(pix:pix+10, pix:pix+10, :), [3,1,2]);
tmp_raw = reshape(tmp_raw, size(tmp_raw, 1), []);
tmp_raw = mean(tmp_raw, 2);

tmp_fo = permute(Fo(pix:pix+10, pix:pix+10, :), [3,1,2]);
tmp_fo = reshape(tmp_fo, size(tmp_fo, 1), []);
tmp_fo = mean(tmp_fo, 2);

tmp_dfof = permute(dFoF(pix:pix+10, pix:pix+10, :), [3,1,2]);
tmp_dfof = reshape(tmp_dfof, size(tmp_dfof, 1), []);
tmp_dfof = mean(tmp_dfof, 2);

figure
hold on,
plot(tmp_raw, 'b')
plot(tmp_fo, 'r', 'linewidth', 2)
plotyy([nan], [nan], 1:numel(tmp_fo), tmp_dfof)
drawnow
