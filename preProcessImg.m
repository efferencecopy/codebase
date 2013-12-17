function out = preProcessImg(in, info, NPIX)

% create the grayscale image, and then do median neighborhood averaging.
if any(regexp(info.Filename, '_red'))
    coeff = [.9, .3, 0];
elseif any(regexp(info.Filename, '_green'))
    coeff = [0.2989    0.5870    0.1140];
end

fun = @(x) median(x(:));
img_gray = my_rgb2gray(in, info, coeff);
img_medFilt = (nlfilter(double(img_gray),[NPIX NPIX],fun));

% make a look up table
maxdac = info.MaxSampleValue(1);
thresh = 0.7
slope = 1.4;
xx = linspace(0, 1, maxdac+1);
LUT = 1-exp(-(xx./thresh).^slope);

% apply the look up table
idx = double(img_medFilt(:)) + 1; % vals b/w 1 and 256 for index to LUT
tmp = round(LUT(idx) .* maxdac);

% reshape and return a uint8 img
out = reshape(tmp, size(in,1), size(in, 2));
out = uint8(out);