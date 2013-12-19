function out = preProcessImg(in, info, NPIX, LUT)

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
xx = linspace(0, 1, maxdac+1);
switch lower(LUT)
    case 'weibull'
        thresh = 0.7
        slope = 1.4;
        LUT = 1-exp(-(xx./thresh).^slope);
    case 'naka'
        n = 4;
        c50 = 0.7;
        rmax = 1.25;
        LUT = rmax.*((xx.^n)./(c50.^n + xx.^n));
        LUT = min(LUT, 1); % don't accept values larger than one...
        img_medFilt = round(img_medFilt ./ 40);
    case 'linear'
        LUT = xx;
end

% apply the look up table
idx = double(img_medFilt(:)) + 1; % vals b/w 1 and 256 for index to LUT
tmp = round(LUT(idx) .* maxdac);

% reshape and return a uint8 img
out = reshape(tmp, size(in,1), size(in, 2));
out = uint8(out);