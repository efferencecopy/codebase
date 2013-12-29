function out = preProcessImg(in, info, NPIX, METHOD)

% create the grayscale image
img_gray = my_rgb2gray(in, info);


% do median neighborhood averaging
if NPIX > 0
    fun = @(x) median(x(:));
    img_gray = (nlfilter(double(img_gray),[NPIX NPIX],fun));
end

switch METHOD
    case 'imadjust'
        out = double(img_gray)./255; % convert to values b/w 0 and 1
        out = imadjust(out, [0 1], [0 1], 1.13);
        out = out.*255; % bring it back to the 8 bit range
    case 'histeq'
        X = 0:255;
        geo = geopdf(X, 0.1);
        geo = geo./sum(geo); % make it integrate to 1;
        hgram = round(numel(in).* geo);
        out = histeq(uint8(img_gray), hgram);
    case 'adapthisteq'
        out = adapthisteq(uint8(img_gray),...
                          'ClipLimit', 1e-6,...
                          'Range', 'full',... % full =  8 bit for uint8
                          'Distribution', 'exponential',...
                          'Alpha', 0.35,...
                          'NBins', 60);
    case 'none'
        out = img_gray;
end


out = uint8(out);

















% % make a look up table
% maxdac = info.MaxSampleValue(1);
% xx = linspace(0, 1, maxdac+1);
% switch lower(LUT)
%     case 'weibull'
%         thresh = 0.7;
%         slope = 1.4;
%         LUT = 1-exp(-(xx./thresh).^slope);
%     case 'naka'
%         n = 4;
%         c50 = 0.7;
%         rmax = 1.25;
%         LUT = rmax.*((xx.^n)./(c50.^n + xx.^n));
%         LUT = min(LUT, 1); % don't accept values larger than one...
%         img_medFilt = round(img_medFilt ./ 40);
%     case 'linear'
%         LUT = xx;
% end
% 
% % apply the look up table
% idx = double(img_medFilt(:)) + 1; % vals b/w 1 and 256 for index to LUT
% tmp = round(LUT(idx) .* maxdac);
% 
% % reshape and return a uint8 img
% out = reshape(tmp, size(in,1), size(in, 2));
% out = uint8(out);