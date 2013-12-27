function plotimg(img, info, scalebarsize, clrmap)



% add the scale bar
if info.XResolution == info.YResolution;
    
    if ~exist('scalebarsize', 'var')
        scalebarsize = 100; % in um
    end
    if ~strcmpi(info.ResolutionUnit, 'centimeter')
        error('Unknown resolution unit')
    end
    
    % convert pix/cm to pix for a specified length (in micrometers)
    nPixForBar = (info.XResolution/1e-2) .* (scalebarsize/1e6);
    nPixForBar = round(nPixForBar);
    
    % cram the error bar into the image mtx
    img(15:20, 40:(40+nPixForBar-1), :) = info.MaxSampleValue(1);

else
    warning('resolution in x and y are different')
end


% generate the plot, assuming that the image is truecolor, forcing the
% scaling to be between the min and max dac values
bits = info.BitsPerSample(1);
maxdac = 2^bits - 1;
imshow(img, [0, maxdac])

if exist('clrmap', 'var')
    colormap(clrmap)
else 
    colormap('jet')
end