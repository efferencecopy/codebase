function hand = plotimg(img, tags, scalebarsize, clrmap)



% add the scale bar
if tags.XResolution == tags.YResolution;
    
    if ~exist('scalebarsize', 'var')
        scalebarsize = 100; % in um
    end
    if ~strcmpi(tags.ResolutionUnit, 'centimeter')
        error('Unknown resolution unit')
    end
    
    % convert pix/cm to pix for a specified length (in micrometers)
    nPixForBar = (tags.XResolution/1e-2) .* (scalebarsize/1e6);
    nPixForBar = round(nPixForBar);
    
    % cram the error bar into the image mtx
    img(15:20, 40:(40+nPixForBar-1), :) = tags.MaxSampleValue(1);

else
    warning('resolution in x and y are different')
end


% generate the plot, assuming that the image is truecolor, forcing the
% scaling to be between the min and max dac values
bits = tags.BitsPerSample(1);
maxdac = 2^bits - 1;
hand = imshow(img, [0, maxdac]);

if exist('clrmap', 'var')
    colormap(clrmap)
else 
    colormap('jet')
end