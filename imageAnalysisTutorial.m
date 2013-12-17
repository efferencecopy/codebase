% An attempt to import TIFF files and process the images automatically.

% clear out the workspace
fin

% load a file
fileName = 'CH_1126_A_p5_s4_2x_red';
imgPath = findfile(fileName, IMGPATH, '.tif');

img = imread(imgPath);
info = imfinfo(imgPath);

% define some things
SCALEBAR = 1000;    % in um
NPIX = 5;          % in pix (for averaging)


%% plot just to confirm

figure
plotimg(img, info, SCALEBAR)

labels = {'red channel', 'green channel', 'blue channel'};
figure
colormap(jet(256))
for a = 1:3
    subplot(2,3,a)
    plotimg(img(:,:,a), info, SCALEBAR)
    title(labels{a})
    colorbar
end

for a = 4:6
    subplot(2,3,a)
    tmp = img(:,:,a-3);
    hist(double(tmp(:)), 256)
    axis tight
    xlim([get(gca, 'xlim') - [6 0]])
end

%% Plot the grayscale version

coeff = {[1,1,1], [.5 .7 .2], [.9 .3 0], [.7 .7 1], [0.2989    0.5870    0.1140]}

for a = 1:numel(coeff)
    img_gray = my_rgb2gray(img, info, coeff{a});
    
    figure
    plotimg(img_gray, info, SCALEBAR, 'gray')
    title(num2str(coeff{a}))
end

figure
img_gray = rgb2gray(img);
plotimg(img_gray, info, SCALEBAR, 'gray')
title('rgb2gray')

%% Now try to do pix averaging, followed by sigmoidal LUT

% create the grayscale image, and then do median neighborhood averaging.
fun = @(x) median(x(:));
img_gray = rgb2gray(img);
img_medFilt = uint8(nlfilter(double(img_gray),[NPIX NPIX],fun));

% apply a lookup table to the smoothed image
maxdac = info.MaxSampleValue(1);
thresh = 0.6;
slope = 1.5;
xx = linspace(0, 1, maxdac+1);
LUT = 1-exp(-(xx./thresh).^slope);
idx = double(img_medFilt(:)) + 1; % vals b/w 1 and 256 for index to LUT
tmp = round(LUT(idx) .* maxdac);
img_filtAndThresh = reshape(tmp, size(img,1), size(img, 2));

figure
subplot(1,3,1)
plotimg(img_gray, info, SCALEBAR, 'gray')

subplot(1,3,2)
plotimg(img_medFilt, info, SCALEBAR, 'gray')

subplot(1,3,3)
plotimg(img_filtAndThresh, info, SCALEBAR, 'gray')





%% Trying to make a merged image.

fileName = 'CH_1126_A_p5_s4_2x';


% deal with the green channel
imgPath = findfile([fileName, '_green'], IMGPATH, '.tif');
img_green = imread(imgPath);
info = imfinfo(imgPath);
ch_green = preProcessImg(img_green, info, NPIX);

% deal with red channel
imgPath = findfile([fileName, '_red'], IMGPATH, '.tif');
img_red = imread(imgPath);
info = imfinfo(imgPath);
ch_red = preProcessImg(img_red, info, NPIX);


% merge the two channels and create a non-sense blue channel
img_merge = cat(3, ch_red, ch_green, zeros(size(ch_red)));
figure
plotimg(img_merge, info)

% compare to simple way:

img_merge_simple = cat(3, rgb2gray(img_red), rgb2gray(img_green), zeros(size(ch_red)));
figure
plotimg(img_merge_simple, info)







