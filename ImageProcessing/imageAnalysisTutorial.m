% An attempt to import TIFF files and process the images automatically.

% clear out the workspace
fin

% load a file
fileName = 'CH_1126_A_p6_s1_10x_red';
idx = regexp(fileName, '_p\d+', 'start');
datPath = [GL_DATPATH, fileName(1:idx-1)];
imgPath = findfile(fileName, datPath, '.tif');

img = imread(imgPath);
info = imfinfo(imgPath);

% define some things
SCALEBAR = 1000;    % in um
NPIX = 0;           % in pix (for averaging)
LUT = 'weibull';    % could be 'weibull', or


%% plot just to confirm

figure
plotimg(img, info, SCALEBAR)

labels = {'red channel', 'green channel', 'blue channel'};
figure
colormap(jet(256))
for a = 1:3
    subplot(2,3,a)
    plotimg(img(:,:,a), info, SCALEBAR')
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
    subplot(2,1,1)
    plotimg(img_gray, info, SCALEBAR, 'gray');
    colorbar
    title(num2str(coeff{a}))
    
    subplot(2,1,2)
    hist(double(img_gray(:)), 256)
    axis tight
    xlim([6 256])
end

figure
img_gray = rgb2gray(img);
plotimg(img_gray, info, SCALEBAR, 'gray');
title('rgb2gray')

%% Now try to do pix averaging, followed by contrast enhancement

% create the grayscale image, and then do median neighborhood averaging
% with different contrast enhancement techniques.

figure
subplot(1,3,1)
tmp = preProcessImg(img, info, NPIX, 'imadjust');
plotimg(tmp, info, SCALEBAR, 'gray');

subplot(1,3,2)
tmp = preProcessImg(img, info, NPIX, 'histeq');
plotimg(tmp, info, SCALEBAR, 'gray');

subplot(1,3,3)
tmp = preProcessImg(img, info, NPIX, 'adapthisteq');
plotimg(tmp, info, SCALEBAR, 'gray');





%% Trying to make a merged image.

close all
fileName = 'CH_1126_A_p8_s2_10x';


% deal with the green channel
idx = regexp(fileName, '_p\d+', 'start');
datPath = [GL_DATPATH, fileName(1:idx-1)];
imgPath = findfile([fileName, '_green'], datPath, '.tif');
img_green = imread(imgPath);
info_green = imfinfo(imgPath);
ch_green = preProcessImg(img_green, info_green, 0, 'imadjust');

% deal with red channel
imgPath = findfile([fileName, '_red'], datPath, '.tif');
img_red = imread(imgPath);
info_red = imfinfo(imgPath);
ch_red = preProcessImg(img_red, info_red, 0, 'imadjust');


% merge the two channels and create a non-sense blue channel
img_merge = cat(3, ch_red, ch_green, zeros(size(ch_red)));
figure
plotimg(img_merge, info);

% compare to simple way:

img_merge_simple = cat(3, my_rgb2gray(img_red, info_red), my_rgb2gray(img_green, info_green), zeros(size(ch_red)));
figure
plotimg(img_merge_simple, info);

%% VIEW A STACK OF IMAGES


params.mouse = 'EMX_1';
params.objective = '2x';
params.contrastMethod = 'none';
params.npix = 0;

[stack, imgTags] = colorMerge(params);
stackViewer(stack, imgTags);


