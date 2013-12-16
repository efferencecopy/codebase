% An attempt to import TIFF files and process the images automatically.

% clear out the workspace
fin

% load a file
imgPath = '~/Crash/Imaging/CH_1126_A/Raw Images/CH_1126_A_p4_s3_2x_green.tif';
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
for a = 1:3
    subplot(2,3,a)
    plotimg(img(:,:,a), info, SCALEBAR)
    title(labels{a})
end

for a = 4:6
    subplot(2,3,a)
    tmp = img(:,:,a-3);
    hist(double(tmp(:)), 256)
    axis tight
    xlim([get(gca, 'xlim') - [6 0]])
end

%% Plot the grayscale version

coeff = {[1,1,1], [.5 .7 .2], [.8 .5 1], [.7 .7 1], [0.2989    0.5870    0.1140]}

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

fun = @(x) median(x(:));
img_medFilt = uint8(zeros(size(img)));
for a = 1:3;
    img_medFilt(:,:,a) = (nlfilter(double(img(:,:,a)),[NPIX NPIX],fun));
end

% plot the results
figure
plotimg(img_medFilt, info, SCALEBAR)
title('Median Filtered')

figure
plotimg(rgb2gray(img_medFilt), info, SCALEBAR, 'gray');

% now for grayscale
img_gray = rgb2gray(img);
img_medFilt = uint8(nlfilter(double(img_gray),[NPIX NPIX],fun));
figure
plotimg(img_medFilt, info, SCALEBAR, 'gray')



















