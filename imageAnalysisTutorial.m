% An attempt to import TIFF files and process the images automatically.

clear all; close all; clc

% load a file
imgPath = '~/Crash/Imaging/CH_1126_A/Raw Images/CH_1126_A_p4_s3_2x_red.tif';
img = imread(imgPath);


%% plot just to confirm

figure
imagesc(img)

figure
for a = 1:3
    subplot(2,3,a)
    imagesc(img(:,:,a))
end

for a = 4:6
    subplot(2,3,a)
    tmp = img(:,:,a-3);
    hist(double(tmp(:)), 256)
    axis tight
    xlim([get(gca, 'xlim') - [6 0]])
end