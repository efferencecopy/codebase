% An attempt to import TIFF files and process the images automatically.

% clear out the workspace
fin

% load a file
fileName = 'CH_1126_A_p6_s1_10x_red';
imgPath = findfile(fileName, GL_IMGPATH, '.tif');

img = imread(imgPath);
info = imfinfo(imgPath);

% define some things
SCALEBAR = 1000;    % in um
NPIX = 5;           % in pix (for averaging)
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
    plotimg(img_gray, info, SCALEBAR, 'gray')
    colorbar
    title(num2str(coeff{a}))
    
    subplot(2,1,2)
    hist(double(img_gray(:)), 256)
    axis tight
    xlim([6 256])
end

figure
img_gray = rgb2gray(img);
plotimg(img_gray, info, SCALEBAR, 'gray')
title('rgb2gray')

%% Now try to do pix averaging, followed by sigmoidal LUT

% create the grayscale image, and then do median neighborhood averaging
% with different contrast enhancement techniques.

figure
subplot(1,3,1)
tmp = preProcessImg(img, info, NPIX, 'imadjust');
plotimg(tmp, info, SCALEBAR, 'gray')

subplot(1,3,2)
tmp = preProcessImg(img, info, NPIX, 'histeq');
plotimg(tmp, info, SCALEBAR, 'gray')

subplot(1,3,3)
tmp = preProcessImg(img, info, NPIX, 'adapthisteq');
plotimg(tmp, info, SCALEBAR, 'gray')





%% Trying to make a merged image.

close all
fileName = 'CH_1126_A_p8_s2_10x';


% deal with the green channel
imgPath = findfile([fileName, '_green'], GL_IMGPATH, '.tif');
img_green = imread(imgPath);
info_green = imfinfo(imgPath);
ch_green = preProcessImg(img_green, info_green, 0, 'imadjust');

% deal with red channel
imgPath = findfile([fileName, '_red'], GL_IMGPATH, '.tif');
img_red = imread(imgPath);
info_red = imfinfo(imgPath);
ch_red = preProcessImg(img_red, info_red, 0, 'imadjust');


% merge the two channels and create a non-sense blue channel
img_merge = cat(3, ch_red, ch_green, zeros(size(ch_red)));
figure
plotimg(img_merge, info)

% compare to simple way:

img_merge_simple = cat(3, my_rgb2gray(img_red, info_red), my_rgb2gray(img_green, info_green), zeros(size(ch_red)));
figure
plotimg(img_merge_simple, info)



%% Automate makeing a stack

fin

mouse = 'CH_1126_A';
objective = '2x';

cd([GL_IMGPATH, filesep, mouse, filesep, 'Raw Images']);

% grab the names in the directory
d = dir;

% initialize the structure of images
[img.green, img.red, info.green, info.red] = deal({});


for a = 1:numel(d);
    
    % display the progress
    if ~rem(a,5)
        fprintf('%d more images to unpack\n', numel(d)-(a-1));
    end
    
    if ~any(strcmp(d(a).name, {'.', '..'})) % skip the hidden files
        
        % make sure the objective used is correct
        sliceObjective = regexpi(d(a).name , '_\d+x', 'match');
        sliceObjective = sliceObjective{1}(2:end);
        if strcmpi(sliceObjective, objective)
            
            
            % id the plate number
            plate = regexpi(d(a).name , '_p\d+', 'match');
            plate = str2double(plate{1}(3:end));
            
            % id the slice number
            slice = regexpi(d(a).name , '_s\d+', 'match');
            slice = str2double(slice{1}(3:end));
            
            
            % red? or green?
            if regexp(d(a).name, '_green'); color = 'green'; end
            if regexp(d(a).name, '_red');   color = 'red'; end
            
            
            % unpack the images
            img_tmp = imread(d(a).name);
            info_tmp = imfinfo(d(a).name);
            img_tmp = preProcessImg(img_tmp, info_tmp, 0, 'adapthisteq');
            
            
            % put the image in a structure according to it's position in
            % the brain, and the color channel
            img = setfield(img, color, {plate, slice},  {img_tmp});
            info = setfield(info, color, {plate, slice}, {info_tmp});
            
        end
    end
end


% combine the images into a merged truecolor RGB. Arrange them in a stack.
idx = 1;
for p = 1:size(img.red, 1)
    for sl = 1:size(img.red, 2)
        
        % find the appropriate images
        ch_green = img.green{p, sl};
        ch_red = img.red{p, sl};
        
        % do some basic error checking
        if isempty(ch_green) && ~isempty(ch_green);
            error('One channel is defined but not the other')
        end
        if isempty([ch_green, ch_red]); continue; end % both channels are undefined... no big deal so no error
        
        % add images to the stack, make up a fake blue channel
        stack{idx} = cat(3, ch_red, ch_green, zeros(size(ch_red)));

        idx = idx + 1;
        
    end
end


imgAlign(stack, info_tmp);
    






