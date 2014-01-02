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



%% Automate makeing a stack

fin

mouse = 'CH_112613_A';
objective = '2x';
contrastMethod = 'none';
NPIX = 0;

% cd to where the images are
cd([GL_DATPATH, filesep, mouse, filesep, 'Histology', filesep, 'Raw Images']);

% grab the names in the directory
d = dir;

% load in (and update) the mouseDB
mdb = initMouseDB;

% initialize the structure of images
[img.green, img.red, img.blue] = deal({});


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
            if regexp(d(a).name, '_red');   color = 'red'; end
            if regexp(d(a).name, '_green'); color = 'green'; end
            if regexp(d(a).name, '_blue');   color = 'blue'; end
            
            
            % unpack the images
            img_tmp = imread(d(a).name);
            info_tmp = imfinfo(d(a).name);
            img_tmp = preProcessImg(img_tmp, info_tmp, NPIX, contrastMethod);
            
            
            % put the image in a structure according to it's position in
            % the brain, and the color channel
            img = setfield(img, color, {plate, slice},  {img_tmp});
            
        end
    end
end




% combine the images into a merged truecolor RGB. Arrange them in a stack.
idx = 1;
numPlates = max([size(img.red,1), size(img.green,1), size(img.blue,1)]);
numSlices = max([size(img.red,2), size(img.green,2), size(img.blue,2)]);

for p = 1:numPlates
    for sl = 1:numSlices
        
        % find the appropriate images
        if ((size(img.red,1)>=p) && (size(img.red,2)>=sl)) && ~isempty(img.red{p,sl})
            ch_red = img.red{p, sl};
        else
            ch_red = zeros(size(img_tmp));
        end
        
        if (size(img.green,1)>=p) && (size(img.green,2)>=sl) && ~isempty(img.green{p,sl})
            ch_green = img.green{p, sl};
        else
            ch_green = zeros(size(img_tmp));
        end
        
        if (size(img.blue,1)>=p) && (size(img.blue,2)>=sl) && ~isempty(img.blue{p,sl})
            ch_blue = img.blue{p, sl};
        else
            ch_blue = zeros(size(img_tmp));
        end
        
        % do some basic error checking
        if isempty(ch_green) && ~isempty(ch_green);
            error('One channel is defined but not the other')
        end
        if all([ch_green(:);ch_red(:);ch_blue(:)]==0); continue; end % all channels are undefined... no big deal so no error
        
        
        % add images to the stack, make up a fake blue channel. During
        % image acquisition, the microscope objective flips the image L/R
        % and U/D, so flip all of them back...
        merge = cat(3, ch_red, ch_green, ch_blue);
        for a = 1:3
            merge(:,:,a) = rot90(merge(:,:,a), 2);
        end
        stack.img{idx} = merge;
        
        % figure out where the slice was in the brain based off of the
        % slice thickness and the slice number
        mdbidx = regexp({mdb(:).name}', mouse);
        mdbidx = ~cellfun(@isempty, mdbidx);
        thickness = mdb(mdbidx).histo.thickness;
        slicesPerPlate = mdb(mdbidx).histo.slicesPerPlate;
        stack.loc(idx) = sum(slicesPerPlate(1:p-1).*thickness) + (sl.*thickness);

        % update the index
        idx = idx + 1;
        
    end
end


imgAlign(stack, info_tmp);
    






