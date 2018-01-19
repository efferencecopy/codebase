%%  create an image library for registration
fin

% ----- Load images ----------
img_store = imageDatastore('\\crash.dhe.duke.edu\charlie\New Histology Photos\test_stitch');

% Display images to be stitched
figure
montage(img_store.Files)
drawnow


% ------------- Identify registration points --------------

% define how features will be detected
detectFeatures = @(x) detectSURFFeatures(x, 'MetricThreshold', 5, 'NumScaleLevels', 4, 'NumOctaves', 1);

% Read the first image from the image set.
grayImage = readimage(img_store, 1);

% Initialize features for I(1)
points = detectFeatures(grayImage);
[features, valid_points] = extractFeatures(grayImage, points);

% Initialize all the transforms to the identity matrix. Note that the
% projective transform is used here because the building images are fairly
% close to the camera. Had the scene been captured from a further distance,
% an affine transform would suffice.
numImages = numel(img_store.Files);

if exist('tforms', 'var')
    clear tforms
end
switch 'affine'
    case 'affine'
        tforms(numImages) = affine2d(eye(3));
    case 'projective'
        tforms(numImages) = projective2d(eye(3));
end

% Iterate over remaining image pairs
for n = 2:numImages

    % Store points and features for I(n-1).
    pointsValidPrevious = valid_points;
    featuresPrevious = features;
    previous_img = grayImage;

    % Read I(n).
    grayImage = readimage(img_store, n);

    % Detect and extract SURF features for I(n).
    points = detectFeatures(grayImage);
    [features, valid_points] = extractFeatures(grayImage, points);

    % Find correspondences between I(n) and I(n-1).
    indexPairs = matchFeatures(features, featuresPrevious);

    matchedPoints = valid_points(indexPairs(:,1));
    matchedPointsPrev = pointsValidPrevious(indexPairs(:,2));


    % Estimate the transformation between I(n) and I(n-1).
    [tmp_tform, inlierpts1, inlierpts2, status]  = estimateGeometricTransform(matchedPoints,...
                                           matchedPointsPrev,...
                                           'similarity',...
                                           'Confidence', 99.9,...  % default 99
                                           'MaxNumTrials', 2000,...  % default 1000
                                           'MaxDistance', 1.5);       % default 1.5 pixels
                                           
%     figure
%     showMatchedFeatures(grayImage, previous_img, inlierpts1, inlierpts2)
%     legend('ptsOne', 'ptsTwo')
%     drawnow
    
    
    % Compute T(n) * T(n-1) * ... * T(1)
    tforms(n) = tmp_tform;
    tforms(n).T = tforms(n).T * tforms(n-1).T;
end


%%

imageSize = size(grayImage);  % all the images are the same size

% Compute the output limits  for each transform
for i = 1:numel(tforms)
    [xlim(i,:), ylim(i,:)] = outputLimits(tforms(i), [1 imageSize(2)], [1 imageSize(1)]);
end

% Next, compute the average X limits for each transforms and find the image that is in the center. Only the X limits are used here because the scene is known to be horizontal. If another set of images are used, both the X and Y limits may need to be used to find the center image.

avgXLim = mean(xlim, 2);

[~, idx] = sort(avgXLim);

centerIdx = floor((numel(tforms)+1)/2);

centerImageIdx = idx(centerIdx);

% Finally, apply the center image's inverse transform to all the others.

Tinv = invert(tforms(centerImageIdx));

for i = 1:numel(tforms)
    tforms(i).T = tforms(i).T * Tinv.T;
end

% Step 3 - Initialize the Panorama
% 
% Now, create an initial, empty, panorama into which all the images are mapped.
% 
% Use the outputLimits method to compute the minimum and maximum output limits over all transformations. These values are used to automatically compute the size of the panorama.

for i = 1:numel(tforms)
    [xlim(i,:), ylim(i,:)] = outputLimits(tforms(i), [1 imageSize(2)], [1 imageSize(1)]);
end

% Find the minimum and maximum output limits
xMin = min([1; xlim(:)]);
xMax = max([imageSize(2); xlim(:)]);

yMin = min([1; ylim(:)]);
yMax = max([imageSize(1); ylim(:)]);

% Width and height of panorama.
width  = round(xMax - xMin);
height = round(yMax - yMin);

% Initialize the "empty" panorama.
panorama = zeros([height width], 'like', grayImage);

% Step 4 - Create the Panorama
% 
% Use imwarp to map images into the panorama and use vision.AlphaBlender to overlay the images together.

blender = vision.AlphaBlender('Operation', 'Binary mask', ...
    'MaskSource', 'Input port');

% Create a 2-D spatial reference object defining the size of the panorama.
xLimits = [xMin xMax];
yLimits = [yMin yMax];
panoramaView = imref2d([height width], xLimits, yLimits);

% Create the panorama.
for i = 1:numImages

    I = readimage(img_store, i);

    % Transform I into the panorama.
    warpedImage = imwarp(I, tforms(i), 'OutputView', panoramaView);

    % Generate a binary mask.
    mask = imwarp(true(size(I,1),size(I,2)), tforms(i), 'OutputView', panoramaView);

    % Overlay the warpedImage onto the panorama.
    panorama = step(blender, panorama, warpedImage, mask);
end

figure
imshow(panorama)
