fin % Clear the workspace/command window so that the user can see the feedback. 

%%%%% LOAD THE RETINOTOPY (ROI) DATA
% Open dialogue box for selecting the processed ROI file
[fname, path] = uigetfile({'*.mat', 'Related Files (*.mat)'},...
                           'Select the Processed ROI file', 'S:\Data');
loadpath = [path, filesep, fname];    

% If no file is selected,
% Present an error message (Code cannot proceed without data)
if isequal(fname, 0)
    errordlg(sprintf(' File not found.\n User selected "Cancel", a file MUST\n be selected for code to progress'),...
                     'File Error')                
% If a file is selected, load the data...
% Also, notify the User when the data has been loaded.
else
    load(loadpath)
    fprintf('Data has been loaded.\n')
end


%%%%% PLOT SELECTED-ROIs OVER VASCULATURE
% Open dialogue box for selecting the raw vasculature image file
[filename, pathname] = uigetfile({'*.tif', 'Related Files (*.tif)'},...
                              'Select the Vasculature file', 'S:\Data');
loadpath = [pathname, filesep, filename];

% Present an error message (Code cannot proceed without data)
if isequal(filename, 0)
    errordlg(sprintf(' File not found.\n User selected "Cancel", a file MUST\n be selected for code to progress'),...
                     'File Error')                
% If a file is selected, load the data, etc...
else  
    % Info is a structure array with one element for each image in the file
    data = imfinfo(loadpath);
    num_images = numel(data); % Total # of Images
    
    VasculatureImg = zeros(data.Width, data.Height, num_images);
    VasculatureImg = imread(loadpath, 'Info', data);
    
    % Read in Intended Image
    figure
    imagesc(VasculatureImg); colormap gray
    
    % Naming the Figure (from the data file)
    figName = regexpi(filename, '(\w\d+)', 'match');
    figName = regexpi(figName{1}, '(\d)', 'match');
    figName = ['k', figName{end-1}, figName{end}];
    set(gcf, 'Name', sprintf('%s Retinotopy', figName), 'NumberTitle', 'off')
    
    % Other figure formatting...
    set(gcf, 'Units', 'normalized', 'Position', [0.2  0.05  0.5  0.75])
    set(gca, 'Position', [0.05  0.05  0.9  0.9], 'XTick', [], 'YTick', []);
    
    maskedImage = repmat(VasculatureImg, [1,1,3]);
    maskedImage = rgb2gray(maskedImage);
    
    % Plot the ROI Boundaries over the Vasculature
    for i_ROI = 1:numel(udat.ROI.RCidx) % For each ROI
        
        rc = udat.ROI.boundaries{i_ROI};
        if ~isempty(rc) && ~any(isnan(rc(:)))
            hold on;	% Don't blow away the image.
            plot(rc(:,2), rc(:,1), 'k-', 'LineWidth', 2);
        end
    end
end


%%%%% SANITY CHECK: Display Selected-ROIs

% Determine Elevation according to selected file
El = regexpi(regexpi(fname, 'el[(.*?)\]', 'match'), '[(.*?)\]', 'match');
El = El{:};
El = regexpi(El, '.*?\d+', 'match');
El = str2num(El{:}{:}(2:end));

% Determine Azimuth according to selected file
Az = regexpi(regexpi(fname, 'az[(.*?)\]', 'match'), '[(.*?)\]', 'match');
Az = Az{:};
Az = regexpi(Az, '.*?\d+', 'match');
Az = str2num(Az{:}{:}(2:end));

% Find the final_img according to Stim-position (specified by fname)
NStim = size(udat.ttypes, 1);
for i_Stim = 1: NStim
    if udat.ttypes(i_Stim) == El & udat.ttypes(i_Stim + NStim) == Az
        Img_ttypes = i_Stim;
    end
end

% Plot the final_img with 'colored' ROIs
figure
set(gcf, 'Name', sprintf('Circled ROIs (Az: %.0f,  El: %.0f)', Az, El), 'NumberTitle', 'off');

plotimg = udat.final_img{Img_ttypes};
plotimg = plotimg + (abs(min(plotimg(:))));
plotimg = plotimg ./ (max(plotimg(:)).*1.2);
plotimg = repmat(plotimg, [1,1,3]);
for i_roi = 1:numel(udat.ROI.RCidx)
    
    rc = udat.ROI.RCidx{i_roi};
    if ~isempty(rc) && ~any(isnan(rc(:)))
        linind = sub2ind(size(plotimg(:,:,1)), rc(:,1), rc(:,2));
        plotimg(linind) = 0;
    end
end

h_img = imshow(plotimg);