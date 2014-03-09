function stackMaker(mName, objective)

%   stackMaker(mname, objective)
%
% Example of what the input structure ('stack') should contain:
%   stack.mouse = 'CH_112613_B';
%   stack.objective = '2x';

% Flow of ideas:
%
% 1) grab all of the monochrome images. Put them in temporary structure
% ("raw"). Have different arrays for the RGB channels. If one channel was
% not aquired, than just shove a bunch of zeros in the array.
%       raw{1}.red.info => one for each channel
%       raw{1}.red.img
%
% 2) Display: 3 columns on the left (one for each RGB channel). In the
% middle, put all the contrast controls. On the right, put the merged
% image, and buttons to advance to the previous/next image. Add buttons for
% flips (UD and LR). Add a button to save export
%
% 3) Save/export should make a .mat struct that can be read in easily by
% other matlab scripts, and should save a tiff for things like imageJ.
% Alternatively, I could write a stand alone script to do this conversion
% (e.g., matStack2Tiff.m)
%
% 4) contrast controls should include:
%    ** adjust peak value (adjust the highest value on the LUT)
%    ** adjust minimum value (adjusts the lowest value of the LUT)
%    ** adjust the "brightness" of the LUT (moves the position of the LUT)
%    ** adjust the "contrast" of the LUT (piviots the LUT around some point
%       and changes the LUT slope)
%    ** down the line I could make some auto adjusters (imadjust,
%       adapthisteq, etc) as radio buttons.
%
% 5) Since there's only one set of contrast controls, I'll need to click on
% a "raw" image to make it "active". The active image should be outlined in
% red or something...
%



%
%  LOAD IN THE RAW IMAGES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pull in the raw images. put the raw images in the userData structure
udat.raw = img_getRaw(mName, objective);

% Initalize the GUI
udat = gui_init(udat);



end



function raw = img_getRaw(mname, objective)

    % cd to where the images are
    global GL_DATPATH
    cd([GL_DATPATH, filesep, mname, filesep, 'Histology', filesep, 'Raw Images']);

    % grab the names in the directory
    d = dir;


    % initialize the structure of images. "raw" will be nPlates x nSlices
    raw = {};


    for a = 1:numel(d);

        % display the progress
        if a == 1 || ~rem(a,5)
            fprintf('%d more images to unpack\n', numel(d)-(a-1));
        end

        if ~any(regexpi(d(a).name, '^[\.]|thumbs')) % skip the hidden files


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
                if regexpi(d(a).name, '_red'); color = 'red'; end
                if regexpi(d(a).name, '_green'); color = 'green'; end
                if regexpi(d(a).name, '_blue');   color = 'blue'; end
                if regexpi(d(a).name, '_white'); continue; end


                % unpack the images
                img = imread(d(a).name);
                img = rot90(img, 2); % rotate the image back from the inverted microscope image
                info = imfinfo(d(a).name);
                switch info.ColorType
                    case 'truecolor'
                        npix = 0;
                        img = preProcessImg(img, info, npix, 'none');
                    case 'grayscale'
                        % no need to do anything, already in grayscale
                end


                % put the image in a structure according to it's position in
                % the brain, and the color channel
                raw{plate, slice}.(color).img = img;
                raw{plate, slice}.(color).info = info;
            end
        end
    end % loop over directory contents
    
    
    % Kill the plate/slice indexing, b/c it isn't necessary anymore, and
    % creates headaches later
    raw = raw(:);
    
    % Roll through the raw images, and create blank entries for the
    % non-existant color channels
    for a = 1:numel(raw)
        colors = {'red', 'green', 'blue'};
        for j = 1:numel(colors)
            if ~isfield(raw{a}, colors{j})
                raw{a}.(colors{j}).img = zeros(size(img));
                raw{a}.(colors{j}).info = info;
            end
        end
    end
    
end




function udat = gui_init(udat)

    % Initialize the GUI, plot the first images to the RGB preview pannels,
    % then plot the merge. By default, the histos will contain the first
    % available channel (likely red or green)

    %
    % open figure
    %
    udat.h.fig = figure;
    set(udat.h.fig, 'units', 'normalized', 'position', [0.0708    0.0044    0.8146    0.8900])

    %
    % initialize three axes for the 'raw' images
    %
    udat.h.axBlue = axes('position', [0.02, 0.02, .31, .31]);
    set(udat.h.axBlue, 'xtick', [],...
                       'ytick', [],...
                       'box', 'on');
    udat.h.imgBlue = imshow(ones(size(udat.raw{1}.blue.img)), 'parent', udat.h.axBlue);
    set(udat.h.imgBlue, 'ButtonDownFcn', {@gui_selectChannel});


    udat.h.axGreen = axes('position', [0.02, 0.35, .31, .31]);
    set(udat.h.axGreen, 'xtick', [],...
                       'ytick', [],...
                       'box', 'on');
    udat.h.imgGreen = imshow(ones(size(udat.raw{1}.green.img)), 'parent', udat.h.axGreen);
    set(udat.h.imgGreen, 'ButtonDownFcn', {@gui_selectChannel});

                   
                   
    udat.h.axRed = axes('position', [0.02, 0.68, .31, .31]);
    set(udat.h.axRed, 'xtick', [],...
                       'ytick', [],...
                       'box', 'on');
    udat.h.imgRed = imshow(ones(size(udat.raw{1}.red.img)), 'parent', udat.h.axRed);
    set(udat.h.imgRed, 'ButtonDownFcn', {@gui_selectChannel});

                   
                   
    
    %
    % initialize the axis for the RGB merged image
    %
    udat.h.axMerge = axes('position', [0.4, 0.48, .5, .5]);
    set(udat.h.axMerge, 'xtick', [], 'ytick', [], 'box', 'on')
    udat.h.imgMerge = imshow(ones(size(udat.raw{1}.red.img)), 'parent', udat.h.axMerge);
    
    
    %
    % initialize the UI controls for various other things
    %
    udat.h.fliplr = uicontrol('style', 'togglebutton',...
                         'units', 'normalized',...
                         'string', 'FLIP L/R',...
                         'Position', [0.72, 0.30, 0.10, 0.05],...
                         'Callback', {@img_flip});
                     
    udat.h.flipud = uicontrol('style', 'togglebutton',...
                         'units', 'normalized',...
                         'string', 'FLIP U/D',...
                         'Position', [0.85, 0.30, 0.10, 0.05],...
                         'Callback', {@img_flip});
                     
    udat.h.slider = uicontrol('style', 'slider',...
                         'units', 'normalized',...
                         'position', [0.73, 0.4, 0.2, 0.05],...
                         'Callback', {@gui_slider},...
                         'Value', 1,...
                         'max', numel(udat.raw),...
                         'min', 1,...
                         'SliderStep', [1./numel(udat.raw), 1./numel(udat.raw).*10]);
                         
    udat.h.textbox = uicontrol('style', 'edit',...
                          'units', 'normalized',...
                          'string', sprintf('1 of %d', numel(udat.raw)),...
                          'position', [0.76, 0.36, 0.15, 0.055],...
                          'callback', {@gui_txtUpdate});



    % Put the udat structure in the UserData field.
    set(udat.h.fig, 'UserData', udat);


end


function gui_selectChannel(hand, ~)
    
    % grab the user data
    udat = get(gcf, 'userdata');
    
    % make all the boxes narrow lines (will bold the current ax later)
    set([udat.h.axRed, udat.h.axGreen, udat.h.axBlue], 'linewidth', 0.5)
    
    % select the channel when the user clicks on an axis.
    switch hand
        case udat.h.axRed
            set(udat.h.axRed, 'linewidth', 12)
            udat.currentColor = 'red';
        case udat.h.axGreen
            set(udat.h.axGreen, 'linewidth', 12)
            udat.currentColor = 'green';
        case udat.h.axBlue
            set(udat.h.axBlue, 'linewidth', 12)
            udat.currentColor = 'blue';
    end
    
    
    % re-set the user data
    set(gcf, 'UserData', udat)
    
end


function gui_slider(varargin)
    
    % grab the data
    udat = get(gcf, 'userdata');
    
    % figure out the value of the slider
    sliceNum = round(get(udat.h.slider, 'Value'));
    
    % update the text field of the textbox
    set(udat.h.textbox, 'String', sprintf('%d of %d', sliceNum, numel(udat.raw)));
        
    % set the user data
    set(gcf, 'userdata', udat);
    
    % update the image
    img_updateAxes
    
end


function gui_txtUpdate(varargin)
    

    % grab the data
    udat = get(gcf, 'userdata');
    
    % grab the value indicated in the text box
    txtString = get(udat.h.textbox, 'String');
    txtString = regexp(txtString, '[\d]+', 'match'); % find the number in the text string
    txtString = txtString{1}; % only consider the first (b/c the second entry could be the total number of slices)
    sliceNum = str2double(txtString);
    
    % do some error checking
    if sliceNum <= 0 || sliceNum > numel(udat.raw)
        oldSliceNum = get(udat.h.slider, 'value');
        set(udat.h.textbox, 'String', sprintf('%d of %d', oldSliceNum, numel(udat.raw)));
        return
    end
    
    % update the slider position and the text of the text box
    set(udat.h.slider, 'Value', sliceNum)
    set(udat.h.textbox, 'String', sprintf('%d of %d', sliceNum, numel(udat.raw)));
    
    % set the user data
    set(gcf, 'userdata', udat);
    
    % update the image
    img_updateAxes
    
end




function img_flip(varargin)
    
    % grab the data
    udat = get(gcf, 'userdata');

    % get the value of the slider
    sliceNum = round(get(udat.h.slider, 'Value'));
    
    %loop over the colors (using dynamic indexing)
    colors = {'red', 'green', 'blue'};
    for a = 1:numel(colors);
        
        subfield = colors{a};
        
        oldImg = udat.raw{sliceNum}.(subfield).img;
        
        switch varargin{1} % the handle to the uibutton
            case udat.h.fliplr
                udat.raw{sliceNum}.(subfield).img = fliplr(oldImg);
            case udat.h.flipud
                udat.raw{sliceNum}.(subfield).img = flipud(oldImg);
            otherwise
                error('IMG_FLIP: unknown UIcontrol handel')
        end
    end
        
        
    % save the new data
    set(gcf, 'userdata', udat);
    
    % update the image
    img_updateAxes
    
end


function img_updateAxes()

    % grab the data
    udat = get(gcf, 'userdata');

    % get the value of the slider
    sliceNum = round(get(udat.h.slider, 'Value'));
    
    % update each raw axis
    set(udat.h.axRed, 'CData', udat.raw{sliceNum}.red.img)
    set(udat.h.axGreen, 'CData', udat.raw{sliceNum}.green.img)
    set(udat.h.axBlue, 'CData', udat.raw{sliceNum}.blue.img)


end



