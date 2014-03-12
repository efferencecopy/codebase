function stackMaker(mName, objective)

%   stackMaker(mname, objective)
%


%  TO DO
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
% 5) Output a mapping between image number and position in the brain (plate
% and slice)


% KNOWN ISSUES
%
% 1) the images from the Nikon camera don't come through, but I'm not sure
% why... 
%
% 2) When a slice does not exist (or there's no picture from that slice) I
% need to omit that 'raw' and 'merge' cell from the analysis.
%


%
%  LOAD IN THE RAW IMAGES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pull in the raw images. put the raw images in the userData structure
udat.raw = img_getRaw(mName, objective);

% initialize the LUTs for the merged image
udat.merge = merge_initLUT(udat);

% Initalize the GUI
udat = gui_init(udat);



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
    tmp_img = udat.raw{1}.blue.img;
    lims = [udat.raw{1}.blue.info.MinSampleValue, udat.raw{1}.blue.info.MaxSampleValue];
    udat.h.axBlue = axes('position', [0.02, 0.02, .31, .31]);
    udat.h.imgBlue = imshow(tmp_img, lims, 'parent', udat.h.axBlue);
    set(udat.h.imgBlue, 'ButtonDownFcn', {@gui_selectChannel});
    set(udat.h.axBlue, 'visible', 'on',...
                       'xtick', [],...
                       'ytick', [],...
                       'box', 'on');

    
    tmp_img = udat.raw{1}.green.img;
    lims = [udat.raw{1}.green.info.MinSampleValue, udat.raw{1}.green.info.MaxSampleValue];
    udat.h.axGreen = axes('position', [0.02, 0.35, .31, .31]);
    udat.h.imgGreen = imshow(tmp_img, lims, 'parent', udat.h.axGreen);
    set(udat.h.imgGreen, 'ButtonDownFcn', {@gui_selectChannel});
    set(udat.h.axGreen, 'visible', 'on',...
                        'xtick', [],...
                        'ytick', [],...
                        'box', 'on');
                   
                   
    tmp_img = udat.raw{1}.green.img;
    lims = [udat.raw{1}.green.info.MinSampleValue, udat.raw{1}.green.info.MaxSampleValue];
    udat.h.axRed = axes('position', [0.02, 0.68, .31, .31]);
    udat.h.imgRed = imshow(tmp_img, lims, 'parent', udat.h.axRed);
    set(udat.h.imgRed, 'ButtonDownFcn', {@gui_selectChannel});
    set(udat.h.axRed,  'visible', 'on',...
                       'xtick', [],...
                       'ytick', [],...
                       'box', 'on');
                   
                   
    
    %
    % initialize the axis for the RGB merged image
    %
    udat.h.axMerge = axes('position', [0.4, 0.48, .5, .5]);
    set(udat.h.axMerge, 'xtick', [], 'ytick', [], 'box', 'on')
    udat.h.imgMerge = imshow(tmp_img, lims, 'parent', udat.h.axMerge);
    
    
    %
    % initialize the UI controls for various other things
    %
    udat.h.fliplr = uicontrol('style', 'togglebutton',...
                         'units', 'normalized',...
                         'string', 'FLIP L/R',...
                         'Position', [0.72, 0.265, 0.10, 0.05],...
                         'Callback', {@img_flip});
                     
    udat.h.flipud = uicontrol('style', 'togglebutton',...
                         'units', 'normalized',...
                         'string', 'FLIP U/D',...
                         'Position', [0.85, 0.265, 0.10, 0.05],...
                         'Callback', {@img_flip});
                     
    udat.h.slider = uicontrol('style', 'slider',...
                         'units', 'normalized',...
                         'position', [0.73, 0.4, 0.2, 0.03],...
                         'Callback', {@gui_slider},...
                         'Value', 1,...
                         'max', numel(udat.raw),...
                         'min', 1,...
                         'SliderStep', [1./numel(udat.raw), 1./numel(udat.raw).*10]);
                         
    udat.h.textbox = uicontrol('style', 'edit',...
                          'units', 'normalized',...
                          'string', sprintf('1 of %d', numel(udat.raw)),...
                          'position', [0.76, 0.34, 0.15, 0.05],...
                          'callback', {@gui_txtUpdate});


    %
    % Create the GUI controls for the contrast adjustments and the
    % histogram viewer
    %
    udat.h.axHist = axes('position', [0.38, 0.29, .25, .15]);
    set(udat.h.axHist, 'tickDir', 'out',...
                       'ytick', [],...
                       'box', 'off')
    set(get(udat.h.axHist, 'title'), 'string', 'Histogram of RAW DAC values')
    set(udat.h.axHist, 'children', bar(1:10, ones(1,10))); % initialize the histogram as a child of the axis
    udat.h.barHist = get(udat.h.axHist, 'children');
    
        
    % add a second axis on top of the histo, for the LUT, make this axis
    % have no linewidth or background.
    udat.h.axLUT = axes('position', get(udat.h.axHist, 'position'));
    set(udat.h.axLUT, 'children', plot(1:10, 1:10, 'b', 'linewidth', 2)); % initialize the histogram as a child of the axis
    udat.h.lineLUT = get(udat.h.axLUT, 'children');
    set(udat.h.axLUT, 'box', 'off',...
                      'xtick', [],...
                      'ytick', [],...
                      'color', 'none',...
                      'ButtonDownFcn', {@img_resetLUT})
    
    udat.h.lutMaxSlider = uicontrol('style', 'slider',...
                         'units', 'normalized',...
                         'position', [0.37, 0.20, 0.27, 0.02],...
                         'Callback', {@img_adjustLUTvals},...
                         'Value', 1,...
                         'max', 10,...
                         'min', 1,...
                         'SliderStep', [0.01, 0.1]);
                     
    udat.h.lutMinSlider = uicontrol('style', 'slider',...
                         'units', 'normalized',...
                         'position', [0.37, 0.15, 0.27, 0.02],...
                         'Callback', {@img_adjustLUTvals},...
                         'Value', 1,...
                         'max', 10,...
                         'min', 1,...
                         'SliderStep', [0.01, 0.1]);

    udat.h.lutDRangeSlider = uicontrol('style', 'slider',...
                         'units', 'normalized',...
                         'position', [0.37, 0.10, 0.27, 0.02],...
                         'Callback', {@img_adjustLUTvals},...
                         'Value', 1,...
                         'max', 10,...
                         'min', 1,...
                         'SliderStep', [0.01, 0.1]);
                  
                  
                  

                      
    % Put the udat structure in the UserData field.
    set(udat.h.fig, 'UserData', udat);
    
    % plot the first set of images
    gui_selectChannel(udat.h.imgRed)
    img_updateHistogram % also updates the raw/merged images by calling img_updateImages
    
end


function gui_selectChannel(hand, ~)
    
    % grab the user data
    udat = get(gcf, 'userdata');
    
    % make all the boxes narrow lines (will bold the current ax later)
    set([udat.h.axRed, udat.h.axGreen, udat.h.axBlue],...
        'ycolor', 'k',...
        'xcolor', 'k',...
        'linewidth', 0.25)
    
    % select the channel when the user clicks on an axis.
    switch hand
        case {udat.h.axRed, udat.h.imgRed}
            set(udat.h.axRed, 'ycolor', 'r', 'xcolor','r', 'linewidth', 3)
            udat.currentColor = 'red';
        case {udat.h.axGreen, udat.h.imgGreen}
            set(udat.h.axGreen, 'ycolor', 'g', 'xcolor','g', 'linewidth', 3)
            udat.currentColor = 'green';
        case {udat.h.axBlue, udat.h.imgBlue}
            set(udat.h.axBlue, 'ycolor', 'b', 'xcolor','b', 'linewidth', 3)
            udat.currentColor = 'blue';
    end
    
    
    % re-set the user data
    set(gcf, 'UserData', udat)
    
    % update the histogram of DAC values and the LUT.
    img_updateHistogram
    gui_updateContrastSliders
    
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
    
    % update the raw images and the histogram viewer. Do this by calling
    % img_updateHistogram (which calls img_updateImages after updating the LUT).
    img_updateHistogram
    gui_updateContrastSliders()
    
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
        oldSliceNum = round(get(udat.h.slider, 'value'));
        set(udat.h.textbox, 'String', sprintf('%d of %d', oldSliceNum, numel(udat.raw)));
        return
    end
    
    % update the slider position and the text of the text box
    set(udat.h.slider, 'Value', sliceNum)
    set(udat.h.textbox, 'String', sprintf('%d of %d', sliceNum, numel(udat.raw)));
    
    % set the user data
    set(gcf, 'userdata', udat);
    
    % update the raw images and the histogram viewer. Do this by calling
    % img_updateHistogram (which calls img_updateImages after updating the LUT).
    img_updateHistogram
    gui_updateContrastSliders()
    
end

function gui_updateContrastSliders(varargin)
    % updates the position of the contrast sliders based off the
    % udat.merge.lut values (which are set by the uicallback function)
    
    % grab the userdata, and some other useful info
    udat = get(gcf, 'userdata');
    sliceNum = round(get(udat.h.slider, 'Value'));
    color = udat.currentColor;
    ch_maxdac = udat.raw{sliceNum}.(color).info.MaxSampleValue;
    ch_mindac = udat.raw{sliceNum}.(color).info.MinSampleValue;
    
    
    % dynamic range slider will be centered in the middle of the LUT
    lut_max = udat.merge{sliceNum}.(color).lut_hi;
    lut_min = udat.merge{sliceNum}.(color).lut_low;
    halfwidth = round(((lut_max-lut_min)/2));
    dRangeVal = halfwidth+lut_min;
    set(udat.h.lutDRangeSlider, 'Value', dRangeVal,...
                                'max', ch_maxdac+halfwidth,...
                                'min', ch_mindac-halfwidth,...
                                'SliderStep', [0.01, 0.15]);
    
    
    % max val slider
    maxval = max(ch_maxdac, dRangeVal+halfwidth);
    ch_setVal = udat.merge{sliceNum}.(color).lut_hi;
    set(udat.h.lutMaxSlider, 'Value', ch_setVal,...
                         'max', maxval,...
                         'min', ch_mindac,...
                         'SliderStep', [0.005, 0.15]);
                     
    % min val slider
    minval = min(ch_mindac, dRangeVal-halfwidth);
    ch_setVal = udat.merge{sliceNum}.(color).lut_low;
    set(udat.h.lutMinSlider, 'Value', ch_setVal,...
                             'max', ch_maxdac,...
                             'min', minval,...
                             'SliderStep', [0.005, 0.15]);
                         
 

    
    
    
    
    
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
                info = imfinfo(d(a).name);
                switch info.ColorType
                    case 'truecolor'
                        npix = 0;
                        img = preProcessImg(img, info, npix, 'none');
                        info.MaxSampleValue = info.MaxSampleValue(1);
                        info.MinSampleValue = info.MinSampleValue(1);
                    case 'grayscale'
                        % no need to do anything, already in grayscale
                end

                %rotate the image back from the inverted microscope image
                img = rot90(img, 2);
                
                % put the image in a structure according to it's position in
                % the brain, and the color channel
                raw{plate, slice}.(color).img = img;
                raw{plate, slice}.(color).info = info;
            end
        end
    end % loop over directory contents
    
    
    % Kill the plate/slice indexing, b/c it isn't necessary anymore, and
    % creates headaches later
    raw = raw';
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


function img_adjustLUTvals(hand, ~)

    % grab the userData
    udat = get(gcf, 'userdata');
    sliceNum = round(get(udat.h.slider, 'Value'));
    color = udat.currentColor;
    
    % gets called when the user toggles a contrast UI slider
    switch hand
        case udat.h.lutMaxSlider
            % grab the desired value
            setVal = get(udat.h.lutMaxSlider, 'value');
            setVal = round(setVal);
            
            % don't let the lut_high val be less than the low val
            if setVal <= udat.merge{sliceNum}.(color).lut_low
                gui_updateContrastSliders
                return
            end
            
            % set the new lut_high value
            udat.merge{sliceNum}.(color).lut_hi = setVal;            
            
            % linear algebra solution to the new slope and yint
            x1 = udat.merge{sliceNum}.(color).lut_low;
            x2 = udat.merge{sliceNum}.(color).lut_hi;            
            y1 = udat.raw{sliceNum}.(color).info.MinSampleValue;
            y2 = udat.raw{sliceNum}.(color).info.MaxSampleValue;
            pred = [x1, 1; x2, 1];
            resp = [y1;y2];
            betas = pred\resp;
            udat.merge{sliceNum}.(color).lut_slope = betas(1);
            udat.merge{sliceNum}.(color).lut_yint = betas(2);
            
            
        case udat.h.lutMinSlider
            
            %grab the desired value
            setVal = get(udat.h.lutMinSlider, 'value');
            setVal = round(setVal);
            
            % make sure that the lut_low val is less than the high val
            if setVal >= udat.merge{sliceNum}.(color).lut_hi
                gui_updateContrastSliders % set things back the way they were
                return
            end
            
            % assign the new value
            udat.merge{sliceNum}.(color).lut_low = setVal;
                       
            % linear algebra solution to the new slope and yint
            x1 = udat.merge{sliceNum}.(color).lut_low;
            x2 = udat.merge{sliceNum}.(color).lut_hi;
            y1 = udat.raw{sliceNum}.(color).info.MinSampleValue;
            y2 = udat.raw{sliceNum}.(color).info.MaxSampleValue;
            pred = [x1, 1; x2, 1];
            resp = [y1;y2];
            betas = pred\resp;
            udat.merge{sliceNum}.(color).lut_slope = betas(1);
            udat.merge{sliceNum}.(color).lut_yint = betas(2);
            
        case udat.h.lutDRangeSlider
            
            % where is the slider now?
            currentPos = round(get(udat.h.lutDRangeSlider, 'value'));
            
            % where was it b/4 the user changed things?
            maxlut = udat.merge{sliceNum}.(color).lut_hi;
            minlut = udat.merge{sliceNum}.(color).lut_low;
            assert(maxlut>minlut, 'LUT ERROR: max val must be greater than min val');
            halfwidth = round(((maxlut-minlut)/2));
            oldval = halfwidth+minlut;
            
            % figure out how much things moved, then adjust the LUT values
            delta = oldval - currentPos;
            maxlut_new = maxlut - delta;
            minlut_new = minlut - delta;
            udat.merge{sliceNum}.(color).lut_hi = maxlut_new;
            udat.merge{sliceNum}.(color).lut_low = minlut_new;
            
            % linear algebra solution to the new slope and yint
            x1 = udat.merge{sliceNum}.(color).lut_low;
            x2 = udat.merge{sliceNum}.(color).lut_hi;
            y1 = udat.raw{sliceNum}.(color).info.MinSampleValue;
            y2 = udat.raw{sliceNum}.(color).info.MaxSampleValue;
            pred = [x1, 1; x2, 1];
            resp = [y1;y2];
            betas = pred\resp;
            udat.merge{sliceNum}.(color).lut_slope = betas(1);
            udat.merge{sliceNum}.(color).lut_yint = betas(2);
            
            
        case udat.h.cntSlopeSlider
    end
    
    % re-set the user data field
    set(gcf, 'userdata', udat);
    
    % update the histogram and LUT, and then show the new (adjusted) images
    img_updateHistogram
    
    % if you change one of the contrast sliders, you could effect the value
    % of the others... Update accordingly.
    gui_updateContrastSliders
    
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
    img_updateImages
    
end


function img_updateImages()

    % grab the data
    udat = get(gcf, 'userdata');

    % get the value of the slider
    sliceNum = round(get(udat.h.slider, 'Value'));
    
    % update each raw axis
    set(udat.h.imgRed, 'CData', img_applyLUT(udat.raw{sliceNum}.red, udat.merge{sliceNum}.red))
    set(udat.h.imgGreen, 'CData', img_applyLUT(udat.raw{sliceNum}.green, udat.merge{sliceNum}.green))
    set(udat.h.imgBlue, 'CData', img_applyLUT(udat.raw{sliceNum}.blue, udat.merge{sliceNum}.blue))

    % update the merged image
    tmp_cat = cat(3, get(udat.h.imgRed, 'CData'), get(udat.h.imgGreen, 'CData'), get(udat.h.imgBlue, 'CData'));
    set(udat.h.imgMerge, 'CData', uint16(tmp_cat)); % needs to be uint16 inorder for imshow to deal with it appropriately
                                     
end

function img_updateHistogram()
    
    % grab the userdat, the slice number, and the color channel
    udat = get(gcf, 'userdata');
    sliceNum = round(get(udat.h.slider, 'Value'));
    color = udat.currentColor;
    
    % set the new LUT values based on the GUI controls
    
    % prepare the histogram
    raw = udat.raw{sliceNum}.(color).img(:);
    minval = udat.raw{sliceNum}.(color).info.MinSampleValue;
    maxval = udat.raw{sliceNum}.(color).info.MaxSampleValue;
    edges = linspace(minval, maxval, 150);
    counts = histc(raw, edges);
    
    % display the histogram
    set(udat.h.barHist, 'xdata', edges, 'ydata', counts);
    set(udat.h.axHist, 'Xlim', [minval-10 maxval],...
                       'Ylim', [0 max(counts)],...
                       'ytick', [],...
                       'box', 'on')
    set(get(udat.h.axHist, 'title'), 'string', 'Histogram of RAW DAC values')
    
    % display the LUT
    xx = minval:maxval;
    m = udat.merge{sliceNum}.(color).lut_slope;
    b = udat.merge{sliceNum}.(color).lut_yint;
    yy = m.*xx+b;
    yy(yy<0) = 0;
    yy(yy>maxval) = maxval;
    set(udat.h.lineLUT, 'xdata', xx, 'ydata', yy)
    set(udat.h.axLUT, 'box', 'off',...
                      'xtick', [],...
                      'ytick', [],...
                      'color', 'none',...
                      'Xlim', [minval-10 maxval],...
                      'Ylim', [0 maxval])
                  
  
    % store the userdata
    set(gcf, 'userdata', udat);
    
    % update the images (raw and merged) using the new LUT
    img_updateImages
    
end


function out = img_applyLUT(raw, merge)
    
    % apply the LUT
    xx = raw.info.MinSampleValue : raw.info.MaxSampleValue;
    LUT = merge.lut_slope .* xx + merge.lut_yint;
    idx = raw.img(:) + 1; % adding one so that indicies go from [1 to maxdac] and can be used as actual indicies
    out = LUT(idx);
    
    % make sure the values are positvie, and not greater than maxdac
    out(out<0) = 0;
    out(out>raw.info.MaxSampleValue) = raw.info.MaxSampleValue;
    
    % reshape to original dims
    out = reshape(out, size(raw.img));
    
end


function img_resetLUT(varargin)

    % check for double clicks
    persistent chk
    
    if isempty(chk)
        chk = 1;
        pause(0.3)
        if chk==1 % single click
            chk = [];
        end
    else %double click
        chk = [];
        
        % grab the userdata
        udat = get(gcf, 'userdata');
        sliceNum = round(get(udat.h.slider, 'Value'));
        color = udat.currentColor;
        
        % reset the LUT values for this color only
        udat.merge{sliceNum}.(color).lut_hi = udat.raw{sliceNum}.(color).info.MaxSampleValue;
        udat.merge{sliceNum}.(color).lut_low = udat.raw{sliceNum}.(color).info.MinSampleValue;
        udat.merge{sliceNum}.(color).lut_slope = 1;
        udat.merge{sliceNum}.(color).lut_yint = 0;
        
        % reset the userdata
        set(gcf, 'userdata', udat);
        
        % reset the sliders, and adjust the images
        img_updateHistogram
        gui_updateContrastSliders
        
    end
    
end



function merge = merge_initLUT(udat)

    for a = 1:numel(udat.raw)
        for color = {'red', 'green', 'blue'}
            merge{a}.(color{1}).lut_hi = udat.raw{a}.(color{1}).info.MaxSampleValue;
            merge{a}.(color{1}).lut_low = udat.raw{a}.(color{1}).info.MinSampleValue;
            merge{a}.(color{1}).lut_slope = 1;
            merge{a}.(color{1}).lut_yint = 0;
        end
    end
    
    % make it a column vector
    merge = merge(:);
    
end


