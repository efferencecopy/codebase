function stackViewer(mouseName, scalebar)

% TO DO

%  Change the call to plotimg to a direct call to imshow

% deal with optional inputs
if ~exist('scalebar', 'var')
    scalebar = 500;
end



%
% IMPORT THE DATA
% if mulitple stacks are present, force the user to choose. If no staks are
% present, than suggest that the user run stackMaker
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global GL_DATPATH
mouseDir = [GL_DATPATH, mouseName, filesep, 'Histology'];
cd(mouseDir)
files = dir('stack*');
if numel(files) == 1
    stackDir = files(1).name;
elseif numel(files) > 1
    stackDir = uigetdir([], 'Select A Stack');
elseif numel(files) == 0
    clc
    fprintf('*** No stacks were found for mouse ''%s'' ***\n', mouseName);
    fprintf('     To create a stack run the stackMaker function:\n');
    fprintf('     stackMaker(''%s'', ''2x'')\n', mouseName);
    return
end

% grab the images
stack = compileStack(stackDir, mouseName, scalebar);
    


%
% SET UP THE GUI BUTTONS
%
%%%%%%%%%%%%%%%%%%
h.fig = figure;
set(h.fig, 'position', [351 87 721 719])
h.fliplr = uicontrol('style', 'togglebutton',...
                     'units', 'normalized',...
                     'string', 'FLIP L/R',...
                     'Position', [0.1, 0.05, 0.25, 0.10],...
                     'Callback', {@img_fliplr});
                 
h.flipud = uicontrol('style', 'togglebutton',...
                     'units', 'normalized',...
                     'string', 'FLIP U/D',...
                     'Position', [0.4, 0.05, 0.25, 0.10],...
                     'Callback', {@img_flipud});
                 
h.slider = uicontrol('style', 'slider',...
                     'units', 'normalized',...
                     'position', [0.7, 0.08, 0.2, 0.05],...
                     'Callback', {@img_change},...
                     'Value', 1,...
                     'max', numel(stack.img),...
                     'min', 1,...
                     'SliderStep', [1./numel(stack.img), 1./numel(stack.img).*10]);
                     
h.textbox = uicontrol('style', 'edit',...
                      'units', 'normalized',...
                      'string', sprintf('1 of %d', numel(stack.img)),...
                      'position', [0.73, 0.02, 0.15, 0.055],...
                      'callback', {@txt_update});
                  


 % plot the first image to the screen. Shove the data into the 'userdata'
 % field
 h.img = plotimg(stack.img{1}, stack.info{1});
 drawnow
 
 %
 % make a plot of all the slices (stem plot) and make the lollipop for the
 % current slice slightly taller than the rest.
 %
 %%%%%%%%%%%%%%%%%%%%%%%%%%
 h.stem = axes;
 hgt = ones(numel(stack.loc),1);
 hgt(1) = 2;
 stem(stack.loc, hgt, 'k')
 set(h.stem, 'position', [.1, .87, .8, .1],...
             'box', 'off',...
             'ycolor', [1 1 1],...
             'ytick', [],...
             'ylim', [0, 2.3]);
 
 
 udat.h = h;
 udat.stack = stack;
 set(h.fig, 'UserData', udat)

end

function stack = compileStack(stackDir, mouseName, scalebar)
    
    % import some things from the MDB
    verbose = false;
    mdb = initMouseDB('update', verbose);
    [~, ind] = mdb.search(mouseName);
    slicesPerPlate = mdb.mice{ind}.histo.slicesPerPlate;
    thickness = mdb.mice{ind}.histo.thickness;
    

    idx = 1;
    cd(stackDir);
    d = dir;
    for a = 1:numel(d);

        % display the progress
        if a == 1 || ~rem(a,5)
            fprintf('%d more images to unpack\n', numel(d)-(a-1));
        end

        if any(regexpi(d(a).name, '^[\.]|thumbs')) % skip the hidden files
            continue
        end

        % id the plate number
        plate = regexpi(d(a).name , '_p\d+', 'match');
        plate = str2double(plate{1}(3:end));

        % id the slice number
        slice = regexpi(d(a).name , '_s\d+', 'match');
        slice = str2double(slice{1}(3:end));
        
        
        % get the merged image
        stack.img{idx} = imread(d(a).name);
        stack.info{idx} = imfinfo(d(a).name);
        
        % deal with the location
        if isscalar(slicesPerPlate)
            stack.loc(idx) = (slicesPerPlate.*(plate-1).*thickness) + ((slice-1).*thickness);
        else
            stack.loc(idx) = sum(slicesPerPlate(1:plate-1).*thickness) + ((slice-1).*thickness);
        end
        
        % update the idx
        idx = idx + 1;
    end
    
    % Make sure the images are in the correct order (poterior to anterior)
    [stack.loc, idx] = sort(stack.loc);
    stack.img = stack.img(idx);
    
    % add the scale bar
    xdim = size(stack.img{1}, 2);
    ydim = size(stack.img{2}, 1);
    pixperum = pixPerMicron(xdim, ydim);
    scalebar_img = ones(10, round(pixperum.*scalebar), 3);
    maxdac = stack.info{1}.MaxSampleValue(1);
    scalebar_img = scalebar_img .* maxdac;
    rows = 25:34;
    cols = 25:size(scalebar_img,2)+24;
    for a = 1:numel(stack.img)
        stack.img{a}(rows, cols, :) = scalebar_img;
    end
    
end


function img_flipud(varargin)

    % grab the data
    udat= get(gcf, 'userdata');

    % get the value of the slider
    sliceNum = round(get(udat.h.slider, 'Value'));

    % flip the data
    for a = 1:3;
        udat.stack.img{sliceNum}(:,:,a) = flipud(udat.stack.img{sliceNum}(:,:,a));
    end
    
    % plot the new image
    set(udat.h.img, 'CData', udat.stack.img{sliceNum})
    drawnow
    
    % save the new data
    set(udat.h.fig, 'userdata', udat);

end

function img_fliplr(varargin)

    % grab the data
    udat= get(gcf, 'userdata');

    % get the value of the slider
    sliceNum = round(get(udat.h.slider, 'Value'));

    % flip the data
    for a = 1:3;
        udat.stack.img{sliceNum}(:,:,a) = fliplr(udat.stack.img{sliceNum}(:,:,a));
    end
    
    % plot the new image
    set(udat.h.img, 'CData', udat.stack.img{sliceNum})
    drawnow
    
    % save the new data
    set(udat.h.fig, 'userdata', udat);


end

function img_change(varargin)
    
    % grab the data
    udat= get(gcf, 'userdata');
    
    % figure out the value of the slider
    sliceNum = round(get(udat.h.slider, 'Value'));
    
    
    % update the text field of the textbox
    set(udat.h.textbox, 'String', sprintf('%d of %d', sliceNum, numel(udat.stack.img)));
    
    % update the image
    set(udat.h.img, 'CData', udat.stack.img{sliceNum})
    drawnow
    
    % update the stem plot
    hand = get(udat.h.stem, 'children');
    hgt = ones(numel(udat.stack.loc),1);
    hgt(sliceNum) = 2;
    set(hand, 'YData', hgt);
    set(udat.h.stem, 'box', 'off');
    drawnow

end

function txt_update(varargin)
    

    % grab the data
    udat= get(gcf, 'userdata');
    
    % grab the value indicated in the text box
    txtString = get(udat.h.textbox, 'String');
    txtString = regexp(txtString, '[\d]+', 'match'); % find the number in the text string
    txtString = txtString{1}; % only consider the first (b/c the second entry could be the total number of slices)
    sliceNum = str2double(txtString);
    
    % do some error checking
    if sliceNum < 0 || sliceNum > numel(udat.stack.img)
        oldSliceNum = get(udat.h.slider, 'value');
        set(udat.h.textbox, 'String', sprintf('%d of %d', oldSliceNum, numel(udat.stack.img)));
        return
    end
    
    % update the slider position and the text of the text box
    set(udat.h.slider, 'Value', sliceNum)
    set(udat.h.textbox, 'String', sprintf('%d of %d', sliceNum, numel(udat.stack.img)));
    
    % plot the new image
    set(udat.h.img, 'CData', udat.stack.img{sliceNum})
    drawnow
    
    % update the stem plot
    hand = get(udat.h.stem, 'children');
    hgt = ones(numel(udat.stack.loc),1);
    hgt(sliceNum) = 2;
    set(hand, 'YData', hgt);
    set(udat.h.stem, 'box', 'off');
    drawnow
    

end

