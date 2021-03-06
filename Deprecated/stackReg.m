
function stackReg(mouseName,scalebar)
%mouseName = 'EB_031714_D';

if ~exist('scalebar', 'var')
    scalebar = 500;
end

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

% base = fullfile('Z:\home\charlie\Data\Mice', mouse, '\Histology');

% nfiles = size(files,1);
% target = imread(fullfile(base,['stack_' mouse '_ver1'],files(3,:).name));


h.fig = figure;
set(h.fig, 'position', [351 87 721 719])
h.rotate = uicontrol('style', 'edit',...
                     'units', 'normalized',...
                     'position', [0.15, 0.875, .1, .05],...
                     'String', 'Rotation (deg)',...
                     'Callback', {@img_rotate});
h.slider = uicontrol('style', 'slider',...
                     'units', 'normalized',...
                     'position', [0.7, 0.08, 0.2, 0.05],...
                     'Callback', {@img_change},...
                     'Value', 1,...
                     'max', numel(stack.img),...
                     'min', 1,...
                     'SliderStep', [1./numel(stack.img), 1./numel(stack.img).*10]);
h.save = uicontrol('style', 'togglebutton',...
                     'units', 'normalized',...
                     'Position', [0.1, 0.05, 0.25, 0.10],...
                     'String', 'Save as .tif',...
                     'Callback', {@img_save});
h.textbox = uicontrol('style', 'edit',...
                      'units', 'normalized',...
                      'string', sprintf('1 of %d', numel(stack.img)),...
                      'position', [0.73, 0.02, 0.15, 0.055],...
                      'callback', {@txt_update});
h.croptop = uicontrol('style', 'togglebutton',...
                     'units', 'normalized',...
                     'Position', [.27, .9, .1, .05],...
                     'String', 'Crop Top',...
                     'Callback', {@img_croptop});
h.cropbot = uicontrol('style', 'togglebutton',...
                     'units', 'normalized',...
                     'Position', [.27, .85,  .1, .05],...
                     'String', 'Crop Bottom',...
                     'Callback', {@img_cropbot}); 
h.cropleft = uicontrol('style', 'togglebutton',...
                     'units', 'normalized',...
                     'Position', [.39, .9,  .1, .05],...
                     'String', 'Crop Left',...
                     'Callback', {@img_cropleft}); 
h.cropright = uicontrol('style', 'togglebutton',...
                     'units', 'normalized',...
                     'Position', [.39, .85,  .1, .05],...
                     'String', 'Crop Right',...
                     'Callback', {@img_cropright});
h.cropsize = uicontrol('style', 'edit',...
                     'units', 'normalized',...
                     'Position', [.51, .875,  .1, .05],...
                     'String', 'Crop size',...
                     'Callback', {@img_cropsize}); 
h.sizerep = uicontrol('style', 'edit',...
                      'units', 'normalized',...
                      'string', sprintf('0 pix for y; 0 pix for x'),...
                      'position', [.27, .825,  .2, .025],...
                      'callback', {@img_sizereport});
h.shiftup = uicontrol('style', 'togglebutton',...
                     'units', 'normalized',...
                     'Position', [.62, .9, .1, .05],...
                     'String', 'Shift Up',...
                     'Callback', {@img_shiftup});
h.shiftdown = uicontrol('style', 'togglebutton',...
                     'units', 'normalized',...
                     'Position', [.62, .85,  .1, .05],...
                     'String', 'Shift Down',...
                     'Callback', {@img_shiftdown}); 
h.shiftleft = uicontrol('style', 'togglebutton',...
                     'units', 'normalized',...
                     'Position', [.74, .9,  .1, .05],...
                     'String', 'Shift Left',...
                     'Callback', {@img_shiftleft}); 
h.shiftright = uicontrol('style', 'togglebutton',...
                     'units', 'normalized',...
                     'Position', [.74, .85,  .1, .05],...
                     'String', 'Shift Right',...
                     'Callback', {@img_shiftright});
h.shiftsize = uicontrol('style', 'edit',...
                     'units', 'normalized',...
                     'Position', [.86, .875,  .1, .05],...
                     'String', 'Shift size',...
                     'Callback', {@img_shiftsize}); 
h.overlay = uicontrol('style', 'togglebutton',...
                     'units', 'normalized',...
                     'Position', [0.45, 0.08, 0.2, 0.05],...
                     'String', 'Overlay',...
                     'Callback', {@img_overlay});
                 
                
% h.register = uicontrol('style', 'togglebutton',...
%                      'units', 'normalized',...
%                      'Position', [0.6, 0.85, 0.25, 0.10],...
%                      'String', 'Register slices to target',...
%                      'Callback', {@img_register});
                     

% plot the first image to the screen. Shove the data into the 'userdata'
% field
h.img = plotimg(stack.img{1}, stack.info{1});
drawnow

udat.h = h; 
udat.stack = stack;
set(h.fig, 'UserData', udat)

end

function img_cropsize(varargin)

    % grab the data
    udat= get(gcf, 'userdata');
    
    % figure out the value of the slider
    sliceNum = round(get(udat.h.slider, 'Value'));
    
    % determine crop size
    crop = str2num(get(udat.h.cropsize, 'String'));
    
end

function img_shiftsize(varargin)

    % grab the data
    udat= get(gcf, 'userdata');
    
    % figure out the value of the slider
    sliceNum = round(get(udat.h.slider, 'Value'));
    
    % determine crop size
    shift = str2num(get(udat.h.shiftsize, 'String'));
    
end

function img_overlay(varargin)

    % grab the data
    udat= get(gcf, 'userdata');
    
    % figure out the value of number slider
    sliceNum = round(get(udat.h.slider, 'Value'));
    
    %make overlay
    falseColorOverlay = imfuse(udat.stack.img{sliceNum}, udat.stack.img{sliceNum-1});
    
    % update the image
    set(udat.h.img, 'CData', falseColorOverlay)
    drawnow
    
    %reset button
    h.redraw = uicontrol('style', 'togglebutton',...
                     'units', 'normalized',...
                     'Position', [0.45, 0.08, 0.2, 0.05],...
                     'String', 'Remove slice',...
                     'Callback', {@img_redraw});
end
function img_redraw(varargin)

    % grab the data
    udat= get(gcf, 'userdata');
    
    % figure out the value of number slider
    sliceNum = round(get(udat.h.slider, 'Value'));
    
    % update the image
    set(udat.h.img, 'CData', udat.stack.img{sliceNum})
    drawnow
        
    %reset button
    h.overlay = uicontrol('style', 'togglebutton',...
                     'units', 'normalized',...
                     'Position', [0.45, 0.08, 0.2, 0.05],...
                     'String', 'Overlay',...
                     'Callback', {@img_overlay});    
end

function img_sizereport(varargin)

    % grab the data
    udat= get(gcf, 'userdata');
    
    % figure out the value of the slider
    sliceNum = round(get(udat.h.slider, 'Value'));
    
    % determine target size
    siz = [udat.stack.info{1}.Height udat.stack.info{1}.Width];

    % determine current size
    siz_rot = size(udat.stack.img{sliceNum});
    
    %determine difference
    y_diff = siz_rot(1)-siz(1);
    x_diff = siz_rot(2)-siz(2);
    
    % update the text field of the textbox
    set(udat.h.sizerep, 'String', sprintf('%d pix for y; %d pix for x', y_diff, x_diff));
end

function img_shiftup(varargin)

    % grab the data
    udat= get(gcf, 'userdata');
    
    % figure out the value of the slider
    sliceNum = round(get(udat.h.slider, 'Value'));
    
    % determine cropsize
    shift = str2num(get(udat.h.shiftsize, 'String'));
    siz = size(udat.stack.img{sliceNum});
    
    %shift
    udat.stack.img{sliceNum}= [udat.stack.img{sliceNum}; zeros(shift,siz(2),siz(3))];
    udat.stack.img{sliceNum}(1:shift,:,:) = [];
    
    % plot the new image
    set(udat.h.img, 'CData', udat.stack.img{sliceNum})
    drawnow
    
    % save the new data
    set(udat.h.fig, 'userdata', udat);
end

function img_shiftdown(varargin)

    % grab the data
    udat= get(gcf, 'userdata');
    
    % figure out the value of the slider
    sliceNum = round(get(udat.h.slider, 'Value'));
    
    % determine cropsize
    shift = str2num(get(udat.h.shiftsize, 'String'));
    siz = size(udat.stack.img{sliceNum});
    
    %shift
    udat.stack.img{sliceNum}(siz(1)+1-shift:siz(1),:,:) = [];
    udat.stack.img{sliceNum}= [zeros(shift,siz(2),siz(3)); udat.stack.img{sliceNum}];
    
    % plot the new image
    set(udat.h.img, 'CData', udat.stack.img{sliceNum})
    drawnow
    
    % save the new data
    set(udat.h.fig, 'userdata', udat);
end

function img_shiftleft(varargin)

    % grab the data
    udat= get(gcf, 'userdata');
    
    % figure out the value of the slider
    sliceNum = round(get(udat.h.slider, 'Value'));
    
    % determine cropsize
    shift = str2num(get(udat.h.shiftsize, 'String'));
    siz = size(udat.stack.img{sliceNum});
    
    %shift
    udat.stack.img{sliceNum}= [udat.stack.img{sliceNum} zeros(siz(1),shift,siz(3))];
    udat.stack.img{sliceNum}(:,1:shift,:) = [];
    
    % plot the new image
    set(udat.h.img, 'CData', udat.stack.img{sliceNum})
    drawnow
    
    % save the new data
    set(udat.h.fig, 'userdata', udat);
end

function img_shiftright(varargin)

    % grab the data
    udat= get(gcf, 'userdata');
    
    % figure out the value of the slider
    sliceNum = round(get(udat.h.slider, 'Value'));
    
    % determine cropsize
    shift = str2num(get(udat.h.shiftsize, 'String'));
    siz = size(udat.stack.img{sliceNum});
    
    %shift
    udat.stack.img{sliceNum}(:,siz(2)+1-shift:siz(2),:) = [];
    udat.stack.img{sliceNum}= [zeros(siz(1),shift,siz(3)) udat.stack.img{sliceNum}];
    
    % plot the new image
    set(udat.h.img, 'CData', udat.stack.img{sliceNum})
    drawnow
    
    % save the new data
    set(udat.h.fig, 'userdata', udat);
end

function img_croptop(varargin)

    % grab the data
    udat= get(gcf, 'userdata');
    
    % figure out the value of the slider
    sliceNum = round(get(udat.h.slider, 'Value'));
    
    % determine cropsize
    crop = str2num(get(udat.h.cropsize, 'String'));
    
    %limit size
    siz = [udat.stack.info{1}.Height udat.stack.info{1}.Width];
    siz_rot = size(udat.stack.img{sliceNum});
    y_diff = siz_rot(1)-siz(1);
    if crop>y_diff
        crop = y_diff;
    end
    
    %crop
    udat.stack.img{sliceNum}(1:crop,:,:) = [];
    
    %pad if negative
    if y_diff < 0
        pad = -1*y_diff;
        udat.stack.img{sliceNum} = [zeros(pad,siz(2),siz(3)); udat.stack.img{sliceNum}];
    end
    
    % determine current size
    siz_rot = size(udat.stack.img{sliceNum});
    
    %determine difference
    y_diff = siz_rot(1)-siz(1);
    x_diff = siz_rot(2)-siz(2);
    
    % update the text field of the textbox
    set(udat.h.sizerep, 'String', sprintf('%d pix for y; %d pix for x', y_diff, x_diff));
    
    % plot the new image
    set(udat.h.img, 'CData', udat.stack.img{sliceNum})
    drawnow
    
    % save the new data
    set(udat.h.fig, 'userdata', udat);
end

function img_cropbot(varargin)

    % grab the data
    udat= get(gcf, 'userdata');
    
    % figure out the value of the slider
    sliceNum = round(get(udat.h.slider, 'Value'));
    
    % determine cropsize
    crop = str2num(get(udat.h.cropsize, 'String'));
    siz_rot = size(udat.stack.img{sliceNum});
    
    %limit size
    siz = [udat.stack.info{1}.Height udat.stack.info{1}.Width];
    siz_rot = size(udat.stack.img{sliceNum});
    y_diff = siz_rot(1)-siz(1);
    if crop>y_diff
        crop = y_diff;
    end
    
    %crop
    udat.stack.img{sliceNum}(1+siz_rot(1)-crop:siz_rot(1),:,:) = [];
    
    %pad if negative
    if y_diff < 0
        pad = -1*y_diff;
        udat.stack.img{sliceNum} = [udat.stack.img{sliceNum}; zeros(pad,siz(2),siz(3))];
    end

    % determine current size
    siz_rot = size(udat.stack.img{sliceNum});
    
    %determine difference
    y_diff = siz_rot(1)-siz(1);
    x_diff = siz_rot(2)-siz(2);
    
    % update the text field of the textbox
    set(udat.h.sizerep, 'String', sprintf('%d pix for y; %d pix for x', y_diff, x_diff));
    
    % plot the new image
    set(udat.h.img, 'CData', udat.stack.img{sliceNum})
    drawnow
    
    % save the new data
    set(udat.h.fig, 'userdata', udat);
end

function img_cropleft(varargin)

    % grab the data
    udat= get(gcf, 'userdata');
    
    % figure out the value of the slider
    sliceNum = round(get(udat.h.slider, 'Value'));
    
    % determine cropsize
    crop = str2num(get(udat.h.cropsize, 'String'));
    
    %limit size
    siz = [udat.stack.info{1}.Height udat.stack.info{1}.Width];
    siz_rot = size(udat.stack.img{sliceNum});
    x_diff = siz_rot(2)-siz(2);
    if crop>x_diff
        crop = x_diff;
    end
    
    %crop
    udat.stack.img{sliceNum}(:,1:crop,:) = [];
    
    %pad if negative
    if x_diff < 0
        pad = -1*x_diff;
        udat.stack.img{sliceNum} = [zeros(siz(1),pad,siz(3)) udat.stack.img{sliceNum}];
    end 

    % determine current size
    siz_rot = size(udat.stack.img{sliceNum});
    
    %determine difference
    y_diff = siz_rot(1)-siz(1);
    x_diff = siz_rot(2)-siz(2);

    % update the text field of the textbox
    set(udat.h.sizerep, 'String', sprintf('%d pix for y; %d pix for x', y_diff, x_diff));
    
    % plot the new image
    set(udat.h.img, 'CData', udat.stack.img{sliceNum})
    drawnow
    
    % save the new data
    set(udat.h.fig, 'userdata', udat);
end

function img_cropright(varargin)

    % grab the data
    udat= get(gcf, 'userdata');
    
    % figure out the value of the slider
    sliceNum = round(get(udat.h.slider, 'Value'));
    
    % determine cropsize
    crop = str2num(get(udat.h.cropsize, 'String'));
    siz_rot = size(udat.stack.img{sliceNum});
    
    %limit size
    siz = [udat.stack.info{1}.Height udat.stack.info{1}.Width];
    siz_rot = size(udat.stack.img{sliceNum});
    x_diff = siz_rot(2)-siz(2);
    if crop>x_diff
        crop = x_diff;
    end
    
    %crop
    udat.stack.img{sliceNum}(:,1+siz_rot(2)-crop:siz_rot(2),:) = [];
    
    %pad if negative
    if x_diff < 0
        pad = -1*x_diff;
        udat.stack.img{sliceNum} = [udat.stack.img{sliceNum} zeros(siz(1),pad,siz(3))];
    end 

    % determine current size
    siz_rot = size(udat.stack.img{sliceNum});
    
    %determine difference
    y_diff = siz_rot(1)-siz(1);
    x_diff = siz_rot(2)-siz(2);
    
    % update the text field of the textbox
    set(udat.h.sizerep, 'String', sprintf('%d pix for y; %d pix for x', y_diff, x_diff));
    
    % plot the new image
    set(udat.h.img, 'CData', udat.stack.img{sliceNum})
    drawnow
    
    % save the new data
    set(udat.h.fig, 'userdata', udat);
end

function img_rotate(varargin)

    % grab the data
    udat= get(gcf, 'userdata');
    
    % figure out the value of the slider
    sliceNum = round(get(udat.h.slider, 'Value'));
    
    % determine rotation
    udat.h.rot = str2num(get(udat.h.rotate, 'String'));
    
     % rotate the data
     udat.stack.img{sliceNum} = imrotate(udat.stack.img{sliceNum},udat.h.rot);
    
    % determine target size
    siz = [udat.stack.info{1}.Height udat.stack.info{1}.Width];

    % determine current size
    siz_rot = size(udat.stack.img{sliceNum});
    
    %determine difference
    y_diff = siz_rot(1)-siz(1);
    x_diff = siz_rot(2)-siz(2);
    
    % update the text field of the textbox
    set(udat.h.sizerep, 'String', sprintf('%d pix for y; %d pix for x', y_diff, x_diff));
    
    % plot the new image
    set(udat.h.img, 'CData', udat.stack.img{sliceNum})
    drawnow
    
    % save the new data
    set(udat.h.fig, 'userdata', udat);

end

function img_save(varargin)
    % grab the data
    udat= get(gcf, 'userdata');
    
    % figure out the value of the slider
    sliceNum = round(get(udat.h.slider, 'Value'));
    
    % save as tif
    
    contents = dir(fileparts(udat.stack.info{sliceNum}.Filename));
    slice_name = contents(sliceNum+2).name;
    mouseName = slice_name(1:11);
    file_name = [slice_name(1:length(slice_name)-4) '_reg.tif'];
    
    if ~exist(fullfile(fileparts(cd),['stack_' mouseName '_reg']))
         mkdir(fullfile(fileparts(cd),['stack_' mouseName '_reg']));
    end
    
    imwrite(udat.stack.img{sliceNum}, fullfile(fileparts(cd),['stack_' mouseName '_reg'],file_name),'TIFF');    
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
end

       
function stack = compileStack(stackDir, mouseName, scalebar);
    
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

% function img_register(varargin)
%     % grab the data
%     udat= get(gcf, 'userdata');
%     
%     % figure out the value of the slider
%     sliceNum = round(get(udat.h.slider, 'Value'));
%     
%     % save first slice
%     if sliceNum == 2
%         % create save path
%         contents = dir(fileparts(udat.stack.info{1}.Filename));
%         slice_name = contents(3).name;
%         mouseName = slice_name(1:11);
%         file_name = [slice_name(1:length(slice_name)-4) '_reg.tif'];
%         udat.stack.img_reg{1} = udat.stack.img_rot{1}; 
%         if ~exist(fullfile(fileparts(cd),['stack_' mouseName '_reg']))
%             mkdir(fullfile(fileparts(cd),['stack_' mouseName '_reg']));
%         end
%         imwrite(udat.stack.img_reg{1}, fullfile(fileparts(cd),['stack_' mouseName '_reg'],file_name),'TIFF');
%     end
%     
%     target = udat.stack.img_reg{sliceNum-1};
%     siz = size(target);
% 
%     img = udat.stack.img_rot{sliceNum};
%     [input_points, base_points] = cpselect(img,target,'Wait', true);         
%     tform= fitgeotrans(input_points, base_points, 'similarity');
%     [registered RA] = imwarp(img,tform);
%     Tinv  = tform.invert.T;
%     ss = Tinv(2,1);
%     sc = Tinv(1,1);
%     scaleRecovered = sqrt(ss*ss + sc*sc);
%     reg_scale = imresize(registered, scaleRecovered);
%         
%     y_offset = ceil(RA.YWorldLimits(1)*scaleRecovered);
%     x_offset = ceil(RA.XWorldLimits(1)*scaleRecovered);
%     if y_offset < 0
%         reg_scale_buffer = reg_scale;
%         y_range = 1-y_offset:siz(1)-y_offset;
%     else
%         reg_scale_buffer = [zeros(y_offset, size(reg_scale,2),size(reg_scale,3)); reg_scale];
%         y_range = 1:siz(1);
%     end
%     if x_offset < 0
%         reg_scale_buffer = reg_scale;   
%         x_range = 1-x_offset:siz(2)-x_offset;
%     else
%         reg_scale_buffer = [zeros(size(reg_scale,1), x_offset,size(reg_scale,3)) reg_scale];
%         x_range = 1:siz(2);
%     end
%     reg = reg_scale(y_range, x_range,:);
%     falsecolorOverlay = imfuse(target,reg);
% 
%     figure; subplot(2,2,1); 
%     imagesq(img); title('img'); 
%     subplot(2,2,2); 
%     imagesq(target); title('target'); 
%     subplot(2,2,3); 
%     imagesq(reg); title('reg'); 
%     subplot(2,2,4); 
%     imagesq(falsecolorOverlay); title('overlay');
%         
%     % create save path
%     contents = dir(fileparts(udat.stack.info{sliceNum}.Filename));
%     slice_name = contents(sliceNum+2).name;
%     mouseName = slice_name(1:11);
%     file_name = [slice_name(1:length(slice_name)-4) '_reg.tif'];
%     
%     udat.stack.img_reg{sliceNum} = reg;
%     imwrite(reg, fullfile(fileparts(cd),['stack_' mouseName '_reg'],file_name),'TIFF');
% end





