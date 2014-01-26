function stackViewer(stack, info)

% roll over the stacks and make sure the images are all in the correct
% orientation
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
                  


 % plot the firt image to the screen. Shove the data into the 'userdata'
 % field
 h.img = plotimg(stack.img{1}, info);
 drawnow
 
 % make a plot of all the slices (stem plot) and make the lollipop for the
 % current slice slightly taller than the rest.
 h.stem = axes;
 hgt = ones(numel(stack.loc),1);
 hgt(1) = 2;
 stem(stack.loc, hgt, 'k')
 set(h.stem, 'position', [.1, .87, .8, .1],...
             'box', 'off',...
             'ycolor', [1 1 1],...
             'ytick', [],...
             'ylim', [0, 2.3]);
 
 
 dat.h = h;
 dat.stack = stack;
 dat.info = info;
 set(h.fig, 'UserData', dat)

end



function img_flipud(varargin)

    % grab the data
    dat = get(gcf, 'userdata');

    % get the value of the slider
    sliceNum = round(get(dat.h.slider, 'Value'));

    % flip the data
    for a = 1:3;
        dat.stack.img{sliceNum}(:,:,a) = flipud(dat.stack.img{sliceNum}(:,:,a));
    end
    
    % plot the new image
    set(dat.h.img, 'CData', dat.stack.img{sliceNum})
    drawnow
    
    % save the new data
    set(dat.h.fig, 'userdata', dat);

end

function img_fliplr(varargin)

    % grab the data
    dat = get(gcf, 'userdata');

    % get the value of the slider
    sliceNum = round(get(dat.h.slider, 'Value'));

    % flip the data
    for a = 1:3;
        dat.stack.img{sliceNum}(:,:,a) = fliplr(dat.stack.img{sliceNum}(:,:,a));
    end
    
    % plot the new image
    set(dat.h.img, 'CData', dat.stack.img{sliceNum})
    drawnow
    
    % save the new data
    set(dat.h.fig, 'userdata', dat);


end

function img_change(varargin)
    
    % grab the data
    dat = get(gcf, 'userdata');
    
    % figure out the value of the slider
    sliceNum = round(get(dat.h.slider, 'Value'));
    
    
    % update the text field of the textbox
    set(dat.h.textbox, 'String', sprintf('%d of %d', sliceNum, numel(dat.stack.img)));
    
    % update the image
    set(dat.h.img, 'CData', dat.stack.img{sliceNum})
    drawnow
    
    % update the stem plot
    hand = get(dat.h.stem, 'children');
    hgt = ones(numel(dat.stack.loc),1);
    hgt(sliceNum) = 2;
    set(hand, 'YData', hgt);
    set(dat.h.stem, 'box', 'off');
    drawnow

end

function txt_update(varargin)
    

    % grab the data
    dat = get(gcf, 'userdata');
    
    % grab the value indicated in the text box
    txtString = get(dat.h.textbox, 'String');
    txtString = regexp(txtString, '[\d]+', 'match'); % find the number in the text string
    txtString = txtString{1}; % only consider the first (b/c the second entry could be the total number of slices)
    sliceNum = str2double(txtString);
    
    % do some error checking
    if sliceNum < 0 || sliceNum > numel(dat.stack.img)
        oldSliceNum = get(dat.h.slider, 'value');
        set(dat.h.textbox, 'String', sprintf('%d of %d', oldSliceNum, numel(dat.stack.img)));
        return
    end
    
    % update the slider position and the text of the text box
    set(dat.h.slider, 'Value', sliceNum)
    set(dat.h.textbox, 'String', sprintf('%d of %d', sliceNum, numel(dat.stack.img)));
    
    % plot the new image
    set(dat.h.img, 'CData', dat.stack.img{sliceNum})
    drawnow
    
    % update the stem plot
    hand = get(dat.h.stem, 'children');
    hgt = ones(numel(dat.stack.loc),1);
    hgt(sliceNum) = 2;
    set(hand, 'YData', hgt);
    set(dat.h.stem, 'box', 'off');
    drawnow
    

end

