function imgAlign(stack, info)




% roll over the stacks and make sure the images are all in the correct
% orientation
h.fig = figure;
set(h.fig, 'position', [351 161 716 645])
h.fliplr = uicontrol('style', 'togglebutton',...
                     'units', 'normalized',...
                     'string', 'FLIP L/R',...
                     'Position', [0.1, 0.05, 0.25, 0.10],...
                     'Callback', {@img_fliplr, h.fig});
h.flipud = uicontrol('style', 'togglebutton',...
                     'units', 'normalized',...
                     'string', 'FLIP U/D',...
                     'Position', [0.4, 0.05, 0.25, 0.10],...
                     'Callback', {@img_flipud});
h.slider = uicontrol('style', 'slider',...
                     'units', 'normalized',...
                     'position', [0.7, 0.01, 0.2, 0.1],...
                     'max', numel(stack),...
                     'min', 1,...
                     'Callback', {@img_change});
h.textbox = uicontrol('style', 'edit',...
                      'units', 'normalized',...
                      'string', sprintf('1 of %d', numel(stack)),...
                      'position', [0.73, 0.007, 0.15, 0.055],...
                      'callback', {@txt_update});
                  


 % plot the firt image to the screen. Shove the data into the 'userdata'
 % field
 plotimg(stack{1}, info);
 dat.h = h;
 dat.stack = stack;
 dat.info = info;
 set(h.fig, 'UserData', dat)

end



function img_flipud(varargin)

dat = get(gcf, 'userdata');

end



function img_fliplr(varargin)

% grab the data
dat = get(gcf, 'userdata');




end

function img_change(varargin)

dat = get(gcf, 'userdata');
end

