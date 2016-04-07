function GetROI(udat)

N_types = numel(udat.AvgTrialON);

% Code to initialize figure
udat.h.fig = figure;
udat.h.fig.Units = 'Normalized';
udat.h.fig.Position = [0.4899   0.2450    0.3211    0.6188];   
udat.h.ax = axes('position', [0.15   0.35   0.70   0.60]);

% Slider
udat.h.slide = uicontrol('style', 'slider',... 
                         'units', 'Normalized',...
                         'position', [0.30   0.235   0.40   0.03],...
                         'Min', 1, 'Max', N_types,...
                         'sliderstep', ones(1,2).*(1/(N_types-1)),...
                         'Value', 1,... 
                         'Interruptible', 'on',...
                         'Callback', @updateImage,...
                         'Enable', 'on',...
                         'ButtonDownFcn', @updateImage);

% Radio Button
udat.h.rButtons = uibuttongroup('Visible', 'on',... 
                                'units', 'Normalized',...
                                'title', 'Use Average Image?',...
                                'position', [0.30   0.155   0.40   0.07],...
                                'SelectionChangedFcn', @rBtnSelection);
% Radio Button Options
udat.h.rBtn1 = uicontrol(udat.h.rButtons,...
                         'Style', 'radiobutton',...
                         'units', 'Normalized',...
                         'String','Yes.',...
                         'Position', [0.17   0.3   0.3   0.6],...
                         'Tag', '1');
udat.h.rBtn2 = uicontrol(udat.h.rButtons,...
                         'Style', 'radiobutton',...
                         'units', 'Normalized',...
                         'String','No.',...
                         'Position', [0.63   0.3   0.3   0.6],...
                         'Tag', '2');

% Pop-up Menu (for specifying the 'selected' Visual Area)
udat.h.menu = uicontrol('style', 'popupmenu',... 
                        'String',{'A';'AL';'AM';'LM';'PM';'RL';'V1'; '  -Select One-'},...                        
                        'units', 'Normalized',...
                        'position', [0.3983   0.115   0.216   0.03],...
                        'Callback', @SpecifyROI);
                                          
% Push Buttons: 
% Call ROI function
udat.h.pButtonROI = uicontrol('style', 'pushbutton',...
                              'units', 'Normalized',...
                              'String', 'Select R.O.I.',...
                              'Position', [0.3786   0.03   0.25   0.05],...
                              'CallBack', @SelectROI);
% Export data
udat.h.pButtonE = uicontrol('style', 'pushbutton',...
                            'units', 'Normalized',...
                            'String', 'EXPORT',...
                            'FontSize', 7.5,...
                            'FontWeight', 'bold',...
                            'Position', [0.823   0.027   0.15   0.045],...
                            'CallBack', @DataExport);
                       
udat.h.fig.UserData = udat;

set(udat.h.menu, 'Value', 8); % Set a default value ('5' being the 2nd on the pop-up list)
set(udat.h.rButtons, 'SelectedObject', udat.h.rBtn2); % Set the Radio Button #2 as default ('No')
fprintf('\n*NOTE: User viewing images in "Default Settings".\n')
    updateImage();
end


function updateImage(one, two)

    udat = get(gcf, 'UserData');

    imgNum = udat.h.slide.Value;
    imgNum = max([1, imgNum]); % slider step requires range, but don't let the img num < 0
    imgNum = round(imgNum);
    axes(udat.h.ax);
    
    imagesc(udat.final_img{imgNum}); colormap('gray');
    udat.h.ax.XTick = [];
    udat.h.ax.YTick = [];
    
    % specify the stimulus position
    elidx = strcmpi(udat.text, 'tGratingElevationDeg');
    azidx = strcmpi(udat.text, 'tGratingAzimuthDeg');
    el = udat.ttypes(imgNum, elidx);
    az = udat.ttypes(imgNum, azidx);
    udat.h.ax.XLabel.String = sprintf('Az: %.0f, El: %.0f', az, el);
    
    udat.sliderValue = imgNum; % 'Artificial' Slider Value for differentiating between images later
    set(gcf, 'UserData', udat);
end

function rBtnSelection(source, callbackdata)
udat = get(gcf, 'UserData');
  
if callbackdata.NewValue.Tag == '1';
    avgImage = {mean(cat(3, udat.final_img{:}),3)};
    
    % Turn off Slider-- there is only one image to display
    set(udat.h.slide, 'Enable', 'Off');
    udat.sliderValue = 0; % 'Artificial' Slider Value for differentiating between images later
    
    imagesc(avgImage{:}); colormap('gray');
    udat.h.ax.XTick = [];
    udat.h.ax.YTick = [];
    
    % Instead of Stimulus Position, display...
    udat.h.ax.XLabel.String = sprintf('Average Image');
       
    udat.ROI.AvgImage = {};
    set(gcf, 'UserData', udat);
     
else callbackdata.NewValue.Tag == '2';
    % Turn on Slider-- there are multiple images to display
    set(udat.h.slide, 'Enable', 'On');
    updateImage();
end
end


function SpecifyROI (hObject, eventdata, handles)
udat = get(gcf, 'UserData');
% Hints: contents = cellstr(get(hObject,'String')) returns contents...
%        contents{get(hObject,'Value')} returns selected item...
items = get(hObject,'String');
index_selected = get(hObject,'Value');
item_selected = items{index_selected};

% Visual Area names to correspond with the ROIc
udat.ROI.VisArea = {'A','AL','AM','LM','PM','RL','V1'};
udat.ROI.idxVisArea = index_selected;

set(gcf, 'UserData',udat);
end
function SelectROI (source, eventdata)
udat = get(gcf, 'UserData');

if udat.sliderValue == 0;
    fprintf('...Selecting R.O.I.s in Average Image.\n')
else udat.sliderValue <= 0;
    imgNum = udat.h.slide.Value;
    fprintf('...Selecting R.O.I.s in Image %d.\n', imgNum)
end

message = sprintf('Left click and hold to begin drawing.\nSimply lift the mouse button to finish');
uiwait(msgbox(message));
hFH = imfreehand(); % Draw...

binaryImage = hFH.createMask(); % Create a binary image ("mask") from the ROI object.

%%%%% Don't necessarily need to create figure each time you select an ROI
% Display the freehand mask.
% udat.H.fig = figure;
% udat.H.fig.Position = [1005  448  256    215];
% imagesc(binaryImage); colormap('gray'); axis off;
% title('Binary "Mask"', 'FontSize', 11);
%%%%%

% Calculate the area, in pixels, that they drew.
DrawnPixelArea_integer = sum(binaryImage(:));
% Another way to calculate it that takes fractional pixels into account.
DrawnPixelArea = bwarea(binaryImage);

% Get coordinates of the points located in the freehand drawn region.
[rows, columns] = find(binaryImage);
ROIc = [columns, rows];
    
% udat.ROI.idxVisArea; Place ROI Coordinates (ROIc) in the index for the specified VisArea
if udat.sliderValue == 0;
    udat.ROI.AvgImage{udat.ROI.idxVisArea} = ROIc;
else udat.sliderValue >= 0;
    imgNum = udat.h.slide.Value;
    udat.ROI.final_img{imgNum}{udat.ROI.idxVisArea} = ROIc;
end

udat.h.fig.UserData = udat;
end


function DataExport (source, eventdata)
udat = get(gcf, 'UserData');

musID = regexpi((regexpi(udat.fid.names_mat{1}, 'data-(\w\d+)', 'match')), '(\w\d+)', 'match');
    musID = musID{1}{1};
Date = regexpi((regexpi(udat.fid.names_mat{1}, '-(\d+)_','match')), '(\d+)', 'match');
    Date = Date{1}{1};
for i_names = 1: numel(udat.fid.names_mat)
Repeat{i_names} = regexpi(udat.fid.names_mat{i_names}, 'Run(\d+)', 'match');
NRepeat(i_names) = regexpi(Repeat{i_names}, '[\d]', 'match');
end

fname = [musID, '_', Date, '_Run[', NRepeat{1}{1}, '-', NRepeat{end}{1}, ']_processedROIs.mat'];

% Save udat file under the predetermined name
[fname, pathname] = uiputfile('.mat', 'Export all Data', fname);
savepath = [pathname, filesep, fname];
    % Decrease file size    
    udat.h = []; 
    udat.sliderValue = [];
    udat.ROI.idxVisArea = [];

save(savepath, 'udat')
end
