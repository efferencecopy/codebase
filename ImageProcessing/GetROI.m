function GetROI(udat)

% clear the command window so that the user can see the feedback. Set the
% 'help' variable so that the help turns off after the first click
udat.h.HELPON = true;
clc

% Specify which brain areas should be considered, and initalize the ROIs
udat.ROI.VisArea = {'A'; 'AL'; 'AM'; 'LM'; 'PM'; 'RL'; 'V1'};
Nareas = numel(udat.ROI.VisArea); 
udat.ROI.RCidx = repmat({[nan,nan]}, 1, Nareas);

% Code to initialize figure
udat.h.fig = figure;
udat.h.fig.Units = 'Normalized';
udat.h.fig.Position = [0.4899   0.2450    0.3211    0.6188];   
udat.h.ax = axes('position', [0.15   0.35   0.70   0.60]);

% Slider
N_types = numel(udat.AvgTrialON);
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

% Radio Button group
udat.h.rButtons = uibuttongroup('Visible', 'on',... 
                                'units', 'Normalized',...
                                'title', 'Use Average Image?',...
                                'position', [0.30   0.155   0.40   0.07],...
                                'SelectionChangedFcn', @rBtnSelection);
udat.h.rBtn1 = uicontrol(udat.h.rButtons,...
                         'Style', 'radiobutton',...
                         'units', 'Normalized',...
                         'String','Yes.',...
                         'Position', [0.17   0.3   0.3   0.6],...
                         'Tag', 'yes');
udat.h.rBtn2 = uicontrol(udat.h.rButtons,...
                         'Style', 'radiobutton',...
                         'units', 'Normalized',...
                         'String','No.',...
                         'Position', [0.63   0.3   0.3   0.6],...
                         'Tag', 'no');
udat.h.rButtons.SelectedObject = udat.h.rBtn2; % Set the Radio Button #2 as default ('No')

                     
% Pop-up Menu (for specifying the 'selected' Visual Area)
udat.h.menu = uicontrol('style', 'popupmenu',... 
                        'String',cat(1,udat.ROI.VisArea, '-Select One-'),...                        
                        'units', 'Normalized',...
                        'position', [0.375   0.115   0.25   0.03],...
                        'Value', 8);
                                          
% Button for 'select ROI'
udat.h.pButtonROI = uicontrol('style', 'pushbutton',...
                              'units', 'Normalized',...
                              'String', 'Select R.O.I.',...
                              'Position', [0.3786   0.03   0.25   0.05],...
                              'CallBack', @SelectROI);
% Button for Export data
udat.h.pButtonE = uicontrol('style', 'pushbutton',...
                            'units', 'Normalized',...
                            'String', 'EXPORT',...
                            'FontSize', 7.5,...
                            'FontWeight', 'bold',...
                            'Position', [0.823   0.027   0.15   0.045],...
                            'CallBack', @DataExport);


udat.h.fig.UserData = udat;
updateImage();
    
end


function updateImage(~,~)

    udat = get(gcf, 'UserData');
    
    % generate the mask of ROIs
    maskImg = zeros(size(udat.final_img{1}));
    for i_roi = 1:numel(udat.ROI.RCidx)
        rc = udat.ROI.RCidx{i_roi};
        if ~isempty(rc) && ~any(isnan(rc(:)))
            linind = sub2ind(size(maskImg), rc(:,1), rc(:,2));
            maskImg(linind) = 1;
        end
    end
    roi_linind = maskImg(:) == 1;


    switch udat.h.slide.Enable
        case 'on' % using the slider to pick a 'final img'
            imgNum = udat.h.slide.Value;
            imgNum = max([1, imgNum]); % slider step requires range, but don't let the img num < 0
            imgNum = round(imgNum);
            axes(udat.h.ax);
            
            pltimg = udat.final_img{imgNum};
            pltimg = pltimg + (abs(min(pltimg(:))));
            pltimg = pltimg ./ (max(pltimg(:)).*1.1);
            pltimg = repmat(pltimg, [1,1,3]);
            pltimg(roi_linind) = 0;

            
            imshow(pltimg); colormap('gray');
            udat.h.ax.XTick = [];
            udat.h.ax.YTick = [];
            
            % specify the stimulus position
            elidx = strcmpi(udat.text, 'tGratingElevationDeg');
            azidx = strcmpi(udat.text, 'tGratingAzimuthDeg');
            el = udat.ttypes(imgNum, elidx);
            az = udat.ttypes(imgNum, azidx);
            udat.h.ax.XLabel.String = sprintf('Az: %.0f, El: %.0f', az, el);
            
            
        case 'off'
            
            axes(udat.h.ax); % set the current axis
            
            pltimg = mean(cat(3, udat.final_img{:}),3);
            pltimg = pltimg + (abs(min(pltimg(:))));
            pltimg = pltimg ./ (max(pltimg(:)).*1.2);
            pltimg = repmat(pltimg, [1,1,3]);
            pltimg(roi_linind) = 0;

            imshow(pltimg); colormap('gray');
            
            udat.h.ax.XTick = [];
            udat.h.ax.YTick = [];
            udat.h.ax.XLabel.String = sprintf('Average Image');
    end
    
    set(gcf, 'UserData', udat);
    
end

function rBtnSelection(~, callbackdata)
    
    % grab the user data
    udat = get(gcf, 'UserData');

    if strcmpi(callbackdata.NewValue.Tag,  'yes'); % use the average image
        
        % Turn off Slider-- there is only one image to display
        udat.h.slide.Enable = 'Off';
        
    else strcmpi(callbackdata.NewValue.Tag,  'no');
        
        % Turn on Slider-- there are multiple images to display
        udat.h.slide.Enable = 'On';
    end
    
    udat.h.fig.UserData = udat;
    updateImage();
end
 

 function SelectROI(~, ~)
 
    udat = get(gcf, 'UserData');
    
    % which brain area is about to be traced?
    brainAreaIdx = udat.h.menu.Value;
    brainAreaName = udat.h.menu.String{brainAreaIdx};
    
    % tell the user which brain area will be traced, and if an existing
    % ROIc entry will be overwritten
    if any(isnan(udat.ROI.RCidx{brainAreaIdx}))
        fprintf('  Drawing a new ROI for area: %s\n', brainAreaName);
    else
        fprintf('  Overwriting the ROI for area: %s\n', brainAreaName);
    end
    
    % only display the help dialog once per session
    if udat.h.HELPON
        message = sprintf('Left click and hold to begin drawing.\nSimply lift the mouse button to finish');
        uiwait(msgbox(message));
        udat.h.HELPON = false;
    end
    
    % draw the ROI and store the XY pairs of the ROI
    hFH = imfreehand('closed', 'true'); % Draw...
    binaryImage = hFH.createMask(); % Create a binary image ("mask") from the ROI object.
    delete(hFH);

    % Get coordinates of the points located in the freehand drawn region.
    [rows, columns] = find(binaryImage);
    ROIc = [rows, columns];

    % Place ROI Coordinates (ROIc) in the index for the specified VisArea
    udat.ROI.RCidx{brainAreaIdx} = ROIc;
    
    % set the user data field
    udat.h.fig.UserData = udat;
    updateImage();
    
end


function DataExport(~,~)

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

    save(savepath, 'udat')
    
end
