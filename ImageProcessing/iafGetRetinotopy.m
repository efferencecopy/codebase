function iafGetRetinotopy(fnames, relevantTrialTypes)


%
%  COMPUTE THE RETINOTOPY
%

% extract the runs, separate the data by trial type
dat = iafExtractRuns(fnames, relevantTrialTypes);
N_types = size(dat.uniqueTrlTypes,1);

% compute the average movie for each trial type
udat.AvgTrialON = {};
udat.AvgTrialOFF = {};
for i_t_types = 1 : N_types    
    udat.AvgTrialON{i_t_types} = iafGetAvgMovie(dat.OnResponse{i_t_types});
    udat.AvgTrialOFF{i_t_types} = iafGetAvgMovie(dat.OffResponse{i_t_types});
end


% make a figure showing the final images.
udat.ttypes = dat.uniqueTrlTypes;
udat.text = dat.trlTypeHeader;
udat.h.fig = figure;
udat.h.fig.Units = 'Normalized';
udat.h.fig.Position = [ 0.2576    0.1867    0.3951    0.6389];
udat.h.ax = axes('position', [0.15 0.28 0.70 0.70]);
udat.h.slide = uicontrol('style', 'slider', 'units', 'normalized', 'position', [0.30, 0.1, 0.40, 0.05],...
                     'Min', 1, 'Max', N_types, 'sliderstep', ones(1,2).*(1/(N_types-1)),...
                     'Value', 1, 'Interruptible', 'on', 'Callback', @updateImage, 'Enable', 'on',...
                     'ButtonDownFcn', @updateImage);
udat.h.artBox = uicontrol('style', 'edit', 'units', 'normalized', 'position', [0.25, 0.02, 0.1, 0.05],...
                          'string', 5, 'Callback', @img_preProcess);
uicontrol('style', 'text', 'units', 'normalized', 'position', [0.14, 0.014, 0.1, 0.05], ...
          'string', 'Artifact', 'FontSize', 12)
udat.h.smoothBox = uicontrol('style', 'edit', 'units', 'normalized', 'position', [0.75, 0.02, 0.1, 0.05],...
                             'string', 1, 'Callback', @img_preProcess);
uicontrol('style', 'text', 'units', 'normalized', 'position', [0.595, 0.014, 0.14, 0.05], ...
          'string', 'Smoothing', 'FontSize', 12)
udat.h.fig.UserData = udat;
img_preProcess();

end



function updateImage(one, two)
    
    udat = get(gcf, 'UserData');
    imgNum = udat.h.slide.Value;
    imgNum = max([1, imgNum]); % slider step requires range, but don't let the img num < 0
    imgNum = round(imgNum);
    axes(udat.h.ax)
    imagesc(udat.final_img{imgNum}); colormap('gray')
    udat.h.ax.XTick = [];
    udat.h.ax.YTick = [];
    
    % specify the stimulus position
    elidx = strcmpi(udat.text, 'tGratingElevationDeg');
    azidx = strcmpi(udat.text, 'tGratingAzimuthDeg');
    el = udat.ttypes(imgNum, elidx);
    az = udat.ttypes(imgNum, azidx);
    udat.h.ax.XLabel.String = sprintf('Az: %.0f, El: %.0f', az, el);
    
    set(gcf, 'UserData', udat);
end

function final_img = img_preProcess(one, two) %AvgTrialON, AvgTrialOFF, filterSD, artifactMultiplier
    
    % grab the udat
    udat = get(gcf, 'userdata');
    
    % figure out the status of all the ui controls
    filterSD = round(str2double(udat.h.smoothBox.String));
    artifactMultiplier = round(str2double(udat.h.artBox.String));
    
    % look for oob errors on the inputs and fix
    if filterSD <= 0
        filterSD = 0;
        udat.h.smoothBox.String = 'off';
    elseif filterSD > 5;
        filterSD = 5;
        udat.h.smoothBox.String = '5';
    end
    
    if artifactMultiplier <= 0
        artifactMultiplier = 0;
        udat.h.artBox.String = 'off';
    end
    
    
    % now do the pre-processing
    N_types = numel(udat.AvgTrialON);
    udat.final_img = {};
    for i_type = 1:N_types
        
        tmpOn = udat.AvgTrialON{i_type};
        tmpOff = udat.AvgTrialOFF{i_type};
        
        if artifactMultiplier > 0
            tmpOn = iafRemoveOutliers(tmpOn, artifactMultiplier);
            tmpOff = iafRemoveOutliers(tmpOff, artifactMultiplier);
        end
        
        if filterSD > 0 % could be different frame nums, so need sequential loops
            for i_frame = 1:size(tmpOn,3);
                tmpOn(:,:,i_frame) = imgaussfilt(tmpOn(:,:,i_frame), filterSD);
            end
            for i_frame = 1:size(tmpOff,3);
                tmpOff(:,:,i_frame) = imgaussfilt(tmpOff(:,:,i_frame), filterSD);
            end
        end
        
        meanAcrossTime_on = mean(tmpOn, 3);
        meanAcrossTime_off = mean(tmpOff, 3);
        
        udat.final_img{i_type} = meanAcrossTime_on - meanAcrossTime_off;
    end
    
    udat.h.fig.UserData = udat;
    updateImage()
end





