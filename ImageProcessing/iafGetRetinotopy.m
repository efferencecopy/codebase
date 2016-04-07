function iafGetRetinotopy(experimentType)


% specify the actual trial types from the experiment type
relevantTrialTypes = trialTypeLibrary(experimentType);
udat.fid.ttypes = relevantTrialTypes;
% has the user specified the file names? if not, call uigetfile and
% figure out which .mat files correspond to which .tiff files. do this in a
% sub-function
[udat.fid.names_mat, udat.fid.names_img, udat.fid.path] = ioGetFnames;


%
%  COMPUTE THE RETINOTOPY
%

% extract the runs, separate the data by trial type
cd(udat.fid.path)
dat = iafExtractRuns(udat.fid.names_mat, udat.fid.names_img, relevantTrialTypes);
N_types = size(dat.uniqueTrlTypes,1);

% save the ImagingRate in the output of this function so that the user can 
% analyze this data later
udat.frameRate = dat.frameRate;

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
udat.h.fig.Position = [ 0.1526    0.2450    0.3211    0.6188];
udat.h.ax = axes('position', [0.15 0.35 0.70 0.60]);
udat.h.slide = uicontrol('style', 'slider', 'units', 'normalized', 'position', [0.30, 0.235, 0.40, 0.03],...
                     'Min', 1, 'Max', N_types, 'sliderstep', ones(1,2).*(1/(N_types-1)),...
                     'Value', 1, 'Interruptible', 'on', 'Callback', @updateImage, 'Enable', 'on',...
                     'ButtonDownFcn', @updateImage);
udat.h.artBox = uicontrol('style', 'edit', 'units', 'normalized', 'position', [0.25, 0.05, 0.1, 0.05],...
                          'string', 5, 'Callback', @img_preProcess);
uicontrol('style', 'text', 'units', 'normalized', 'position', [0.12, 0.04, 0.12, 0.05], ...
          'string', 'Artifact', 'FontSize', 10)
udat.h.smoothBox = uicontrol('style', 'edit', 'units', 'normalized', 'position', [0.75, 0.05, 0.1, 0.05],...
                             'string', 1, 'Callback', @img_preProcess);
uicontrol('style', 'text', 'units', 'normalized', 'position', [0.573, 0.04, 0.17, 0.05], ...
          'string', 'Smoothing', 'FontSize', 10)
udat.h.PushBtn = uicontrol('style', 'pushbutton',...
                           'units', 'Normalized',...
                           'String', 'Get R.O.I.',...
                           'Position', [0.3786   0.15   0.25   0.05],...
                           'Callback', @get_ROI);
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

function img_preProcess(one, two) %AvgTrialON, AvgTrialOFF, filterSD, artifactMultiplier
    
    % grab the udat
    udat = get(gcf, 'UserData');
    
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
    
    % include the pre-processed data in the output of this function so that
    % the user can analyze these data later
    udat.preProcessed_ON{i_type} = tmpOn;
    udat.preProcessed_OFF{i_type} = tmpOff;
    end
    
    udat.h.fig.UserData = udat;
    updateImage()
end

function [names_mat, names_img, path] = ioGetFnames()

    % call uigetfile. with default directory. Save
    [names, path] = uigetfile({'*.mat;*.tif',...
                               'Related Files (*.mat;*.tif)'},...
                               'Select the Related MATLAB & TIFF files',...
                               'MultiSelect', 'on',...
                               'S:\Data');
    % Display choice
    if isequal(names,0)
        errordlg(sprintf(' File not found.\n User selected "Cancel".'), 'File Error')
    end
    
    names_mat = names(strmatch('data', names));
    names_img = names(strmatch('tiff', names));

    % Attempt to check that the MATLAB files match up with the TIFF files
    Nnames = numel(names_mat);
    for i_names = 1:Nnames
        tag_mat{i_names} = regexpi(names_mat{i_names}, 'run[\d].', 'match');
        tag_img{i_names} = regexpi(names_img{i_names}, 'run[\d].', 'match');
        if strcmpi(tag_mat{i_names}, tag_img{i_names}) == 1;
        else
            errordlg ('Your MATLAB files may not match up with your TIFF files!')
        end
    end
end

function relevantTrialTypes = trialTypeLibrary(experimentType)

 switch lower(experimentType)
     case 'retinotopy'
         relevantTrialTypes = {'tGratingElevationDeg', 'tGratingAzimuthDeg'};
     case 'sfvstf'
         relevantTrialTypes = {'tGratingSpatialFreqCPD', 'tGratingTemporalFreqCPS'};
     otherwise
         error('experiment type not yet defined')
 end

end


function get_ROI (source, eventdata)

% grab the udat
udat = get(gcf, 'UserData');

GetROI(udat);
end