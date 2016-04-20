function iafGetRetinotopy(experimentType)


% specify the actual trial types from the experiment type
udat.fid.ttypes = trialTypeLibrary(experimentType);

% has the user specified the file names? if not, call uigetfile and
% figure out which .mat files correspond to which .tiff files. do this in a
% sub-function
[udat.fid.names_mat, udat.fid.names_img, udat.fid.path] = ioGetFnames;


%
%  COMPUTE THE RETINOTOPY
%

% extract the runs, separate the data by trial type
cd(udat.fid.path)
dat = iafExtractRuns_offPeriod(udat.fid.names_mat, udat.fid.names_img, udat.fid.ttypes);
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

%
% Initialize the GUI window and display the final images
%
udat.ttypes = dat.uniqueTrlTypes;
udat.text = dat.trlTypeHeader;
udat.h.fig = figure;
udat.h.fig.Units = 'Normalized';
udat.h.fig.Position = [0.1526    0.2450    0.3211    0.6188];
udat.h.ax = axes('position', [0.15 0.35 0.70 0.60]);
if N_types > 1
    udat.h.slide = uicontrol('style', 'slider',...
                             'units', 'normalized',...
                             'position', [0.30, 0.2, 0.40, 0.05],...
                             'Min', 1, 'Max', N_types,...
                             'sliderstep', ones(1,2).*(1/(N_types-1)),...
                             'Value', 1,...
                             'Interruptible', 'on',...
                             'Callback', @updateImage,...
                             'Enable', 'on',...
                             'ButtonDownFcn', @updateImage);
end
udat.h.artBox = uicontrol('style', 'edit',...
                          'units', 'normalized',...
                          'position', [0.25, 0.11, 0.1, 0.05],...
                          'string', 'off',...
                          'Callback', @img_preProcess);
uicontrol('style', 'text',...
          'units', 'normalized',...
          'position', [0.12, 0.1, 0.12, 0.05], ...
          'string', 'Artifact',...
          'FontSize', 10);
udat.h.smoothBox = uicontrol('style', 'edit',...
                             'units', 'normalized',...
                             'position', [0.75, 0.11, 0.1, 0.05],...
                             'string', 1,...
                             'Callback', @img_preProcess);
uicontrol('style', 'text',...
          'units', 'normalized',...
          'position', [0.56, 0.1, 0.18, 0.05], ...
          'string', 'Smoothing',...
          'FontSize', 10)
udat.h.PushBtn = uicontrol('style', 'pushbutton',...
                           'units', 'Normalized',...
                           'String', 'Get R.O.I.',...
                           'Position', [[0.375   0.02   0.25   0.05]],...
                           'Callback', @get_ROI);
udat.h.fig.UserData = udat;
img_preProcess();

end



function updateImage(~,~)
    
    udat = get(gcf, 'UserData');
    
    N_types = numel(udat.final_img);
    if N_types > 1
        imgNum = udat.h.slide.Value;
        imgNum = max([1, imgNum]); % slider step requires range, but don't let the img num < 0
        imgNum = round(imgNum);
    else
        imgNum = 1;
    end
    
    axes(udat.h.ax)
    imagesc(udat.final_img{imgNum}); colormap('gray');
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

function img_preProcess(~,~)
    
    % grab the udat
    udat = get(gcf, 'UserData');
    
    % figure out the status of all the ui controls
    filterSD = round(str2double(udat.h.smoothBox.String));
    artifactMultiplier = str2double(udat.h.artBox.String);
    
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
    
    global GL_DATPATH
    
    % call uigetfile. with default directory. Save
    [names, path] = uigetfile({'*.mat;*.tif',...
                               'Related Files (*.mat;*.tif)'},...
                               'Select the Related MATLAB & TIFF files',...
                               'MultiSelect', 'on',...
                               '\\crash.dhe.duke.edu\data\home\kyra\Data');
    
    % check to make sure the user selected a valid set of files
    assert(iscell(names) || names~=0, 'ERROR: No files were selected.')
    assert(numel(names)>=2 && rem(numel(names),2)==0, 'ERROR: Must select at least 2 files, and each .mat file must have associated .tiff file');
    
    
    names_mat = names(strncmpi('data', names, 4));
    names_img = names(strncmpi('tiff', names, 4));
    assert(numel(names_mat) == numel(names_img), 'ERROR: each .mat file must have associated .tiff file')

    % make sure that the order of .mat files corresponds to the order of
    % .tiff files: make sure each file comes from the same run.
    Nnames = numel(names_mat);
    for i_names = 1:Nnames
        
        tag_mat = regexpi(names_mat{i_names}, 'run[\d].', 'match');
        tag_img = regexpi(names_img{i_names}, 'run[\d].', 'match');
        
        assert(strcmpi(tag_mat, tag_img),...
            'The MATLAB files do not match up with the TIFF files!')
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

function get_ROI(~,~)

    % grab the udat
    udat = get(gcf, 'UserData');

    % strip out the udat gui handels, so that the new gui has a clean slate
    udat.h = [];

    % call the gui that selects ROIs
    GetROI(udat);

end
