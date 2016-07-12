function iafGetRetinotopy(experimentType, METHOD)

% store the inputs for later
udat.experimentType = experimentType;
udat.method = METHOD;

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
switch lower(METHOD)
    case 'intrinsic'
        dat = iafExtractRuns(udat.fid.names_mat, udat.fid.names_img, udat.fid.ttypes);
        artifactDefault = 4;
        smoothingDefault = 2;
    case 'calcium'
        dat = iafExtractRuns_offPeriod(udat.fid.names_mat, udat.fid.names_img, udat.fid.ttypes);
        artifactDefault = 'off';
        smoothingDefault = 'off';
end
N_types = size(dat.uniqueTrlTypes,1);

% save the ImagingRate in the output of this function so that the user can 
% analyze this data later
udat.frameRate = dat.frameRate;

% store only the average across trials (separated into ON trials OFF trials, and by stimulus
% type)
for i_type = 1 : N_types
    udat.AvgTrialON{i_type} = iafGetAvgMovie(dat.OnResponse{i_type});
    udat.AvgTrialOFF{i_type} = iafGetAvgMovie(dat.OffResponse{i_type});
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

% Name the figure according to filename
musID = regexpi((regexpi(udat.fid.names_mat{1}, 'data-(\w\d+)', 'match')), '(\d+)', 'match');
set(udat.h.fig, 'Name', sprintf('%s %s', musID{1}{1}, experimentType), 'NumberTitle', 'off');

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
                          'string', artifactDefault,...
                          'Callback', @img_preProcess);
uicontrol('style', 'text',...
          'units', 'normalized',...
          'position', [0.12, 0.1, 0.12, 0.05], ...
          'string', 'Artifact',...
          'FontSize', 10);
udat.h.smoothBox = uicontrol('style', 'edit',...
                             'units', 'normalized',...
                             'position', [0.75, 0.11, 0.1, 0.05],...
                             'string', smoothingDefault,...
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
    
    % pull out the correct 'final image'
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
    
    if strcmpi(udat.method, 'intrinsic')
        preProcessed_ON = cellfun(@(x) iafRemoveOutliers(x, artifactMultiplier), udat.AvgTrialON, 'uniformoutput', false);
        preProcessed_OFF = cellfun(@(x) iafRemoveOutliers(x, artifactMultiplier), udat.AvgTrialOFF, 'uniformoutput', false);
    else
        preProcessed_ON = udat.AvgTrialON;
        preProcessed_OFF =  udat.AvgTrialOFF;
    end
    
    
    % now do the pre-processing, start by filtering the average images
    Nttypes = numel(udat.AvgTrialON);
    if filterSD > 0
        udat.h.ax.Title.String = 'Filtering the raw data'; drawnow
        preProcessed_ON = cellfun(@(x) imgaussfilt(x, filterSD), preProcessed_ON, 'uniformoutput', false);
        preProcessed_OFF = cellfun(@(x) imgaussfilt(x, filterSD), preProcessed_OFF, 'uniformoutput', false);
        udat.h.ax.Title.String = ''; drawnow
    end
    
    % Compute an image to show in the gui
    for i_type = 1 : Nttypes

        meanAcrossTime_on = nanmean(preProcessed_ON{i_type}, 3);
        meanAcrossTime_off = nanmean(preProcessed_OFF{i_type}, 3);
        
        if strcmpi(udat.method, 'calcium')
            final_img = meanAcrossTime_on;
        else
            final_img = meanAcrossTime_on - meanAcrossTime_off;
        end
        
        
        % deal with outliers in the calcium imaging data
        if strcmpi(udat.method, 'calcium')
            l_oob = abs(final_img)>artifactMultiplier;
            replacementVal = nanmean(final_img(~l_oob));
            if artifactMultiplier > 0 && ~isnan(artifactMultiplier)
                final_img(l_oob) = replacementVal;
            end
            
            % figure out how to deal with NaNs
            l_nans = isnan(final_img);
            final_img(l_nans) = replacementVal;
        end
        
        
        % store the final image
        udat.final_img{i_type} = final_img;
        
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
     case 'size'
         relevantTrialTypes = {'tGratingDiameterDeg'};
     case 'annulus'
         relevantTrialTypes = {'tAnnulusGratingDiameterDeg', 'tGratingDiameterDeg'};
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
