function iafGetRetinotopy

%
%  SET UP THE IMPORTANT PARAMETERS THAT DEFINE THE ANALYSIS
%

fnames = {'data-i539-160303-1130.mat', '160303_i539_lateral_az[0 30]el[10 -10]_run1_1_MMStack.tif';...
    'data-i539-160303-1207.mat', '160303_i539_lateral_az[0 30]el[10 -10]_run2_1_MMStack.tif'};


relevantTrialTypes = {'tGratingElevationDeg', 'tGratingAzimuthDeg'}; % options below
% tBaseGratingDirectionDeg
% tGratingAzimuthDeg
% tGratingElevationDeg
% tGratingDirectionDeg
% tGratingContrast
% tGratingDiameterDeg
% tGratingSpatialFreqCPD
% tGratingTemporalFreqCPS
% tGratingStartingPhaseDeg
% tGratingSpeedDPS


%
%  COMPUTE THE RETINOTOPY
%

% extract the runs, separate the data by trial type
dat = iafExtractRuns(fnames, relevantTrialTypes);
N_types = size(dat.uniqueTrlTypes,1);

% compute the average movie for each trial type
AvgTrialON = {};
AvgTrialOFF = {};
for i_t_types = 1 : N_types    
    AvgTrialON{i_t_types} = iafGetAvgMovie(dat.OnResponse{i_t_types});
    AvgTrialOFF{i_t_types} = iafGetAvgMovie(dat.OffResponse{i_t_types});
end


% Average across time, for each (average) ttype movie. Then combine the on
% and off responses
final_img = {};
for i_type = 1:N_types
    tmp = AvgTrialON{i_type};
    for i_frame = 1:size(tmp,3);
        tmp(:,:,i_frame) = imgaussfilt(tmp(:,:,i_frame), 1);
    end
    meanAcrossTime_on = mean(tmp, 3);
    
    tmp = AvgTrialOFF{i_type};
    for i_frame = 1:size(tmp,3);
        tmp(:,:,i_frame) = imgaussfilt(tmp(:,:,i_frame), 1);
    end
    meanAcrossTime_off = mean(tmp, 3);
    
    final_img{i_type} = meanAcrossTime_on - meanAcrossTime_off;
end


% make a figure showing the final images.
h.fig = figure;
h.fig.Position = [847   199   454   475];
h.fig.Units = 'Normalized';
h.ax = axes('position', [0.125 0.200 0.750    0.75]);
h.dat = final_img;
h.ttypes = dat.uniqueTrlTypes;
h.text = dat.trlTypeHeader;
h.slide = uicontrol('style', 'slider', 'units', 'normalized', 'position', [0.25, 0.05, 0.50, 0.05],...
                     'Min', 1, 'Max', numel(final_img), 'sliderstep', ones(1,2).*(1/(numel(final_img)-1)),...
                     'Value', 1, 'Interruptible', 'on', 'Callback', @updateImage, 'Enable', 'on',...
                     'ButtonDownFcn', @updateImage);
h.fig.UserData = h;
updateImage();

end

function updateImage(one, two)
    h = get(gcf, 'UserData');
    imgNum = h.slide.Value;
    imgNum = max([1, imgNum]); % slider step requires range, but don't let the img num < 0
    imgNum = round(imgNum);
    axes(h.ax)
    imagesc(h.dat{imgNum}); colormap('gray')
    h.ax.XTick = [];
    h.ax.YTick = [];
    
    % specify the stimulus position
    elidx = strcmpi(h.text, 'tGratingElevationDeg');
    azidx = strcmpi(h.text, 'tGratingAzimuthDeg');
    el = h.ttypes(imgNum, elidx);
    az = h.ttypes(imgNum, azidx);
    h.ax.XLabel.String = sprintf('Az: %.0f, El: %.0f', az, el);
end








