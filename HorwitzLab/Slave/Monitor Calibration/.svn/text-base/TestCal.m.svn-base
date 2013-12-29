function TestCal(RGBs, cals, nrepeats)

% function TestCal(RGBs, cals[, nrepeats])
%
% This function is going to test the ability of the calibration files
% to predict the spectra obtained by various combinations of R,G,B
%
% If RGBs is a scalar, we pick the correspoding number of random RGB
% values.  If it is an Nx3 matrix, we use the passed in values as RGB
% triplets.
%
% cals should be a cell array of calibration structures.  Once the
% measurement is complete, we use these calibration structures to make and
% assess the quality of the spectrum predictions.
%
% The optional argument, nrepeats, gives the number of times to repeat
% measurement on a single RGB triplet (the mean if RGBs is a scalar, or
% the first row if RGBs is a matrix).  This give an estimate of the
% variance that cannot possibly be overcome by a better calibration
% structure.
%
% GDLH 3/2/07, ZALB 1/5/12 - PR-705 version
% ZALB 05/15/13 - Standardized code to support all devices

if nargin < 2
    error('Please pass in the # of test colors and the calibration data!');
elseif nargin < 3
    nrepeats = 0;
end

display_device = input('Which high color depth device are you using? [0=None, 1=Bits++Colour++, 2=ViewPixx, 3=Native10Bit] ');
if isempty(display_device) || ~ismember(display_device, 0:3)
    return
elseif display_device == 2
    scanning_mode = input('Do you want scanning backlight mode? [0=Off, 1=On] ');
    if isempty(scanning_mode) || ~ismember(scanning_mode, 0:1), scanning_mode = 1; end
else
    scanning_mode = -1;
end

meter_device = input('Which device are you using to calibrate your display? [1=PR650, 6=PR705] ');
if isempty(meter_device) || ~ismember(meter_device, [1 6]), meter_device = 0; end
switch meter_device
    case 1
        S_from = [380 4 101];
    case 6
        S_from = [380 2 201];
end

boxSize = 400;
Xoffset = 0;
Yoffset = 0;

% Setting up rgbs if they have not been explictly passed in.
if isscalar(RGBs)
    RGBs = rand(RGBs, 3);
end

% Optionally repeating a particular RGB triplet
RGBs = [repmat(RGBs(1,:), nrepeats, 1); RGBs];

disp(RGBs);

disp('Turn on Photometer and hit <return>');
disp('Then focus on white square and hit <return> again');
input('');

oldverbosity = Screen('Preference', 'Verbosity', 2);
oldsynclevel = Screen('Preference', 'SkipSyncTests', 2);

PsychImaging('PrepareConfiguration');
switch display_device
    case 0
    case 1, PsychImaging('AddTask', 'General', 'EnableBits++Color++Output', 1);
    case 2, PsychImaging('AddTask', 'General', 'EnableDataPixxC48Output', 1);
        % untested
    case 3, PsychImaging('AddTask', 'General', 'EnableNative10BitFramebuffer');
end
[window,screenRect] = PsychImaging('OpenWindow', 0, cals{end}.bgColor);

if display_device == 2
    SetOrToggleScanningBacklight(scanning_mode);
elseif display_device == 0
    Screen('ColorRange', window, 1); % notify Screen that RGB intensities max out at 1
    LoadIdentityClut(window);
end

boxRect = [0 0 boxSize boxSize];
boxRect = CenterRect(boxRect, screenRect);
boxRect = boxRect + [Xoffset -Yoffset Xoffset -Yoffset];
Screen('FillRect', window, [1 1 1], boxRect);
Screen('Flip', window);

input('');
% Open the serial port and set up the device
CMCheckInit(meter_device);

pause(10);  % Time to leave room
displayOrder = randperm(size(RGBs, 1));

spdMat = zeros(size(RGBs, 1), S_from(3));
for i = 1:size(RGBs, 1)
    Screen('FillRect', window, RGBs(displayOrder(i),:), boxRect);
    Screen('Flip', window);
    spd = MeasSpd(S_from, meter_device); % splining to S_to happens here
    if isempty(spd), break; end
    spdMat(displayOrder(i),:) = spd(:);
end

Screen('Preference', 'Verbosity', oldverbosity);
Screen('Preference', 'SkipSyncTests', oldsynclevel);
sca();
CMClose(meter_device);

if ~isempty(spd)
    try
        save_file_path = sprintf('%s-TestCalOut.mat', datestr(now, ...
            'yyyymmddTHHMMSS'));
        save(save_file_path, 'spdMat', 'RGBs', 'cals', 'scanning_mode', 'S_from');
    catch E
        if strcmp(E.identifier, 'MATLAB:save:permissionDenied')
            fprintf('cd to a directory with write permissions and execute a ''return''\n');
            keyboard
            save(save_file_path, 'spdMat', 'RGBs', 'cals', 'scanning_mode', 'S_from');
        else
            rethrow(E);
        end
    end
    
    %%% Analysis
    
    % do this in case someone executes this analysis code manually
    % using old data (before June 2013) where RGBs are 8 bit
    if find(RGBs > 1, 1)
        RGBs = RGBs / 255;
    end
    
    % spline all spectra to S_to, the most common S in the calibration file
    Cschemes = cellfun(@(x) x.S_device, cals, 'unif', 0);
    [nil,nil,scheme_idxs] = unique(vertcat(Cschemes{:}), 'rows');
    S_to = Cschemes{mode(scheme_idxs)};
    
    maxabsY = 0;
    prevexistingfigh = max(get(0, 'Children'));
    if isempty(prevexistingfigh)
        prevexistingfigh = 0;
    end
    for calIdx = 1:length(cals)
        cal = cals{calIdx};
        % map the RGBs to the nearest gamma table index
        [nil,GT_idxs] = histc(RGBs(:), ...
            [-inf (cal.gammaInput(1:end-1)' + cal.gammaInput(2:end)') / 2 inf]);
        GT_idxs = reshape(GT_idxs, size(RGBs));
        % this gives us linear indices so we can index the GT directly
        GT_idxs = bsxfun(@plus, GT_idxs, ...
            0:size(cal.gammaTable,1):numel(cal.gammaTable)-1);
        
        rgbbkgnd = FindModelWeights(cal.P_ambient, cal.P_device);
        % "actual" is a bit of a misnomer since we're finding the
        % coefficients on each of the guns by regression.
        actual = FindModelWeights(SplineSpd(S_from, spdMat', S_to), ...
            SplineSpd(cal.S_device, cal.P_device, S_to));
        
        figure(calIdx + prevexistingfigh);
        set(gcf, 'DefaultAxesColorOrder', eye(3));
        title(sprintf('%s: Calibration #%d', cals{1}.describe.monitor, calIdx));
        ylabel('Actual-pred');
        xlabel('Requested');
        plot(RGBs, actual' - bsxfun(@plus, cal.gammaTable(GT_idxs), rgbbkgnd'), '.', 'markersize', 10)
        maxabsY = max([maxabsY max(abs(get(gca, 'Ylim')))]);
    end
    for i = 1:length(cals)
        figure(i + prevexistingfigh);
        set(gca, 'Ylim', [-maxabsY maxabsY]);
    end
else
    warning('There was an empty measurement. No data saved!');
end
