function out_rgb = AdaptColSearch(target_xy, meter_device, display_device)

% rgb = AdaptColSearch(target_xy [, meter_device] [, display_device])
%
% Finds rgb values that provide a light with the desired CIE x,y
% chromaticity coordinates that is closest to r = g = b = 0.5 (all
% phosphors at half maximum).
%
% The point is to find an rgb triplet that we can use in CalibrateMonSpd,
% thereby calibrating against a background that is something known and not
% specific to the device.
%
% ZALB 05/09/13
%
% Letting the user specify which display and PR device to use. I've also
% updated the display setup code to use `PsychImaging`.
%
% GDLH 5/6/07
% ZALB 7/13/12, 05/09/13

boxSize = 200;
MAXITER = 10;
TOLVAL = .001;
s = load('T_xyz1964');
T_xyz1964 = s.T_xyz1964;
S_xyz1964 = s.S_xyz1964;
clear('s');

scanning_mode = -1;

wait_for_user_input = 0;
if nargin < 2 || isempty(meter_device)
    wait_for_user_input = 1;
    meter_device = input('Which PR device is hooked up? [1=PR-650, 6=PR-705] ');
    if isempty(meter_device) || ~ismember(meter_device, [1 6])
        return
    end
    CMCheckInit(meter_device);
    if meter_device == 6
        PR705config('Cycles', 5); % lets use the average of 5 measurements
    end
end

if nargin < 3 || isempty(display_device)
    wait_for_user_input = 1;
    display_device = input('Which high color depth device are you using? [0=None, 1=Bits++, 2=ViewPixx, 3=Native10Bit] ');
    if isempty(display_device) || ~ismember(display_device, 0:3)
        return
    elseif display_device == 2
        scanning_mode = input('Do you want scanning backlight mode? [0=Off, 1=On] ');
    end
end

primaries = zeros(S_xyz1964(end), 3);

if wait_for_user_input
    disp('Turn on Photometer and hit <return>');
    disp('Then focus on white square and hit <return> again');
    input('');
end

oldverbosity = Screen('Preference', 'Verbosity', 2);
oldsynclevel = Screen('Preference', 'SkipSyncTests', 2);

if display_device > 0
    PsychImaging('PrepareConfiguration');
    switch display_device
        case 1, PsychImaging('AddTask', 'General', 'EnableBits++Color++Output', 1);
        case 2, PsychImaging('AddTask', 'General', 'EnableDataPixxC48Output', 1);
        % untested
        case 3, PsychImaging('AddTask', 'General', 'EnableNative10BitFramebuffer');
    end
    [window,screenRect] = PsychImaging('OpenWindow', 0, .6);
    if scanning_mode > -1, SetOrToggleScanningBacklight(scanning_mode); end
    primaryrgbs = eye(3);
else
    [window,screenRect] = Screen('OpenWindow', 0, [150 150 150]);
    Screen('LoadNormalizedGammaTable', window, linspace(0,1,256)' * ones(1,3));
    primaryrgbs = 255*eye(3);
end
HideCursor();

screenWidth = screenRect(3);
screenHeight = screenRect(4);

boxRect = boxSize * [-1 -1 1 1] / 2;
boxRect = CenterRectOnPoint(boxRect, screenWidth/2, ...
    screenHeight/2);

if wait_for_user_input
    Screen('FillRect', window, (254 * (display_device == 0) + 1) * [1 1 1], boxRect);
    Screen('Flip', window);
    input('');
    pause(2);
end

for i = 1:3
    Screen('FillRect', window, primaryrgbs(i,:));
    Screen('Flip', window);
    [spd,nil] = MeasSpd(S_xyz1964, meter_device);
    if isempty(spd), break; end
    primaries(:,i) = spd;
end

if ~isempty(spd)
    XYZ = T_xyz1964 * primaries;
    maxlum = sum(XYZ(2,:));
    x = XYZ(1,:) ./ sum(XYZ);
    y = XYZ(2,:) ./ sum(XYZ);
    figure; axes; hold on;
    plot([x x(1)], [y y(1)]);
    plot(target_xy(1), target_xy(2), 'm*');
    
    % Making sure requested color is in the monitor gamut
    hullpts = convhull([x target_xy(1)], [y target_xy(2)]);
    if any(hullpts == 4)
        error('The requested x,y coordinates are not inside monitor gamut.')
    end
    
    % Conversion from x,y,luminance to XYZ
    targetlum = maxlum / 2;
    XplusZ = targetlum * (1 - target_xy(2)) / target_xy(2);
    factor1 = target_xy(1) * targetlum / (1 - target_xy(1));
    factor2 = 1 + ((1 - target_xy(2)) / target_xy(2));
    factor3 = 1 + (target_xy(1) / (1 - target_xy(1)));
    X = factor1 * factor2 / factor3;
    Z = XplusZ - X;
    
    % Sanity checks
    disp([target_xy(1) X / (X + targetlum + Z)]);
    disp([target_xy(2) targetlum / (X + targetlum + Z)]);
    
    % The vector in rgb space that corresponds to the
    % target x and y chromaticity coordinates (at half max lum.)
    predrgb = XYZ \ [X targetlum Z]';
    
    % Taking the projection onto predrgb vector from [.5 .5 .5]
    % to find the rgbs that are on the appropriate line and
    % closest to [.5 .5 .5]
    projrgb = predrgb .* ([.5 .5 .5] * predrgb) / (predrgb' * predrgb);
    
    % Iterative search to find the DAC values.
    tol = sum(projrgb .* TOLVAL); % Change TOLVAL for ViewPixx?
    if display_device > 0
        guess = projrgb;
        guesses = [1 1 1; 0 0 0];
    else
        guess = round(projrgb * 255);
        guesses = [255 255 255; 0 0 0];
    end
    err = 1;
    iter = 0;
    allweights = [1 1 1; 0 0 0];
    
    while all(err > tol)
        Screen('FillRect', window, guess);
        Screen('Flip', window);
        
        [spd,nil] = MeasSpd(S_xyz1964, meter_device);
        if isempty(spd), break; end
        
        weights = regress(spd, primaries);
        err = norm(weights - projrgb);
        allweights = [allweights; weights'];
        guesses = [guesses; guess'];
        
        try
            if display_device > 0
                % do I need to worry about ties?
                guess = [interp1(allweights(:,1),guesses(:,1),projrgb(1));
                    interp1(allweights(:,2),guesses(:,2),projrgb(2));
                    interp1(allweights(:,3),guesses(:,3),projrgb(3))];
            else
                noise = normrnd(0,1e-9,size(guesses,1),3); % avoiding problems with ties
                guess = [round(interp1(allweights(:,1),guesses(:,1)+noise(:,1),projrgb(1)));
                    round(interp1(allweights(:,2),guesses(:,2)+noise(:,2),projrgb(2)));
                    round(interp1(allweights(:,3),guesses(:,3)+noise(:,3),projrgb(3)))];
            end
        catch
            sca
            keyboard
        end
        update = guesses(end,:) - guess';
        
        guess = max(guess, [0 0 0]');
        if display_device > 0
            guess = min(guess, [1 1 1]');
        else
            guess = min(guess, [255 255 255]');
        end
        disp(guess');
        if iter == MAXITER || ...
                (display_device > 0 && all(update == 0)) || ...
                (display_device == 0 && all(round(update) == 0))
            err = 0;
        end
        iter = iter + 1;
    end
    
    if ~isempty(spd)
        newXYZ = T_xyz1964*spd;
        x = newXYZ(1) / sum(newXYZ);
        y = newXYZ(2) / sum(newXYZ);
        
        figure; axes;
        set(gcf, 'DefaultAxesColorOrder', [1 0 0; 0 1 0; 0 0 1]);
        plot(guesses);
        
        out_rgb = guess;
        fprintf('The final RGBs are: %.5g %.5g %.5g\n', guess');
        fprintf('target (x,y): (%.5g,%.5g), actual (x,y): (%.5g,%.5g)\n', target_xy(:)', [x y]);
        fprintf('target (r,g,b): (%.5g,%.5g,%.5g), actual (r,g,b): (%.5g,%.5g,%.5g)\n', projrgb(:)', weights(:)');
    end
end

if display_device == 0
    Screen('LoadNormalizedGammaTable', window, oldClut);
end

if wait_for_user_input
    CMClose(meter_device);
end

Screen('Preference', 'Verbosity', oldverbosity);
Screen('Preference', 'SkipSyncTests', oldsynclevel);
sca();
