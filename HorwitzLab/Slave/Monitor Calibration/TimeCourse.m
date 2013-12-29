% TimeCourse.m
%
% Looking at how the intensity of the three monitor phosphors vary across time.
% Waiting for things to stabilize, then changing the background and waiting
% for things to stabilize again.  The question is - when I change the
% background does the monitor take a while to "adjust" or does it just need
% to warm up when I first turn it on and from then on I can change the
% background as much as I want without having to worry that I should let it
% stabilize again.
%
% GDLH 05/23/07
%
% This should work with the Bits++ in colour mode.  We're just looking at
% the stability of the phosphor intensites (so the precise voltages to the
% guns aren't critical as long as they are fixed).  All that's going to
% happen is the high byte will be the same as the low byte, which is fine.
% The flag for COLOURMODE, below, is just passed on to AdaptColSearch.m,
% for which it is more important.
%
% ZALB 05/09/13
%
% Letting the user specify which display and PR device to use. I've also
% updated the display setup code to use `PsychImaging`.
%
% GDLH 06/07/07
% ZALB 07/13/12, 05/09/13

% Setting stuff up
boxSize = 200;
Xoffset = 0;
Yoffset = 0;
TOL = 0.005; % Proportion tolerance for knowing when we've stabilized
PERIODMINS = 5; % Sampling period in minutes

S_to = [380 5 81];

display_device = input('Which high color depth device are you using? [0=None, 1=Bits++, 2=ViewPixx, 3=Native10Bit] ');
if isempty(display_device) || ~ismember(display_device, 0:3)
    return
elseif display_device == 2
    scanning_mode = input('Do you want scanning backlight mode? [0=Off, 1=On] ');
    if isempty(scanning_mode) || ~ismember(scanning_mode, 0:1), scanning_mode = 1; end
else
    scanning_mode = -1;
end

meter_device = input('Which PR device is hooked up? [1=PR-650, 6=PR-705] ');
if isempty(meter_device) || ~ismember(meter_device, [1 6])
    return
end

BKGNDRGBS = [190 190 190; 100 100 100; 255 0 0; 0 255 0; 0 0 255];
BOXRGBS = [255 0 0; 0 255 0; 0 0 255; 255 255 255];

if display_device > 0
    BKGNDRGBS = BKGNDRGBS / 255;
    BOXRGBS = BOXRGBS / 255;
end

disp('Turn on Photometer and hit <return>');
disp('Then focus on white square and hit <return> again');
input('');

CMCheckInit(meter_device);
if meter_device == 6
    PR705config('Cycles', 3); % lets use the average of 5 measurements
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
    [window,screenRect] = PsychImaging('OpenWindow', 0, BKGNDRGBS(1,:));
    if display_device == 2
        scanning_mode = SetOrToggleScanningBacklight(scanning_mode);
    end
else
    [window,screenRect] = Screen('OpenWindow', 0, BKGNDRGBS(1,:));
    LoadIdentityClut(window);
end
HideCursor();

boxRect = boxSize * [-1 -1 1 1] / 2;
boxRect = CenterRectOnPoint(boxRect, screenRect(3)/2 + Xoffset, ...
    screenRect(4)/2 - Yoffset);

Screen('FillRect', window, BOXRGBS(end,:), boxRect);
Screen('Flip', window);

input('');
pause(2);

Rspd = []; Gspd = []; Bspd = []; Whitespd = []; timevect = []; rgbs = [];
bkgndidx = 1;
while bkgndidx <= size(BKGNDRGBS,1)
    Screen('FillRect', window, BKGNDRGBS(bkgndidx,:), screenRect);    
    Screen('FillRect', window, BOXRGBS(1,:), boxRect);
    Screen('Flip', window);
    
    [spd,nil] = MeasSpd(S_to, meter_device);
    if isempty(spd), break; end
    Rspd = [Rspd;spd'];
    
    Screen('FillRect', window, BOXRGBS(2,:), boxRect);
    Screen('Flip', window);
    
    [spd,nil] = MeasSpd(S_to, meter_device);
    if isempty(spd), break; end
    Gspd = [Gspd;spd'];
    
    Screen('FillRect', window, BOXRGBS(3,:), boxRect);
    Screen('Flip', window);
    
    [spd,nil] = MeasSpd(S_to, meter_device);
    if isempty(spd), break; end
    Bspd = [Bspd;spd'];
    
    Screen('FillRect', window, BOXRGBS(4,:), boxRect);
    Screen('Flip', window);
    
    [spd,nil] = MeasSpd(S_to, meter_device);
    if isempty(spd), break; end
    Whitespd = [Whitespd;spd'];
    
    tmp = clock;
    timevect = [timevect; tmp(4) tmp(5)];
    pause(PERIODMINS*60);
    
    rgbs = [rgbs; BKGNDRGBS(bkgndidx,:)];
    
    % A little online data analysis...
    if size(Rspd,1) > 1 % if we have at least two measurements
        changes = [(sum(Rspd(end,:))-sum(Rspd(end-1,:)))./sum(Rspd(end,:)),...
            (sum(Gspd(end,:))-sum(Gspd(end-1,:)))./sum(Gspd(end,:)),...
            (sum(Bspd(end,:))-sum(Bspd(end-1,:)))./sum(Bspd(end,:))];
        if all(abs(changes) < TOL)
            bkgndidx = bkgndidx+1;
        end
    end
end

Screen('Preference', 'Verbosity', oldverbosity);
Screen('Preference', 'SkipSyncTests', oldsynclevel);
sca();

if ~isempty(spd)
    save_str = sprintf('%s-TimeCourse.mat', datestr(now, 'yyyymmddTHHMMSS'));
    save(save_str, 'Rspd', 'Gspd', 'Bspd', 'Whitespd', 'timevect', 'rgbs', 'scanning_mode');
    
    %%% Analysis
    avgRspd = mean(Rspd);
    avgGspd = mean(Gspd);
    avgBspd = mean(Bspd);
    avgWspd = mean(Whitespd);
    
    figure; hold on;
    plot(avgRspd,'r-');
    plot(avgGspd,'g-');
    plot(avgBspd,'b-');
    plot(avgWspd,'k-');
    
    tmp = diff(timevect(:,1));
    daychanges = find(tmp == -23); % Day changes
    for i = 1:length(daychanges)
        timevect(daychanges(i)+1:end,1) = ...
            timevect(daychanges(i)+1:end,1) + 24;
    end
    
    figure; hold on;
    weights = [];
    for i = 1:size(Rspd,1)
        weights = [weights regress(avgRspd',Rspd(i,:)')];
    end
    plot(timevect(:,1)+timevect(:,2)*(1/60),weights,'r.-');
    weights = [];
    for i = 1:size(Gspd,1)
        weights = [weights regress(avgGspd',Gspd(i,:)')];
    end
    plot(timevect(:,1)+timevect(:,2)*(1/60),weights,'g.-');
    weights = [];
    for i = 1:size(Bspd,1)
        weights = [weights regress(avgBspd',Bspd(i,:)')];
    end
    plot(timevect(:,1)+timevect(:,2)*(1/60),weights,'b.-');
    
    ylabel('Weight to get to average');
    xlabel('Time (hours)');
    
    %%% Test of linearity for white
    pred = avgRspd+avgGspd+avgBspd;
    figure; axes; hold on;
    plot(pred,'b-');
    plot(avgWspd,'k-');
    [weights] = regress(avgWspd',[avgRspd' avgGspd' avgBspd']);
    gunnames = {'red','green','blue'};
    for i = 1:3
        if weights(i) < 1
            percentchange = 100*(1-weights(i));
            fprintf('%s gun is %0.2g%% too dim in white\n', gunnames{i}, percentchange);
        else
            percentchange = 100*(weights(i)-1);
            fprintf('%s gun is %0.2g%% too bright in white\n', gunnames{i}, percentchange);
        end
    end
    
    %%%% More analysis
    figure; axes; hold on;
    newtimevect = cumsum(diff(timevect(:,1)+timevect(:,2)*(1/60)));
    changes = any(diff(rgbs)~=0,2);
    
    plot(newtimevect,diff(sum(Rspd,2))./sum(Rspd(2:end,:),2),'r.-');
    plot(newtimevect,diff(sum(Gspd,2))./sum(Gspd(2:end,:),2),'g.-');
    plot(newtimevect,diff(sum(Bspd,2))./sum(Bspd(2:end,:),2),'b.-');
    plot(newtimevect(changes == 1),0,'k*');
    plot([min(newtimevect) max(newtimevect)],[TOL TOL],'k:');
    plot([min(newtimevect) max(newtimevect)],[-TOL -TOL],'k:');
    title('(x(2)-x(1))/x(2)');
    xlabel('Time (hours)');
    ylabel('prop. change from prev. measurement');
    
    background_rgb = AdaptColSearch([.33 .33], meter_device, display_device);
    save(save_str, 'background_rgb', '-append');
else
    warning('No data saved: there was a measurement error!');
end

CMClose(meter_device);
