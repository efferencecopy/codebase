%#ok<*DEFNU,*ASGLU,*NASGU,*USENS>
function DTslave_scot()
% This slave program should accept commands from rex and present a flash of
% color in a particular location in space for a particular amount of
% time.
%
% CAH 11/07
% CAH 06/08 Major overhaul. Took out training content, converted to
%           nested functions and added functionality for adaptive proceedures.
% CAH 10/08 Changed trial type allocation (to a 4d matrix), ballanced
%           allocations to T1 and T2, enabled asymetric target locations.
% CAH 02/09 Added functionality for QUEST adaptive proceedures.
% CAH 04/09 Added "catch trials" where the online gui can force the
%           presentation of the optimized drifting grating.
%
% ZLB 02/12 A scotopic detection task with colored squares on black
%           background. Later added ViewPixx functionality.

%***************  SETUP GLOBALS AND CONSTANTS *************%
% UDP COMMUNICATION INFO
udpCom.sock = [];
udpCom.port = 6665;
udpCom.rexip = '192.168.1.120';

KbName('UnifyKeyNames');
key_esc = KbName('escape');

% MONITOR CALIBRATION INFO
ptb.w = [];
ptb.mon_width = [];
ptb.mon_height = [];
ptb.mon_center = []; % (x,y) of the screen center
ptb.mon_spd = [];
ptb.pixperdeg = [];
ptb.oldGamma = [];
ptb.refreshRate = [];
ptb.maxdac = (2^16) - 1; % 65536 dac values in colour mode
ptb.invGamma = [];
ptb.bkgndrgb = [];
ptb.gammaTable = [];
ptb.exitNow = 0;
ptb.vpixx = ConnectedToViewPixx();
ptb.ccmode = 1;

% STIMULUS INFO
stim.blockcount = []; % gets decremented on each block transition
stim.counters = []; % counters for trial types
stim.done = 0;
stim.nframes = 0;
stim.framecount = 0; % number of frames that flash has been on screen.
stim.rect = []; % the texture drawing rectangle
stim.rgb = [0 0 0]; % the actual RGB of the stimulus
stim.halfTrials = []; % the number of trials per block per target i.e 0.5(ntrials per block)
stim.fprect = []; % drawing rectangle for the FP
stim.position = [];
stim.t1rect = [];
stim.t2rect = [];
stim.texture = [];
stim.image = [];
stim.targrgb = [1 1 1]; % white saccade targets
stim.rgb_fp = [255 255 255];

% QUEST PARAMETERS
q.muPrior = [0 0 0]; % best guesses for threshold (1 x nColors)
q.sdPrior = []; % stand dev common to all priors of thresholds
q.modBeta = []; % slope of psy fun common to all stim conditions
q.Q = {}; % quest functions (one for each stim condition)
q.domain = {}; % range of possible thresholds (one for each color)

% TRIAL INFO
trial.drawFP = 0;
trial.drawStim = 0;
trial.drawTargs = 0;

%**********************************************************%
% Start by opening a connectionless udp socket
udpCom.sock = pnetStart(udpCom.port);
if udpCom.sock == -1
    error('UDP not properly initialized');
end

pnet(udpCom.sock, 'setreadtimeout', 0);
pnet(udpCom.sock, 'setwritetimeout', 0);

messageIsAvailable = 0;
% main experimental loop
while true
    % check for udp messages
    messageIsAvailable = pnet(udpCom.sock, 'readpacket', ...
        1000, 'noblock');
    
    if messageIsAvailable
        dealWithUdpMessages(messageIsAvailable);
        messageIsAvailable = 0;
    end
    
    % draw a stimulus
    if ptb.w
        if trial.drawFP
            % for DrawLines: [xstart xend; ystart yend] and remember y
            % increases going down the monitor (opposite of cartesian in y).
            spokelength = 5; % in degrees I think
            screenDelta = repmat([ptb.fpPos(1); -ptb.fpPos(2)] / spokelength / 10, 1, 2);
            Screen('DrawLines', ptb.w, spokelength * ptb.pixperdeg * ([-1  1; -1 1] + screenDelta), 10, stim.rgb_fp, ptb.mon_center);
            Screen('DrawLines', ptb.w, spokelength * ptb.pixperdeg * ([ 1 -1; -1 1] + screenDelta), 10, stim.rgb_fp, ptb.mon_center);
            Screen('FillRect', ptb.w, [0 0 0], stim.fprect)
            annurect = pos2rect(ptb.fpPos, 2 * ptb.fpSize);
            Screen('FrameOval', ptb.w, stim.rgb_fp, annurect, round((annurect(3)-annurect(1))/20));
        end
        
        if trial.drawTargs
            if stim.correct_only && stim.location_idx == 1
                Screen('FillRect', ptb.w, stim.targrgb, stim.t1rect);
            elseif stim.correct_only && stim.location_idx == 2
                Screen('FillRect', ptb.w, stim.targrgb, stim.t2rect);
            else
                Screen('FillRect', ptb.w, stim.targrgb, stim.t1rect);
                Screen('FillRect', ptb.w, stim.targrgb, stim.t2rect);
            end
        end
        
        if trial.drawStim
            % increment the counter. It starts at zero at the beginning of
            % each trial. Then show the appropriate frame of the movie
            stim.framecount = stim.framecount + 1;
            
            Screen('DrawTexture', ptb.w, stim.texture, [], stim.rect, [], 0);
            
            % send the '1st frame' msg before the flip. if you wait till
            % after the flip the time of the stimulus presentation will be
            % misrepresented
            if stim.framecount == 1
                pnet(udpCom.sock, 'write', 'MACstimFirstFrame>> >>');
                pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
            end
        end
        
        if trial.drawStim && stim.framecount == stim.nframes
            Screen('Close', stim.texture);
            stim.texture = [];
        end
        
        Screen('Flip', ptb.w);
        
        % send the 'last frame' message after the last flip. otherwise the
        % time of stimulus presentation will be misestimated by Rex.
        if trial.drawStim
            if stim.framecount == stim.nframes
                pnet(udpCom.sock, 'write', 'MACstimLastFrame>> >>');
                pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
                trial.drawStim = 0;
            end
        end
    end
    
    % if rex wants to exit back to master slave
    if ptb.exitNow, return; end
    
    % press ESC to shutdown
    [keyisdown,secs,keycode] = KbCheck();
    if keyisdown && keycode(key_esc) % if you press escape
        if ptb.w > 0
            pnet(udpCom.sock, 'close');
			sca();
			ptb.exitNow = 1;
        end
        return
    end
end % main while loop

% begin nested functions
    function fpToggle(onOff)
        trial.drawFP = onOff;
    end

    function stimToggle(onOff)
        trial.drawStim = onOff;
    end

    function targToggle(onOff)
        trial.drawTargs = onOff;
    end

    function setupColors(rex_intensities, num_trialsperblock, num_blocks, correct_only, stim_placement)
        % initialize these parameters so that clicking reset states makes
        % the program run smoothly.
        stim.done = 0;
        stim.counters = [];
        stim.blockcount = num_blocks;
        stim.halfTrials = ceil(num_trialsperblock / 2); % so that trials get allocated to each target evenly
        
        stim.correct_only = correct_only;
        
        % rex_intensities should be a 9 element vector. Transform this so that
        % there are three rows and each triplet goes down a column
        rex_intensities = reshape(rex_intensities, 3, 3);
        noColor = find(sum(abs(rex_intensities)) == 0);
        if noColor
            rex_intensities(:,noColor) = []; % delete LMS triplets that are all zeros
            q.muPrior(noColor) = [];    % delete these too!
        end
        
        norms = sqrt(diag(rex_intensities' * rex_intensities));
        units = rex_intensities' ./ repmat(norms, 1, 3); % triplets go across rows!
        
        % set up the first quest function for each color.
        q.colorDirs = units;
        for c = 1:size(q.colorDirs, 1);
            q.domain{c} = linspace(.001, 1, 1000);
            prior = 1 / (q.sdPrior * sqrt(2*pi)) * exp(-(q.domain{c} - q.muPrior(c)).^2 / (2 * q.sdPrior^2));
            q.Q{c} = log(prior);
        end
        
        % set up the trial type matrix
        num_colors = size(rex_intensities, 2);
        counters = fullfact([num_colors 2]); % 2, one for each target
        switch stim_placement
            case 1
                counters(:,end) = 2;
            case 2
                counters(:,end) = 1;
        end
        stim.counters = [counters  ones(size(counters, 1), 1) * stim.halfTrials]; % make a column of counters.
    end

    function setupTrial(stimPos, stimOnTime, stimSize, fpPos, fpSize, targSize, T1Dist, T2Dist, lastTrialGood, lastTrialCorrect)
        [stim.rgbTrial, stim.colordir_idx, stim.location_idx] = ...
            Quest(lastTrialGood, lastTrialCorrect);
        
        if any(isnan(stim.rgbTrial)), return; end
        
        RGB = round(ptb.maxdac * stim.rgbTrial) + 1;
        stim.rgb = [ptb.invGamma(RGB(1),1) ...
                      ptb.invGamma(RGB(2),2) ...
                      ptb.invGamma(RGB(3),3)];
        
        if stim.location_idx == 1
            stim.position = stimPos;
        else
            stim.position = -stimPos;
        end
        
        stim.fprect = pos2rect(fpPos, fpSize);
        ptb.fpPos = fpPos;
        ptb.fpSize = fpSize;
        stim.rect = pos2rect(stim.position, stimSize);
        T1PosInDeg = round(stimPos * T1Dist);
        T2PosInDeg = round(stimPos * T2Dist);
        stim.t1rect = pos2rect(T1PosInDeg, targSize);
        stim.t2rect = pos2rect(-T2PosInDeg, targSize);
        stim.nframes = ceil(ptb.refreshRate * stimOnTime / 1000);
        stim.framecount = 0;
        
        % make square texture
        sizeinpix = (stim.rect(3) - stim.rect(1));

        stim.image = zeros(sizeinpix, sizeinpix/2, 3);
        
        stim.image = bsxfun(@plus, stim.image, reshape(stim.rgb, [1 1 3]));
        stim.texture = Screen('MakeTexture', ptb.w, stim.image, [], [], 2);
    end

    function [rgbTrial, colordir, location] = Quest(lastTrialGood, lastTrialCorrect)
        % deal with the previous trial's counter and quest functions
        if lastTrialGood
            % update the counter
            [t, idx] = ismember([stim.colordir_idx stim.location_idx], ...
                stim.counters(:,1:end-1), 'rows');
            
            stim.counters(idx, end) = stim.counters(idx, end) - 1;
            
            % update the quest function
            normIntensity = norm((stim.rgbTrial(:)));
            liklihood = .98 - 0.48 * exp(-(normIntensity ./ q.domain{stim.colordir_idx}) .^ q.modBeta);
            if ~lastTrialCorrect
                liklihood = 1 - liklihood;
            end
            q.Q{stim.colordir_idx} = q.Q{stim.colordir_idx} + log(liklihood);
        end
        
        % obtain stimulus parameters for the upcoming trial
        [colordir, location] = trialAllocator();
        if ~any([colordir location])
            rgbTrial = NaN; % crappy way to deal with the end of an expt
            return
        end
        
        %mean_estimate = sum(q.Q{colordir} .* q.domain{colordir}) / sum(q.Q{colordir});
        
        % variance = sum((q.domain{colordir} - mean_estimate) .^ 2 .* q.Q{colordir}) / sum(q.Q{colordir});   
        %rgbTrial = q.colorDirs(colordir,:) * mean_estimate;
        [val, idx] = max(q.Q{colordir});
        rgbTrial = q.colorDirs(colordir,:) * q.domain{colordir}(idx);
        
    end

    function [colordir, location] = trialAllocator()
        % first, figure out which trial types are fair game. Don't advance
        % to the next block of trials until all trial types have been
        % exhausted for the current block
        typesAvailable = find(stim.counters(:,end));
        
        % if no types are available than advance to the next block.
        % Randomize the first trial type. If the last block just finished
        % than send the "all done" signal to REX.
        if isempty(typesAvailable)
            if stim.blockcount == 1 % all done?
                pnet(udpCom.sock, 'write', 'MACexptdone>> >>');
                pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
                stim.done = 1;
                [colordir, location] = deal(NaN); %return something so that the program doesn't bonk
                return
            else % decrement the block counter, reinitialize the tType counters
                stim.blockcount = stim.blockcount - 1;
                stim.counters(:,end) = stim.halfTrials;
            end
            
            % now call the trial allocator recursively to get the new trial type
            [colordir, location] = trialAllocator();
            return % is this necessary for proper recursive calling????
        end
        
        % at long last, find a trial type and contrast level. Store these
        % indicies away so that on the next trial you can decrement the
        % appropriate counter!
        randomDraw = unidrnd(length(typesAvailable));
        tTypeChoiceInd = typesAvailable(randomDraw);
        colordir = stim.counters(tTypeChoiceInd,1);
        location = stim.counters(tTypeChoiceInd,2);
    end

%
%  SET QUEST PARAMETERS
%
% Simply sets the local quest parameters
% based on those specified by REX. This function gets called before
% setupColorContrasts.
%  priorSD      => scalar in units of gun intensity??
%  beta         => scalar; slope of the models psyFun
%  alpha        => vector; guesses of alpha for each color
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function setQuestParams(priorSD, beta, alpha)
        q.muPrior = alpha'; % only one estimate per gun
        q.sdPrior = priorSD;
        q.modBeta = beta;
    end
    
    function setupMonitor(mondist, screenwidth, calFilePath)        
        %load in the calibration data and fill up the ptb structure
        
        load(calFilePath);
        calData = cals{end}; %#ok<USENS>
        
        ptb.bkgndrgb = [0 0 0];
        ptb.realbkgndRGB = round(255 .* calData.bgColor);
        ptb.realbkgndrgb = [calData.gammaTable(ptb.realbkgndRGB(1)+1, 1), calData.gammaTable(ptb.realbkgndRGB(2)+1, 2), ...
            calData.gammaTable(ptb.realbkgndRGB(3)+1, 3)];
        
        ptb.gammaTable = calData.gammaTable;
        ptb.mon_spd = calData.P_device;      
        ptb.invGamma = InvertGammaTable(calData.gammaInput, calData.gammaTable, ptb.maxdac+1);      

        % start up the imaging pipeline
        if ~isempty(Screen('Windows'))
            ptb.w = max(Screen('Windows'));
            Screen('FillRect', ptb.w, ptb.bkgndrgb);
        else
            PsychImaging('PrepareConfiguration');
            if ptb.vpixx
                PsychImaging('AddTask', 'General', 'EnableDataPixxC48Output', ptb.ccmode); % in mode '1': every 2nd column of pixels is ignored
            else
                PsychImaging('AddTask', 'General', 'EnableBits++Color++Output', ptb.ccmode);
            end
            ptb.w = PsychImaging('OpenWindow', 0, ptb.bkgndrgb);
        end
        
        HideCursor();
        [ptb.mon_width, ptb.mon_height] = Screen('WindowSize', ptb.w);
        theta = atand(screenwidth / 2 / mondist);
        ptb.pixperdeg = ptb.mon_width / 2 / theta;
        ptb.refreshRate = Screen('NominalFrameRate', ptb.w, 1);
        ptb.mon_center = [ptb.mon_width ptb.mon_height] / 2; % (x,y) at the center of the screen
        
        % inform REX that the monitor setup is complete
        pnet(udpCom.sock, 'write', 'MonitorSetupComplete>> >>');
        pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
    end

    function allOff()
        trial.drawFP = 0;
        trial.drawStim = 0;
        trial.drawTargs = 0;
        if ~isempty(stim.texture)
            Screen('Close', stim.texture);
            stim.texture = [];            
        end
    end

%
%  SEND QUEST STATS TO REX
%
%   this fxn only gets called by REX during quest experiments and sends
%   accross the modal value of the appropriate posterior distribution. Due
%   to the current architecture of Rex's state set (and how I update the
%   quest stats) this number actually reflects the modal value of the most
%   recent trial NOT the current trial. The current trials results won't be
%   known until 'setupTrial' gets called for the next trial.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function updateQuestStats()
        %         [m, idx] = max(q.Q{stim.colordir_idx});
        %         pThreshMode = q.domain{stim.colordir_idx}(idx);
        %pThreshMean = sum(q.Q{stim.colordir_idx} .* q.domain{stim.colordir_idx}) / sum(q.Q{stim.colordir_idx});
        
        % use send to rex to transmit the mean value of the posterior
        %sendToRex(udpCom, pThreshMean, 'double', 'updateQuestStats();');
        
        
        [m, idx] = max(q.Q{stim.colordir_idx});
        pThreshMode = q.domain{stim.colordir_idx}(idx);
        
        % use send to rex to transmit the modal value of the posterior
        sendToRex(udpCom, pThreshMode, 'double', 'updateQuestStats();');
   
        
    end

%
%   SEND QUEST COLOR RANGES TO REX
%
%  When using quest it's necessary to represent the domain of each quest
%  function (the same domain is used for the liklihood functions). Instead
%  of sending the entire string of numbers I'll just send the first, last,
%  and stepsize. I'll send over the data for each color direction tested as
%  a vector of numbers. Each successive triplet is for a new color. Since
%  Rex is expecting a 9 vector i'll pad with zeros for the cases where I
%  use fewer than 3 colors.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function sendQuestColorRanges()
        msg = zeros(1, 9);
        nColors = length(q.domain);
        for i = 1:nColors;
            msg(i * 3 - 2) = min(q.domain{i});
            msg(i * 3 - 1) = max(q.domain{i});
            msg(i * 3) = q.domain{i}(2) - q.domain{i}(1); % pretty kludgy
        end
        
        % now send the message using sendToRex
        sendToRex(udpCom, msg, 'double', 'sendQuestColorRanges();');
    end

    function rect = pos2rect(centPos, rect_size)
        % centPos and size come from rex in units of tenths of degrees so convert to
        % pixels first
        xyPosX = round(centPos(1) / 10 * ptb.pixperdeg);
        xyPosY = round(centPos(2) / 10 * ptb.pixperdeg);
        % use 2 * pixperdeg for y component on ViewPixx
        rect_size = round(rect_size / 10 * ptb.pixperdeg);
        
        % make sure that the width of the rect is even (b/c of colour mode)
        if rem(rect_size, 2), rect_size = rect_size + 1; end
        
        rect_template = rect_size * [-1 -1 1 1] / 2;
        % now translate the rectangle into place
        rect = CenterRectOnPoint(rect_template, ptb.mon_center(1) + xyPosX, ...
            ptb.mon_center(2) - xyPosY);
%         halfWidth = rect_size / 2;
%         rect = [ptb.mon_center(1) + xyPosX - floor(halfWidth/(1+ptb.vpixx)) ...
%             ptb.mon_center(2) - xyPosY - floor(halfWidth) ...
%             ptb.mon_center(1) + xyPosX + ceil(halfWidth/(1+ptb.vpixx)) ...
%             ptb.mon_center(2) - xyPosY + ceil(halfWidth)];
        
        % because we're using colourmode the drawing rectangle must start on an
        % even (high order) byte pixel. If by chance the rect starts on and odd
        % pixel nudge the x cordinates of the UL & LR corners over one.
        if rem(rect(1), 2)
            rect = rect + [1 0 1 0];
        end
    end

    function dealWithUdpMessages(msgSize)
        % open the message
        message = pnet(udpCom.sock, 'read', msgSize, 'char');
        
        % see if we should return to MasterSlave
        if strncmpi(message, 'return', 6)
            stk = dbstack();
            if ~strcmp(stk(end).name, mfilename)
                ptb.exitNow = 1;
            end
        end
        
        try
            eval(message);
        catch ME
            fprintf('Trouble with message: "%s"\n', message);
            disp(getReport(ME));
        end
    end

end % DTslave_scot
