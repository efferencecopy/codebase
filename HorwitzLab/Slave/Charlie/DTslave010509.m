function DTslave

% This slave program should accept commands from rex and present a flash of
% color in a particular location in space for a particular amount of
% time.
%
% CAH 11/07
% CAH 06/08 Major overhaul. Took out training content, converted to
%           nested functions and added functionality for adaptive proceedures.
% CAH 10/08 Changed trial type allocation (to a 4d matrix), ballanced
%           allocations to T1 and T2, enabled asymetric target locations. 

%***************  SETUP GLOBALS AND CONSTANTS *************%
%UDP COMMUNICATION INFO
udpCom.sock = [];
udpCom.port = 6665;
udpCom.rexip = '10.208.108.60';
udpCom.plexip = '10.208.108.61';

%MONITOR CALIBRATION INFO
ptb.w = [];
ptb.monWidth = [];
ptb.monHeight = [];
ptb.monCent = [];       %(x,y) of the screen center
ptb.pixperdeg = [];
ptb.oldGamma = [];
ptb.refreshRate = [];
ptb.maxdac = (2^16)-1;  %65535 dac values in colour mode
ptb.monSpd = [];
ptb.M = [];
ptb.invGamma = [];
ptb.bkgndrgb = [];
ptb.bkgndRGB = [];
ptb.gammaTable = [];
ptb.exitNow = 0;          %a flag that gets set to one when rex forces an exit back to master slave

%STIMULUS INFO
stim.blockCounter = [];         % gets decrimented on each block transition
stim.counters = [];             % counters for trial types
stim.exptDone = 0;
sitm.exptMethod = 0;            % 1=MoCS, 2=Staircasing
stim.flashFrameCounter = [];    % number of frames that flash has been on screen.
stim.flashTimeProf = [];        % temporal profile of the stimulus
stim.flashNumFrames = [];       % number of frames the flash is on the monitor
stim.flashRect = [];            % the texture drawing rectangle
stim.flashDeltaPhase = [];
stim.gaborTheta = 0;           % a float: the gabor's orientation in radians
stim.gaborLambda = 0;          % an int; wavelength in pixels
stim.gaborPhase = 0;           % a float; in radians
stim.gaborSigma = 0;           % a float; standard dev of gabor (comes in 10ths of degs, but converted to pix)
stim.gaborGamma = 0;           % a float; 1 give circles, else elipses
stim.halfTrials = [];           % the number of trials per block per target i.e 0.5(ntrials per block)
stim.fpRect = [];               % drawing rectangle for the FP
stim.flashPos = [];
stim.frameRectIpsi = [];        % drawing Rect for wire frame surrounding the flash
stim.frameRectContra = [];
stim.frameRGB = [];             % gets set up during monitor setup
stim.RFPos = [0 0];            % this gets set by updateGaborParams
stim.rgbLevels = [];            % the rgb triplets for each stimulus type
stim.spPeriodLevels = [];       % a vector of possible spatial periods (specifed by rex as spt. freqs and converted below)
stim.T1Rect = [];               
stim.T2Rect = [];

%TRIAL INFO
trial.drawFP = 0;
trial.drawFrame = 0;
trial.drawFlash = 0;
trial.drawTarg = 0;


%**********************************************************%

%Start by opening a connectionless udp socket
udpCom.sock = pnetStart(udpCom.port);
if (udpCom.sock == -1)
    error('UDP not properly initialized')
end
pnet(udpCom.sock, 'setreadtimeout', 0);
pnet(udpCom.sock, 'setwritetimeout', 0);



%main experimental loop
while(1)
    %check for udp messages
    dealWithUdpMessages();

    %draw a stimulus
    if (ptb.w)
        if(trial.drawFP)
            Screen('FillRect', ptb.w, [0 0 0], stim.fpRect);
        end

        if(trial.drawFrame)
            Screen('FrameRect', ptb.w, stim.frameRGB, stim.frameRectIpsi, 1);
            Screen('FrameRect', ptb.w, stim.frameRGB, stim.frameRectContra, 1);
        end

        if(trial.drawTarg)
            Screen('FillRect', ptb.w, stim.frameRGB, stim.T1Rect);
            Screen('FillRect', ptb.w, stim.frameRGB, stim.T2Rect);
        end

        if(trial.drawFlash)
            drawFlash();

            if(stim.flashFrameCounter == stim.flashNumFrames);
                pnet(udpCom.sock, 'write', 'MACstimLastFrame>> >>');
                pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
                trial.drawFlash = 0;
            end

            if(stim.flashFrameCounter == 1);
                pnet(udpCom.sock, 'write', 'MACstimFirstFrame>> >>');
                pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
            end

        end

        Screen('Flip', ptb.w);
    end
    
    %if rex wants to exit back to master slave:
    if ptb.exitNow
        return
    end
    
    %press ESC to shutdown
    [keyisdown,secs,keycode] = KbCheck();
    if (keycode(41)) %if you press escape or rex want's to bail
        if (ptb.w > 0);
            pnet(udpCom.sock, 'close');
            Screen('LoadNormalizedGammaTable', ptb.w, ptb.oldGamma);
            Screen('CloseAll');
            ShowCursor;
        end
        return
    end
end %main while loop



%
%       NESTED SUBFUNCTIONS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%*******************%
% STIUMLUS TOGGLES  %
%*******************%

%fixation point
    function fpToggle(onOff)
        trial.drawFP = onOff;
    end

%wireframe
    function wireFrameToggle(onOff)
        trial.drawFrame = onOff;
    end

%stimulus flash
    function flashToggle(onOff)
        trial.drawFlash = onOff;
    end

%targets
    function targToggle(onOff)
        trial.drawTarg = onOff;
    end




%**************************%
% MONITOR AND TRIAL SETUP  %
%**************************%

    function setupTrial(RFpos, flashLength, fpPos, fpSize, frameOffset, targSize, T1Dist, T2Dist, gaborTheta, gaborGamma, gaborSigma, gaborSpeed, nSigmas, lastTrialGood, lastTrialCorrect)
        %accepts parameters from rex then sets up the stimuli accordingly.
        %nSigmas, fpSize and frameOffset should be in tenths of degrees.
        %fpPos should be with respect to the center of the
        %screen in tenths of degrees. flashLength should be in units of
        %milli-seconds. target distances should be b/w 0&1 signifying the distance b/w
        %the fp and the flash represented in thousandths.

        %pick a color direction, contrast, and spatial frequency to present. This will
        %depend on weather your using an adaptive proceedure or method of
        %constant stimuli.
        switch stim.exptMethod
            case 1
                [stim.rgbTrial, stim.gaborLambda, stim.tColor, stim.tSfs, stim.tContrast, stim.tStimLoc] = MoCS(lastTrialGood);
            case 2
                [stim.rgbTrial, stim.gaborLambda, stim.tColor, stim.tSfs, stim.tContrast, stim.tStimLoc] = Staircase(lastTrialGood, lastTrialCorrect);
        end
        
        %convert the rgb for each trial into DAC values so that they can
        %get dropped into the data file
        if any(isnan(stim.rgbTrial))
            disp('End of Experiment');
            return
        end
        RGB = round(ptb.maxdac .* stim.rgbTrial) + 1;
        stim.flashRGB(1:3) = round(ptb.maxdac .* [ptb.invGamma(RGB(1), 1), ptb.invGamma(RGB(2), 2), ptb.invGamma(RGB(3), 3)]);
        
        %reset the frame counter on every trial
        stim.flashFrameCounter = 0;
        
        %set up the spatial location of the FP, wire frame, stimulus flash, and
        %saccade targets
        stim.fpRect = pos2rect(fpPos, fpSize);

        %deal with the flash location, which is chosen randomly by the
        %trialAllocator subfunction
        switch stim.tStimLoc
            case 1
                stim.flashPos = RFpos;
            case 2
                stim.flashPos = RFpos .* -1;
        end
        sizeInTenthsDeg = gaborSigma .* nSigmas .* 2; %multiply by two so that +/-nSigmas mass is shown
        stim.flashRect = pos2rect(stim.flashPos, sizeInTenthsDeg);
        stim.frameRectIpsi = pos2rect(stim.flashPos, sizeInTenthsDeg+frameOffset);
        stim.frameRectContra = pos2rect((stim.flashPos.*-1), sizeInTenthsDeg+frameOffset);

        %determine the location of the targets. Remeber that two right
        %triangles that have the same angles but different side lengths have
        %the same ratio of side lengths.
        T1PosInDeg = round(RFpos .* T1Dist);
        stim.T1Rect = pos2rect(T1PosInDeg, targSize);
        T2PosInDeg = round(-RFpos .* T2Dist);
        stim.T2Rect = pos2rect(T2PosInDeg, targSize);

        %deal with the time weighting function
        stim.flashNumFrames = ceil(ptb.refreshRate .* (flashLength./1000));
        rampLength = ceil(stim.flashNumFrames ./ 4);
        ramp = linspace(0, 1, rampLength);  %ramp is 1/6th of the total duration on either side
        plateau = ones(1,stim.flashNumFrames - (rampLength .* 2));
        stim.flashTimeProf = [ramp, plateau, fliplr(ramp)];

        %setup the gabor, get the spatial period (lambda) below
        stim.gaborTheta = gaborTheta;
        stim.gaborSigma = round(ptb.pixperdeg .* gaborSigma ./ 10);
        stim.gaborGamma = gaborGamma;
        stim.flashDeltaPhase = gaborSpeed .* 2 .* pi .* (1 ./ ptb.refreshRate); %the amount to advance each frame
    end

%METHOD OF CONSTANT STIMULI CONTRAST GENERATOR
    function [rgbTrial, lambda, color, sfs, contrast, targ] = MoCS(lastTrialGood)

        %decriment the appropriate counter if the last trial was a good one
        if lastTrialGood
            ind = find((stim.counters(:,1) == stim.tColor) & (stim.counters(:,2) == stim.tSfs) & (stim.counters(:,3) == stim.tContrast) & (stim.counters(:,4) == stim.tStimLoc));
            stim.counters(ind, end) = stim.counters(ind, end) - 1;
        end

        [color, sfs, contrast, targ] = trialAllocator();
        if ~any([color, sfs, contrast, targ])
            [rgbTrial, lambda] = deal(NaN); % crappy way to deal with the end of an expt
            return
        end


        %updates these things before the function returns
        rgbTrial = stim.rgbLevels{color}(contrast,:); %the rgb triplet for the trial type
        lambda = stim.spPeriodLevels(sfs);
    end

% STAIRCASING CONTRAST GENERATOR
    function [rgbTrial, lambda, color, sfs, contrast, targ] = Staircase(lastTrialGood, lastTrialCorrect)
        if (lastTrialGood)
            %decrement the appropriate counter
            ind = sub2ind(size(stim.counters), stim.tColor, stim.tSfs, stim.tContrast, stim.tStimLoc);
            stim.counters(ind) = stim.counters(ind) - 1;
            
            %update the stats for the most recent trial type
            if lastTrialCorrect
                %update the stats for that trial type
                stim.stairStats{stim.tColor, stim.tSfs}(1) = 1;
                stim.stairStats{stim.tColor, stim.tSfs}(1, 2:4) = ...
                    stim.stairStats{stim.tColor, stim.tSfs}(1, 2:4) - stim.stairDecLMS(:,stim.tColor)';
            else
                %explicitly set the stats for each trial outcome
                stim.stairStats{stim.tColor, stim.tSfs}(1) = 0;
                stim.stairStats{stim.tColor, stim.tSfs}(1, 2:4) = ...
                    stim.stairStats{stim.tColor, stim.tSfs}(1, 2:4) + stim.stairIncLMS(:,stim.tColor)';
            end
        end
        
        %now set up the next trial
        [color, sfs, contrast, targ] = trialAllocator();
        if ~any([color, sfs, contrast, targ])
            [rgbTrial, lambda, indicies] = deal(NaN); % crappy way to deal with the end of an expt
            return
        end
        trialLMS = stim.stairStats{color, sfs}(1, 2:4);

        %now check for out of bounds errors. don't let the contrast (vector
        %norm) go below bkgnd, and don't let the unit vectors change sign
        %(with respect to the one specified by rex)
        signFlips = any((trialLMS .* stim.stairMaxLMS(:, color)') < 0);
        if any(signFlips);
            trialLMS(1, 1:3) = 0; %if any of the CC fall below bkgnd, set them back to bkgnd
        end

        if norm(trialLMS) > norm(stim.stairMaxLMS(:,color));
            trialLMS(1, 1:3) = stim.stairMaxLMS(:,color)';
        end

        %reset the LMS triplet in the stim.stairStats vector.
        stim.stairStats{color, sfs}(1, 2:4) = trialLMS;

        %updates these things before the function returns
        %the rgb triplet for the trial type
        rgbTrial = [inv(ptb.M) * ((1+(trialLMS'./100)) .* ptb.bkgndLMS)]';
        indicies = [color, sfs, contrast, targ];
        lambda = stim.spPeriodLevels(sfs);
    end

% trial type allocator
    function [color, sfs, contrast, targ] = trialAllocator()
            %first, figure out which trial types are fair game. Don't advance
        %to the next block of trials until all trial types have been
        %exhausted for the current block
        typesAvailable = find(stim.counters(:,end));
        

        %if no types are available than advance to the next block.
        %Randomize the first trial type. If the last block just finished
        %than send the "all done" signal to REX.
        if isempty(typesAvailable)
            if (stim.blockCounter == 1) %all done?
                pnet(udpCom.sock, 'write', 'MACexptdone>> >>');
                pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
                stim.exptDone = 1;
                [color, sfs, contrast, targ] = deal(NaN); %return something so that the program doesn't bonk
                return
            else % decrement the block counter, reinitialize the tType counters
                stim.blockCounter = stim.blockCounter - 1;
                stim.counters(:,end) = stim.halfTrials;
            end

            %now call the trial allocator recursively to get the new trial type
            [color, sfs, contrast, targ] = trialAllocator();
            return %is this necessary for proper recursive calling????
        end

        %at long last, find a trial type and contrast level. Store these
        %indicies away so that on the next trial you can decriment the
        %appropriate counter!
        randomDraw = unidrnd(length(typesAvailable));
        tTypeChoiceInd = typesAvailable(randomDraw);
        color = stim.counters(tTypeChoiceInd, 1);
        sfs = stim.counters(tTypeChoiceInd, 2);
        contrast = stim.counters(tTypeChoiceInd, 3);
        targ = stim.counters(tTypeChoiceInd, 4);
    end


%MONITOR SETUP
    function setupMonitor(mondist, screenwidth)

        %load in the calibration data and fill up the ptb structure
        S = load('T_cones_smj.mat');
        load('/Monitor Calibration/Monitor data/Dell 4/Dell4BitsCal')
        calData = cals{end};
        fundamentals = S.T_cones_smj;
        ptb.bkgndRGB = round(255 .* calData.bgColor); %bkgnd voltages discritized b/w 0&255
        ptb.bkgndrgb = [calData.gammaTable(ptb.bkgndRGB(1)+1, 1), calData.gammaTable(ptb.bkgndRGB(2)+1, 2), ...
            calData.gammaTable(ptb.bkgndRGB(3)+1, 3)]; %add one to create and index b/w 1 & 256
        stim.frameRGB = ptb.bkgndRGB - 30;


        ptb.gammaTable = calData.gammaTable;
        ptb.monSpd = SplineSpd(calData.S_device, calData.P_device, S.S_cones_smj);
        ptb.M = fundamentals * ptb.monSpd;
        ptb.invGamma = InvertGamma(calData, 1);

        %open a ptb Screen window and determine various conversion factors
        if isempty(Screen('windows'));
            ptb.w = Screen('OpenWindow', 0, ptb.bkgndRGB);
        else
            ptb.w = Screen('windows');
            Screen('FillRect', ptb.w, ptb.bkgndRGB);
        end
        
            
        warning('off', 'MATLAB:conversionToLogical');
        HideCursor;
        
        [ptb.monWidth, ptb.monHeight] = Screen('WindowSize', ptb.w);
        theta = atand((screenwidth./2) ./ mondist);
        ptb.pixperdeg = (ptb.monWidth ./ 2) ./ theta;
        ptb.oldGamma = Screen('LoadNormalizedGammaTable', ptb.w, repmat(linspace(0,1,256)', 1, 3));
        ptb.refreshRate = 1 ./ Screen('GetFlipInterval', ptb.w); %the measured refresh rate not the reported rate (in Hz)
        ptb.monCent = [ptb.monWidth, ptb.monHeight] ./ 2; %(x,y) at the center of the screen

        %determine what the bkgndlms excitations are
        ptb.bkgndLMS = ptb.M * ptb.bkgndrgb';


        %set these up for the first trial. On subsequent trials they'll be set
        %up by AllOff()
        stim.gaborPhase = 0;
        stim.flashMaxInc = 0; %for debugging differences b/w expected and presented cone contrasts

    end


    function drawFlash()

        %increment the counter
        stim.flashFrameCounter = stim.flashFrameCounter+1;
        Screen('Close'); %clear the textures from the last refresh

        %start by making the gabor, then multiplying by the increment from
        %background and by the temporal weighting fxn.
        gabor = makeGabor(stim.gaborTheta, stim.gaborLambda, stim.gaborPhase, stim.gaborSigma, stim.gaborGamma);
        rgb_increment = stim.rgbTrial - ptb.bkgndrgb;
        flashImg = repmat(gabor, [1, 1, 3]);
        flashImg(:,:,1) = (flashImg(:,:,1) .* stim.flashTimeProf(stim.flashFrameCounter) .* rgb_increment(1)) + ptb.bkgndrgb(1);
        flashImg(:,:,2) = (flashImg(:,:,2) .* stim.flashTimeProf(stim.flashFrameCounter) .* rgb_increment(2)) + ptb.bkgndrgb(2);
        flashImg(:,:,3) = (flashImg(:,:,3) .* stim.flashTimeProf(stim.flashFrameCounter) .* rgb_increment(3)) + ptb.bkgndrgb(3);

        % this is a hacknied way dealing with LMS triplets that are out of the
        % monitor gammut, but will work for now.
        try
            %now convert to DAC values
            imgSize = size(gabor);
            flashrgb = round(ptb.maxdac .* flashImg) + 1; %intensities b/w 1 & maxdac +1
            flashRGB = ones(size(flashImg)); %preallocate so that indexing works


            flashRGB(:,:,1) = reshape(round(ptb.maxdac .* ptb.invGamma(flashrgb(:,:,1), 1)), imgSize(1), imgSize(2));
            flashRGB(:,:,2) = reshape(round(ptb.maxdac .* ptb.invGamma(flashrgb(:,:,2), 2)), imgSize(1), imgSize(2));
            flashRGB(:,:,3) = reshape(round(ptb.maxdac .* ptb.invGamma(flashrgb(:,:,3), 3)), imgSize(1), imgSize(2));
        catch
            fprintf('Out of bounds error. Are any of the Low Contrast Scalars Zero??\n');
            keyboard
        end

        flashTexture = Screen('MakeTexture', ptb.w, TranslateToColourMode(flashRGB));
        Screen('DrawTexture', ptb.w, flashTexture, [], stim.flashRect, [], 0);

        %update the position of the gabor for the next frame refresh
        stim.gaborPhase = stim.gaborPhase + stim.flashDeltaPhase;
        
        %find the index to the pixel with the largest contrast. at the end
        %of the trial, the RGB for this pixel will get dropped in the data
        %file.
        norms = sqrt(sum([flashImg .* flashImg], 3)); %in normalized gun space
        frameMaxInc = max(max(norms));
        if frameMaxInc > stim.flashMaxInc;
            [r, c] = find(norms == frameMaxInc, 1);
            stim.flashRGB(4) = flashRGB(r,c,1);  % DAC values b/w 0&maxdac
            stim.flashRGB(5) = flashRGB(r,c,2);
            stim.flashRGB(6) = flashRGB(r,c,3);
            stim.flashMaxInc = frameMaxInc;
        end
    end


%**************************%
% ACCESSORY FUNCTIONS      %
%**************************%
%BLANK THE SCREEN
    function allOff()

        trial.drawFP = 0;
        trial.drawFrame = 0;
        trial.drawFlash = 0;
        trial.drawTarg = 0;

        %reset these gloabal parameters for the next trial
        stim.gaborPhase = 0;
        Screen('Close') %I think that this closes/deletes all offscreen textures. redundant with above.
        stim.flashMaxInc = 0; %for debugging differences b/w expected and presented cone contrasts
    end

%MAKE A DRAWING RECTANGLE
    function rect = pos2rect(centPos, size)
        %draws a square of side length "size" (in pix) arround a point on the
        %monitor specified by centPos (in (x,y) cordinates) size is in
        %units of tenths of degrees

        %centPos and size come from rex in units of tenths of degrees so convert to
        %pixels first
        xyPosX = round(centPos(1) ./ 10 .* ptb.pixperdeg);
        xyPosY = round(centPos(2) ./ 10 .* ptb.pixperdeg);
        size = round(size ./ 10 .* ptb.pixperdeg);

        %make sure that the width of the rect is even (b/c of colour mode)
        if rem(size, 2);
            size = size + 1;
        end

        % now compute the rectangle.
        halfWidth = size./2;
        rect = [(ptb.monCent(1) + xyPosX-floor(halfWidth)), (ptb.monCent(2) - xyPosY-floor(halfWidth)), ...
            (ptb.monCent(1) + xyPosX+ceil(halfWidth)), (ptb.monCent(2) - xyPosY+ceil(halfWidth))];

        %because we're using colourmode the drawing rectangle must start on an
        %even (high order) byte pixel. If by chance the rect starts on and odd
        %pixel nudge the x cordinates of the UL & LR corners over one.
        if rem(rect(1), 2)
            rect = rect + [1, 0, 1, 0];
        end

    end

% DETERMNINE THE COLOR CONTRASTS
    function setupColorContrasts(LMSfromREX, LMSlowScalar, sptFreqfromREX, nSptFreqs, nContrasts, nTrialsPerBlock, nBlocks, expMethod)

        stim.exptDone = 0; %so that clicking reset states make this run o.k.
        stim.counters = []; % so that hitting reset states clears this out
        stim.blockCounter = nBlocks;
        stim.counters = []; % so that hitting reset states clears this out
        stim.halfTrials = ceil(nTrialsPerBlock./2); %so that trials get allocated to each target evenly
        
        

        %determine which experimental method to use. MoCS == 1, staircasing
        % == 2.
        stim.exptMethod = expMethod;

        %LMSfromREX should be a 9 element vector. Transform this so that
        %there are three rows and each triplet goes down a column
        LMSfromREX = reshape(LMSfromREX, 3,3);
        noColor = find(sum(abs(LMSfromREX)) == 0);
        if noColor
            LMSfromREX(:,noColor) = []; %delete LMS triplets that are all zeros
            LMSlowScalar(noColor) = []; %don't forget to delete the scalars
        end
        numColors = size(LMSfromREX,2);



        %now set up the possible spatial frequencies. Log space the
        %frequencies between the upper and lower bound.
        spFreqLevels = logspace(log10(sptFreqfromREX(1)), log10(sptFreqfromREX(2)), nSptFreqs);
        degPerCyc = 1./spFreqLevels;
        stim.spPeriodLevels = round(ptb.pixperdeg .* degPerCyc); %should I always round up to prevent cases where period would = 0?


        %set up the color triplets and convert to gun intensity space
        switch stim.exptMethod
            case 1 %MoCS
                for a = 1:numColors;
                    scaleVector = logspace(log10(LMSlowScalar(a)), log10(1), nContrasts);
                    coneLevels = LMSfromREX(:,a) * scaleVector;
                    
                    % Notice that I'm explicitly setting the 1st contrast to zero!!!
                    coneLevels(:,2:nContrasts+1) = coneLevels;
                    coneLevels(:,1) = 0;
                    
                    % convert from CC to cone intensities.
                    coneLevels = ((coneLevels ./ 100) + 1) .* repmat(ptb.bkgndLMS, 1, size(coneLevels, 2));

                    %lastly, convert to gun intensity space.
                    stim.rgbLevels{a} = [inv(ptb.M) * coneLevels]';
                end
            case 2 %staircasing
                %the staircasing will be conducted in CC space (i.e. on
                %each trial a certain amount of CC will be added to each of
                %the cones (which means I'll need to find the rgb triplets
                %later).
                
                nContrasts = 1; %just in case it was set to something else on rex.

                %start by finding the inc and dec CC values
                decval = 0.05; %2 thousandths of a %CC
                percentAtThresh = 0.816;
                incval = (decval .* percentAtThresh) ./ (1 - percentAtThresh);

                stim.stairStats = {}; % so that hitting reset states clears this out
                for a = 1:numColors;
                    %set up a stats matrix that keeps track or the
                    %performance on the last trial as well as the most
                    %recent LMS CC values
                    [stim.stairStats{a,1:nSptFreqs}] = deal([0, LMSfromREX(:,a)']);

                    %now determine the inc and dec sizes for each color
                    %direction such that the ratio L:M:S stays constant
                    %throughout the experiment
                    unitLMS = LMSfromREX(:,a) ./ norm(LMSfromREX(:,a));
                    stim.stairIncLMS(:,a) = unitLMS .* incval;
                    stim.stairDecLMS(:,a) = unitLMS .* decval;

                    %find the highest contrast that's still in the monitor
                    %gamut.
                    guess = [0; 0; 0]; %start at the white point (in CC space) and work out
                    while(1)
                        %deal with the peak of the gabor
                        peaklms = (1 + (guess./100)) .* ptb.bkgndLMS;
                        peakrgb = inv(ptb.M) * peaklms;
                        peakrgb = round(ptb.maxdac .* peakrgb) + 1; %indicies into invGamma.
                        
                        %deal with the trough
                        troughlms = (1 - (guess./100)) .* ptb.bkgndLMS;
                        troughrgb = inv(ptb.M) * troughlms;
                        troughrgb = round(ptb.maxdac .* troughrgb) + 1; %indicies into invGamma.
                        
                        %figure out what to do
                        peakInGamut = all(peakrgb < ptb.maxdac+1) && all(peakrgb > 0);
                        troughInGamut = all(troughrgb < ptb.maxdac+1) && all(troughrgb > 0);
                        
                        if (peakInGamut && troughInGamut)
                            guess = guess + stim.stairIncLMS(:,a);
                        else
                            guess = guess - stim.stairIncLMS(:,a);
                            break
                        end
                    end
                    stim.stairMaxLMS(:,a) = guess;
                end
        end
                
        %set up the trial type matrix differently for staircasing
        %and MoCS
        counters = fullfact([numColors, nSptFreqs, nContrasts, 2]); %2, one for each target
        switch stim.exptMethod
            case 1
                counters(:,3) = counters(:,3) + 1; % reserve '1' for the zero contrast condition
                counters((end+1):(end+2),:) = [1,1,1,1;...
                                           1,1,1,2];  %take the zero contrast at each of the spatial locations
                stim.counters = [counters, ones(size(counters,1), 1).*stim.halfTrials]; %make a column of counters.
            case 2
                stim.counters = [counters, ones(size(counters,1), 1).*stim.halfTrials]; %make a column of counters.
        end
    end



% CREATE A GABOR
    function gabor = makeGabor(theta, lambda, phi, sigma, gamma)
        halfSize = (stim.flashRect(3) - stim.flashRect(1)) ./ 2;  %base the size of the gabor off of the flash's texture window so that there are no mistakes
        row = -halfSize:2:halfSize-1; %subtract one so that you don't overlap the texture window;
        col = -halfSize:halfSize-1;
        [X, Y] = meshgrid(row, col);
        xprime = X .* cos(-theta) + Y .* sin(-theta);
        yprime = -X .* sin(-theta) + Y .* cos(-theta);
        gabor = exp(-(xprime.^2 + gamma.^2 .* yprime.^2) ./ (2.* sigma.^2)) .* cos(2 .* pi .* yprime ./ lambda + phi);
    end


% GET/SET GABOR GLOBALS
    function loadGaborFits(fitMethod);
        %fitMethod (1) = Grating, (2) = WN online

        global gl;


        %pull in the params and convert to units that REX uses. The logic
        %here is that REX feeds these parameters to the Mac on every trial.
        %So I'll send these to REX at the beginning of a run, and REX will
        %send them back when the time comes.
        if fitMethod == 1;
            if ~isfield(gl, 'grating');
                return
            end
            stim.gaborTheta = gl.grating.preforient;                     % radians
            stim.gaborLambda = (1./gl.grating.prefsf);                   % deg/cyc
            stim.gaborPhase = 0;
            stim.gaborSigma = (gl.grating.prefdiam ./ 6) .* 10;          % 10ths of degs. dividing by six is arbitratry scaling of sigma
            stim.gaborGamma = 1;                                     % no other choice unless we modle the WN data
            stim.RFPos = [gl.grating.x, gl.grating.y] .* 10;         % comes in dva. local standard is 10ths dva
        elseif (fitMethod == 2)
            if ~isfield(gl, 'gabor');
                return
            end
            stim.gaborTheta = gl.gabor.theta;                           % radians
            stim.gaborLambda = gl.gabor.lambda;                         % deg/cyc
            stim.gaborPhase = 0;
            stim.gaborSigma = gl.gabor.sigma .* 10;                     % DTdefault is 10ths of degrees but WN fits with degrees
            stim.gaborGamma = gl.gabor.gamma;                           % 1 gives circles, else elipses

            %now determine the X,Y for the stimulus center. I'm going to send
            %these values to REX, so they need to be in tenths of degrees.
            x = [gl.stim.x + gl.gabor.xoffset]  .* 10; %these globals are in dva, so multiply by 10 to obtain
            y = [gl.stim.y + gl.gabor.yoffset]  .* 10;  %the local standard of 10ths of dva
            stim.RFPos = [x,y]; %needs to be in tenths of degrees for use locally
        end


    end


%HANDLE UDP COMMUNICATIONS
    function dealWithUdpMessages()

        msgSize = pnet(udpCom.sock, 'readpacket', 260, 'no block');
        if ~msgSize
            return
        end

        try
            %open the message
            message = pnet(udpCom.sock, 'read', msgSize, 'char');
            
            %see if we should return to MasterSlave
            if strcmpi('return', message(1:6));
                ptb.exitNow = 1;
            end
            
            %try to execute the message
            %disp(message); %for debugging
            eval(message);
            
        catch
            fprintf('Unknown message: %s\n', message);
            error = lasterror;
            disp(error.message);
            disp(error.identifier);
            for a = 1:size(error.stack,1)
                disp(error.stack(a,1));
            end
        end
    end

end %DTslave


%%
%******* NOTES AND QUESTIONS **************
% 1) previous editions had a "blank texture" do I need one here?

%3) optional arg to framerate?

%11) check method of creating targRects (plus they seem to be offset)

%13) out of bounds errrors currently catch CC lower than zero, but there
%may be occasions where negative CC are desireable??? If so, than I need
%different error checking.

% out or bounds errors occur when the LMS triplet is a NaN. The occurs b/c
% rex desires O CC conditions (log(0) is a NaN).


% does stim.RFPos differ from the RFpos dumped by REX on each trial? I'm
% nervous that recasting to an int is a bad idea
