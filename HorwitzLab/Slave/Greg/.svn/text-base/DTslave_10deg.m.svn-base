function DTslave_10deg

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

global DT_movie;


%***************  SETUP GLOBALS AND CONSTANTS *************%
%UDP COMMUNICATION INFO
udpCom.sock = [];
udpCom.port = 6665;
udpCom.rexip = '192.168.1.120';
udpCom.plexip = '192.168.1.121';

%MONITOR CALIBRATION INFO
ptb.w = [];
ptb.monWidth = [];
ptb.monHeight = [];
ptb.monCent = [];               %(x,y) of the screen center
ptb.pixperdeg = [];
ptb.oldGamma = [];
ptb.refreshRate = [];
ptb.maxdac = (2^16)-1;          %65536 dac values in colour mode
ptb.monSpd = [];
ptb.M = [];
ptb.invGamma = [];
ptb.bkgndlms = [];
ptb.bkgndrgb = [];
ptb.bkgndRGB = [];
ptb.gammaTable = [];
ptb.exitNow = 0;                %a flag that gets set to one when rex forces an exit back to master slave

%STIMULUS INFO
stim.blockCounter = [];         % gets decrimented on each block transition
stim.counters = [];             % counters for trial types
stim.exptDone = 0;
sitm.exptMethod = 0;            % 1=MoCS, 2=Staircasing
stim.flashFrameCounter = 0;     % number of frames that flash has been on screen.
stim.flashTimeProf = [];        % temporal profile of the stimulus
stim.flashNumFrames = [];       % number of frames the flash is on the monitor
stim.flashRect = [];            % the texture drawing rectangle
stim.flashDeltaPhase = [];
stim.flashMovie = [];
stim.flashMaxInc = 0;
stim.gaborTheta = 0;            % a float: the gabor's orientation in radians
stim.gaborLambda = 0;           % an int; wavelength in pixels
stim.gaborPhase = 0;            % a float; in radians
stim.gaborSigma = 0;            % a float; standard dev of gabor (comes in 10ths of degs, but converted to pix)
stim.gaborGamma = 0;            % a float; 1 give circles, else elipses
stim.gratingColor = [];         % the color for the grating catch trials (color of highest firing rate in grating.d)
stim.halfTrials = [];           % the number of trials per block per target i.e 0.5(ntrials per block)
stim.fpRect = [];               % drawing rectangle for the FP
stim.flashPos = [];
stim.frameRectIpsi = [];        % drawing Rect for wire frame surrounding the flash
stim.frameRectContra = [];
stim.frameRGB = [];             % gets set on each trial (frames can either be present or absent)
stim.RFPos = [0 0];             % this gets set by loadGaborFits
stim.rgbLevels = [];            % the rgb triplets for each stimulus type
stim.spPeriodLevels = [];       % a vector of possible spatial periods (specifed by rex as spt. freqs and converted below)
stim.T1Rect = [];               
stim.T2Rect = [];
stim.targRGB = [];              % gets set during monior setup.


% QUEST PARAMETERS
q.muPrior = [0 0 0];            % best guesses for threshold (1 x nColors)
q.sdPrior = [];                 % stand dev common to all priors of thresholds
q.modBeta = [];                 % slope of psy fun common to all stim conditions
q.Q = {};                       % quest functions (one for each stim condition)
q.range = {};                   % range of possible thresholds (one for each color); given by monitor gamut


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
            Screen('FillRect', ptb.w, [0,0,0], stim.fpRect);
        end

        if(trial.drawFrame)
            Screen('FrameRect', ptb.w, stim.frameRGB, stim.frameRectIpsi, 1);
            Screen('FrameRect', ptb.w, stim.frameRGB, stim.frameRectContra, 1);
        end

        if(trial.drawTarg)
            Screen('FillRect', ptb.w, stim.targRGB, stim.T1Rect);
            Screen('FillRect', ptb.w, stim.targRGB, stim.T2Rect);
        end

        if(trial.drawFlash)
            %increment the counter. It starts at zero at the beginning of
            %each trial. Then show the appropriate frame of the movie
            stim.flashFrameCounter = stim.flashFrameCounter+1;
            Screen('DrawTexture', ptb.w, stim.flashMovie(stim.flashFrameCounter), [], stim.flashRect, [], [0]);
            Screen('Close', stim.flashMovie(stim.flashFrameCounter));
            
            % send the '1st frame' msg before the flip. if you wait till
            % after the flip the time of the stimulus presentation will be
            % misrepresented
            if (stim.flashFrameCounter == 1);
                pnet(udpCom.sock, 'write', 'MACstimFirstFrame>> >>');
                pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
            end
        end

        Screen('Flip', ptb.w);
        
        % send the 'last frame' message after the last flip. otherwise the
        % time of stimulus presentation will be misestimated by Rex.
        if trial.drawFlash
            if(stim.flashFrameCounter == stim.flashNumFrames);
                pnet(udpCom.sock, 'write', 'MACstimLastFrame>> >>');
                pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
                trial.drawFlash = 0;
            end
        end

    end
    
    %if rex wants to exit back to master slave:
    if ptb.exitNow
        return
    end
    
    %press ESC to shutdown
    [keyisdown,secs,keycode] = KbCheck();
    if (keycode(41)) %if you press escape
        if (ptb.w > 0);
            pnet(udpCom.sock, 'close');
            Screen('LoadNormalizedGammaTable', ptb.w, ptb.oldGamma);
            Screen('CloseAll');
            ShowCursor;
        end
        return
    end
end %main while loop


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    function wireFrameToggle(onOff, showFrames)
        trial.drawFrame = onOff;
        if showFrames
            stim.frameRGB = ptb.bkgndRGB - 30;
        else
            stim.frameRGB = ptb.bkgndRGB;
        end
    end

%stimulus flash
    function flashToggle(onOff)
        trial.drawFlash = onOff;
    end

%targets
    function targToggle(onOff)
        trial.drawTarg = onOff;
    end




%***************************************%
% MONITOR, EXPERIMENT, AND TRIAL SETUP  %
%***************************************%

%
%         SETUP THE COLOR DIRECTIONS AND TRIAL COUNTERS
%
% Rex specifies all the color directions and spatial frequencies. Take
% these parameters and convert to units used locally on the slave program.
% Lastly, set up a counters matrix that keeps track of the trial types
% presented. Input arguments have the following form:
%  LMSfromREX       => 9 vector of LMS CC triplets (up to 3 dirs)
%  LMSlowScalar     => 3 vector, each element scales 1 of the above colors
%  sptFreqfromREX   => 2 vector, low and high sfs to test
%  nSptFreqs        => scalar, how many sfs to test
%  nContrasts       => slave adds the "zero" contrast (n+1) cntrsts total
%  nTrialsPerBlock  => half this number are allocated to each targ location
%  nBlocks          => number of block per experiment
%  expMethod        => (mocs = 1) (stairs = 2)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function setupColorContrasts(LMSfromREX, LMSlowScalar, sptFreqfromREX, nSptFreqs, nContrasts, nTrialsPerBlock, nBlocks, expMethod)
        
        % initialize these parameters so that clicking reset states makes
        % the program run smoothly.
        stim.exptDone = 0;
        stim.counters = [];
        stim.blockCounter = nBlocks;
        stim.halfTrials = ceil(nTrialsPerBlock./2); %so that trials get allocated to each target evenly
        stim.exptMethod = expMethod;
        
        %set up the possible spatial frequencies. Log space the
        %frequencies between the upper and lower bound. Don't let the user
        %specify a spt freq of zero
        sptFreqfromREX(~sptFreqfromREX) = eps; %no zeros allowed.
        spFreqLevels = logspace(log10(sptFreqfromREX(1)), log10(sptFreqfromREX(2)), nSptFreqs);
        degPerCyc = 1./spFreqLevels;
        stim.spPeriodLevels = round(ptb.pixperdeg .* degPerCyc); %in pix

        %LMSfromREX should be a 9 element vector. Transform this so that
        %there are three rows and each triplet goes down a column
        LMSfromREX = reshape(LMSfromREX, 3,3);
        noColor = find(sum(abs(LMSfromREX)) == 0);
        if noColor
            LMSfromREX(:,noColor) = []; %delete LMS triplets that are all zeros
            LMSlowScalar(noColor) = []; %don't forget to delete the scalars
            q.muPrior(noColor) = [];    %delete these too!
        end
        numColors = size(LMSfromREX,2);
        
        
        %determine the extent of the monitor gamut for the specified color
        %directions
        norms = sqrt(diag(LMSfromREX' * LMSfromREX));
        units = LMSfromREX' ./ repmat(norms, 1, 3); %triplets go across rows!
        for c = 1:size(units, 1);
            guess = 0;
            incSize = 0.1; %1 tenth of a %CC
            while(1)
                tmpLMS = units(c,:) .* guess;

                %peak of the gabor:
                peakrgb = inv(ptb.M)*((1+tmpLMS./100)' .* ptb.bkgndlms);
                peakRGB = round(ptb.maxdac .* peakrgb) + 1;
                peakInGamut = all(peakRGB < (ptb.maxdac+1)) && all(peakRGB > 0);

                %trough of the gabor:
                troughrgb = inv(ptb.M)*((1-tmpLMS./100)' .* ptb.bkgndlms);
                troughRGB = round(ptb.maxdac .* troughrgb) + 1;
                troughInGamut = all(troughRGB < (ptb.maxdac+1)) && all(troughRGB > 0);

                if (peakInGamut && troughInGamut)
                    guess = guess + incSize;
                else
                    maxContrast(c) = guess - incSize;
                    break
                end
            end
        end

        %set up the color triplets and convert to gun intensity space
        switch stim.exptMethod
            case 1 %MoCS
                for a = 1:numColors;
                    
                    %make a scaleingVector. Don't let the LMSlowScalar go
                    %above 1 or below 0 b/c this might cause out of gamut problems
                    %(and it defeats the purpose).
                    LMSlowScalar(LMSlowScalar > 1) = 1;
                    LMSlowScalar(LMSlowScalar <= 0) = eps;
                    scaleVector = logspace(log10(LMSlowScalar(a)), log10(1), nContrasts);
                    
                    %set the high end to that specified by Rex or, if this
                    %value is out of the gammut, than set to the max
                    %possible
                    if norm(LMSfromREX(:,a)) > maxContrast(a);
                        coneLevels = (units(a,:)' .* maxContrast(a)) * scaleVector;
                    else
                        coneLevels = LMSfromREX(:,a) * scaleVector;
                    end
                    
                    % Notice that I'm explicitly setting the 1st contrast to zero!!!
                    coneLevels(:,2:nContrasts+1) = coneLevels;
                    coneLevels(:,1) = 0;
                    
                    % convert from CC to cone intensities.
                    coneLevels = ((coneLevels ./ 100) + 1) .* repmat(ptb.bkgndlms, 1, size(coneLevels, 2));

                    %lastly, convert to gun intensity space.
                    stim.rgbLevels{a} = [inv(ptb.M) * coneLevels]';
                end
            case 2 %staircasing
                %staircasing will be conducted in CC space (i.e. on
                %each trial a certain amount of CC will be added to each of
                %the cones (which means I'll need to find the rgb triplets
                %later).
                
                nContrasts = 1; %just in case it was set to something else on rex.

                %start by finding the inc and dec CC values
                decval = 0.05; 
                percentAtThresh = 0.816;
                incval = (decval .* percentAtThresh) ./ (1 - percentAtThresh);

                stim.stairStats = {}; % so that hitting reset states clears this out
                for a = 1:numColors;
                    %set up a stats matrix that keeps track of the
                    %performance on the last trial as well as the most
                    %recent LMS CC values
                    [stim.stairStats{a,1:nSptFreqs}] = deal([0, LMSfromREX(:,a)']);

                    %now determine the inc and dec sizes for each color
                    %direction such that the ratio L:M:S stays constant
                    %throughout the experiment
                    unitLMS = LMSfromREX(:,a) ./ norm(LMSfromREX(:,a));
                    stim.stairIncLMS(:,a) = unitLMS .* incval;
                    stim.stairDecLMS(:,a) = unitLMS .* decval;
                    stim.stairMaxLMS(:,a) = unitLMS .* maxContrast(a);
                end
            case 3 %QUEST
                %set up the first quest function for each color/sfs type.
                %Each type will have it's own range of contrasts given by
                %the monitor gamut for that color direction.
                
                q.colorDirs = units; %as determined at the top of the subfunction
                for c = 1:size(q.colorDirs, 1);
                    q.range{c} = 0.01:0.05:maxContrast(c);
                    % set up the priors
                    for s = 1:length(stim.spPeriodLevels);
                        prior = 1./(q.sdPrior.*sqrt(2.*pi)) .* exp(-(q.range{c} - q.muPrior(c)).^2 ./ (2.*q.sdPrior.^2));
                        q.Q{c,s} = log(prior);
                    end
                        
                end
        end

                
        %set up the trial type matrix
        switch stim.exptMethod
            case 1 %MoCS
                counters = fullfact([numColors, nSptFreqs, nContrasts, 2]); %2, one for each target
                counters(:,3) = counters(:,3) + 1; % reserve '1' for the zero contrast condition
                counters((end+1):(end+2),:) = [1,1,1,1;...
                                           1,1,1,2];  %take the zero contrast at each of the spatial locations
                stim.counters = [counters, ones(size(counters,1), 1).*stim.halfTrials]; %make a column of counters.
            case 2 %Stairs
                % staircasing is not currently supported
            case 3 %quest
                counters = fullfact([numColors, nSptFreqs, 2]); %2, one for each target
                counters(:,4) = counters(:,3); %confroming to the template mandated by MoCS
                counters(:,3) = 0; %just a place holder
                stim.counters = [counters, ones(size(counters,1), 1).*stim.halfTrials]; %make a column of counters.          
        end
    end

%
%   SETUPTRIAL
%
% This fxn accepts stimulus parameters from Rex and sets up the following
% trial according to the following inputs
%  RFpos            => a 2 vector, the (x,y) of the RF in tenths of degrees
%  flashLength      => in milliseconds
%  fpPos            => in tenths of degrees
%  fpSize           => in tenths of degrees
%  frameOffset      => tenths of degrees, distance b/w the frame & gabor
%  targSize         => tenths of degrees
%  T1 & T2Dist      => from fp in tenths of degrees
%  gaborTheta       => in radians
%  gaborGamma       => a floating point number
%  gaborSigma       => one sided standard dev in tenths of degrees
%  gaborSpeed       => temporal freq in cyc/sec
%  nSigmas          => number of SDs to present (effectively the gabor size)
%  lastTrialGood    => binary flag
%  lastTrialCorrect => binary flag used in staircasing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function setupTrial(RFpos, flashLength, fpPos, fpSize, frameOffset, targSize, T1Dist, T2Dist, gaborTheta, gaborGamma, gaborSigma, gaborSpeed, nSigmas, lastTrialGood, lastTrialCorrect, trialType)
        
        %Many of the gabor parameters are specified by rex. Convert those
        %to the appropriate units, and then select the remaining parameters
        %randomly (below)
        stim.gaborTheta = gaborTheta;
        stim.gaborSigma = round(ptb.pixperdeg .* gaborSigma ./ 10);
        stim.gaborGamma = gaborGamma;
        stim.flashDeltaPhase = gaborSpeed .* 2 .* pi .* (1 ./ ptb.refreshRate); %the amount to advance each frame
      
        %pick a color direction, contrast, and spatial frequency to present. This will
        %depend on weather your using an adaptive proceedure or method of
        %constant stimuli.
        switch stim.exptMethod
            case 1 %MoCS
                [stim.rgbTrial, stim.gaborLambda, stim.tColor, stim.tSfs, stim.tContrast, stim.tStimLoc] = MoCS(lastTrialGood, trialType);
            case 2 %Stairs
                [stim.rgbTrial, stim.gaborLambda, stim.tColor, stim.tSfs, stim.tContrast, stim.tStimLoc] = Staircase(lastTrialGood, lastTrialCorrect);
            case 3 %Quest
                [stim.rgbTrial, stim.gaborLambda, stim.tColor, stim.tSfs, stim.tContrast, stim.tStimLoc] = Quest(lastTrialGood, lastTrialCorrect);
        end
        
        %convert the rgb for each trial into DAC values
        if any(isnan(stim.rgbTrial))
            disp('End of Experiment');
            return
        end
        RGB = round(ptb.maxdac .* stim.rgbTrial) + 1;
        stim.flashRGB(1:3) = round(ptb.maxdac .* [ptb.invGamma(RGB(1), 1), ptb.invGamma(RGB(2), 2), ptb.invGamma(RGB(3), 3)]);
       
        
        %spatial location of the FP, wire frame, stimulus flash, and
        %targets. On grating trials stim.tStimLoc is initialized to an
        %empty vector (to avoid book keeping issue) but I still want to
        %present it in the RF
        if isempty(stim.tStimLoc) || (stim.tStimLoc == 1)
                stim.flashPos = RFpos;
        elseif stim.tStimLoc == 2
                stim.flashPos = -RFpos;
        end
        
        %make the stimulus drawing rectangle
        stim.fpRect = pos2rect(fpPos, fpSize);
        if trialType == 1 %grating
            stim.flashRect = pos2rect(stim.flashPos, stim.prefDiam.*10);
        else %detection
            sizeInTenthsDeg = gaborSigma .* nSigmas .* 2; %multiply by two so that +/-nSigmas mass is shown
            stim.flashRect = pos2rect(stim.flashPos, sizeInTenthsDeg);
            stim.frameRectIpsi = pos2rect(stim.flashPos, sizeInTenthsDeg+frameOffset);
            stim.frameRectContra = pos2rect((stim.flashPos.*-1), sizeInTenthsDeg+frameOffset);
            T1PosInDeg = round(RFpos .* T1Dist);
            stim.T1Rect = pos2rect(T1PosInDeg, targSize);
            T2PosInDeg = round(-RFpos .* T2Dist);
            stim.T2Rect = pos2rect(T2PosInDeg, targSize);
        end


        %deal with the time weighting function
        stim.flashNumFrames = ceil(ptb.refreshRate .* (flashLength./1000));
        rampLength = ceil(stim.flashNumFrames ./ 4);
        ramp = linspace(0, 1, rampLength);  %ramp is 1/4th of the total duration on either side
        plateau = ones(1,stim.flashNumFrames - (rampLength .* 2));
        stim.flashTimeProf = [ramp, plateau, fliplr(ramp)]; 
        
        %make the stimulus movie ahead of time.% zero this out the
        %flashFrameCounter here (and not in AllOff. Otherwise the number of
        %frames won't be droped with the trialParams.
        stim.flashFrameCounter = 0;
        switch trialType
            case 0 % detection trial
                halfSize = (stim.flashRect(3) - stim.flashRect(1)) ./ 2;  %base the size of the gabor off of the flash's texture window so that there are no mistakes
                row = -halfSize:2:halfSize-1; %subtract one so that you don't overlap the texture window;
                col = -halfSize:halfSize-1;
                [X, Y] = meshgrid(row, col);
                stim.flashMovie = makeGaborMovie(X,Y);
            case 1 %grating (fixation) trial
                stim.flashMovie = makeGratingMovie();
        end
                
    end


%
%   METHOD OF CONSTANT STIMULI CONTRAST GENERATOR
%
% Determines which stimulus to present next and decrements the appropriate
% counter based on the previous trial's results. 
%
%  rgbTrial     => the rgb triplet for the next trial
%  lambda       => spatial period in pix
%  color        => scalar index
%  sfs          => scalar index
%  contrast     => scalar index
%  targ         => (1 == ipsi to RF) (2 == contr to RF)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [rgbTrial, lambda, color, sfs, contrast, targ] = MoCS(lastTrialGood, trialType)

        %decriment the appropriate counter if the last trial was a good one
        if lastTrialGood
            [tf, idx] = ismember([stim.tColor, stim.tSfs, stim.tContrast, stim.tStimLoc], stim.counters(:,1:4), 'rows');
            stim.counters(idx, end) = stim.counters(idx, end) - 1;
        end
        
        switch trialType
            case 0 %detection trials
                [color, sfs, contrast, targ] = trialAllocator();
                if ~any([color, sfs, contrast, targ])
                    [rgbTrial, lambda] = deal(NaN); % crappy way to deal with the end of an expt
                    return
                end
                rgbTrial = stim.rgbLevels{color}(contrast,:); %the rgb triplet for the trial type
                lambda = stim.spPeriodLevels(sfs);
            case 1 % grating (fixation) trials
                tLMS = stim.gratingColor(1,1:3)+1; %in CC
                tlms = tLMS' .* ptb.bkgndlms; %in cone excitations
                rgbTrial = (inv(ptb.M) * tlms)'; %in gun intensities
                lambda = stim.prefLambda;
                color = [];   %setting these to empty vectors should obviate book keeping issues when it comes time to decriment the "appropriate" trial type on the next trial
                sfs = [];
                contrast = [];
                targ = [];
        end
        
    end


%
%     STAIRCASING CONTRAST GENERATOR
%
% Chooses the parameters for the next stimulus during staircasing
% experiments. This function also decrements the appropriate counters based
% on the previous trial's outcome and keeps track of the subjects
% performance for each of the stimulus types
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        rgbTrial = [inv(ptb.M) * ((1+(trialLMS'./100)) .* ptb.bkgndlms)]';
        indicies = [color, sfs, contrast, targ];
        lambda = stim.spPeriodLevels(sfs);
    end


%
%    QUEST TRIAL GENERATOR
%
% This function updates quest functions (stored seperately for each
% color/sfs combo) on each trial. It also selects the next trial's contrast
% according to the modal value of the quest function. The intial values of
% the 'output' arguments (which are nested globals) reflect the previous
% trial's params. Update the appropriate quest function then draw new
% parameters using the trialAllocator function. See notes from MoCS for
% conventions on output argument types.
%
% q.range       => range of intensities in CC (0to100%, not 0to1)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [rgbTrial, lambda, color, sfs, contrast, targ] = Quest(lastTrialGood, lastTrialCorrect);
        %deal with the previous trial's counter and quest functions
        if lastTrialGood
            %update the counter
            [t, loc] = ismember([stim.tColor, stim.tSfs, 0, stim.tStimLoc], stim.counters(:,1:4), 'rows');
            stim.counters(loc, end) = stim.counters(loc, end) -1;
            
            %update the quest function
            tConeContrast = norm(((ptb.M*stim.rgbTrial(:))-ptb.bkgndlms) ./ ptb.bkgndlms);
            tConeContrast = tConeContrast.*100; %convert to %CC 1to100;
            liklihood = 1-0.5.*exp(-(tConeContrast./q.range{stim.tColor}).^q.modBeta);
            if ~lastTrialCorrect
                liklihood = 1-liklihood;
            end
            q.Q{stim.tColor, stim.tSfs} = q.Q{stim.tColor, stim.tSfs} + log(liklihood);
        end

        %obtain stimulus parameters for the upcomming trial
        [color, sfs, contrast, targ] = trialAllocator();
        if ~any([color, sfs, contrast, targ])
            [rgbTrial, lambda] = deal(NaN); % crappy way to deal with the end of an expt
            return
        end
        
        [val, idx] = max(q.Q{color, sfs});
        coneContrast = q.range{color}(idx);
        LMS = (q.colorDirs(color,:) .* coneContrast)./100; % in %CC (0to1)
        rgbTrial = [inv(ptb.M) * ((1+LMS)' .* ptb.bkgndlms)]';
        lambda = stim.spPeriodLevels(sfs);
    end

%
%   TRIAL TYPE ALLOCATOR
%
% This function randomly selects a stimulus type to present on the next
% trial and keeps track of how many blocks of trials have been presented.
% When all trials from all blocks have been presented this function returns
% NANs, which are interpreted by the trial setup function to signify the
% end of an experiment
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%
%  SET QUEST PARAMETERS
%
% Called only when exptMeth == 3. Simply sets the local quest parameters
% based on those specified by REX. This function gets called before
% setupColorContrasts. 
%  priorSD      => scalar in units of CC
%  beta         => scalar; slope of the models psyFun
%  alpha_x      => vector; guesses of alpha for each color/sfs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function setQuestParams(priorSD, beta, alpha_1, alpha_2, alpha_3);
        q.muPrior = [alpha_1; alpha_2; alpha_3]; %only one estimate per color
        q.sdPrior = priorSD;
        q.modBeta = beta;
    end


%
%  MONITOR SETUP
%
% Computes all the relevant monitor information including the gamma table,
% monitor spectra and M matrix.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function setupMonitor(mondist, screenwidth, calFilePath)

        %load in the calibration data and fill up the ptb structure
        S = load('T_cones_smj10.mat');
        load(calFilePath)
        calData = cals{end};
        fundamentals = S.T_cones_smj;
        ptb.bkgndRGB = round(255 .* calData.bgColor); %bkgnd voltages discritized b/w 0&255
        ptb.bkgndrgb = [calData.gammaTable(ptb.bkgndRGB(1)+1, 1), calData.gammaTable(ptb.bkgndRGB(2)+1, 2), ...
            calData.gammaTable(ptb.bkgndRGB(3)+1, 3)]; %add one to create and index b/w 1 & 256
        stim.targRGB = ptb.bkgndRGB - 30;


        ptb.gammaTable = calData.gammaTable;
        ptb.monSpd = SplineSpd(calData.S_device, calData.P_device, S.S_cones_smj);
        ptb.M = fundamentals * ptb.monSpd;
        ptb.invGamma = InvertGamma(calData, 1);
        ptb.bkgndlms = ptb.M * ptb.bkgndrgb';

        %open a ptb Screen window and determine various conversion factors.
        %Don't open window if MasterSlave has already done so
        if isempty(Screen('windows'));
            ptb.w = Screen('OpenWindow', 0, ptb.bkgndRGB);
        else
            ptb.w = Screen('windows');
            Screen('FillRect', ptb.w, ptb.bkgndRGB);
        end
        
            
        HideCursor;
        [ptb.monWidth, ptb.monHeight] = Screen('WindowSize', ptb.w);
        theta = atand((screenwidth./2) ./ mondist);
        ptb.pixperdeg = (ptb.monWidth ./ 2) ./ theta;
        ptb.oldGamma = Screen('LoadNormalizedGammaTable', ptb.w, repmat(linspace(0,1,256)', 1, 3));
        ptb.refreshRate = Screen('NominalFrameRate', ptb.w, 1);
        ptb.monCent = [ptb.monWidth, ptb.monHeight] ./ 2; %(x,y) at the center of the screen
        
        %inform REX that the monitor setup is complete
        pnet(udpCom.sock, 'write', 'MonitorSetupComplete>> >>');
        pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
    end

%
%          MAKE THE STIMULUS MOVIE
%
% This function makes all frames of the stimulus so that the timing of
% gabor presentation can be tightly controlled. It also keeps track of the
% exptected and actual cone contrast of the gabor 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function movie = makeGaborMovie(X,Y);
        DT_movie = {};
        xprime = X .* cos(-stim.gaborTheta) + Y .* sin(-stim.gaborTheta);
        yprime = -X .* sin(-stim.gaborTheta) + Y .* cos(-stim.gaborTheta);
        rgb_increment = stim.rgbTrial - ptb.bkgndrgb;
        
        for a = 1:stim.flashNumFrames;

            %start by making the gabor, then multiplying by the increment from
            %background and by the temporal weighting fxn.
            gabor = exp(-(xprime.^2 + stim.gaborGamma.^2 .* yprime.^2) ./ (2.* stim.gaborSigma.^2)) .* cos(2 .* pi .* yprime ./ stim.gaborLambda + stim.gaborPhase);
            flashImg = repmat(gabor, [1, 1, 3]);
            flashImg(:,:,1) = (flashImg(:,:,1) .* stim.flashTimeProf(a) .* rgb_increment(1)) + ptb.bkgndrgb(1);
            flashImg(:,:,2) = (flashImg(:,:,2) .* stim.flashTimeProf(a) .* rgb_increment(2)) + ptb.bkgndrgb(2);
            flashImg(:,:,3) = (flashImg(:,:,3) .* stim.flashTimeProf(a) .* rgb_increment(3)) + ptb.bkgndrgb(3);

            %now convert to DAC values
            imgSize = size(gabor);
            flashrgb = round(ptb.maxdac .* flashImg) + 1; %intensities b/w 1 & maxdac +1
            flashRGB = ones(size(flashImg)); %preallocate so that indexing works
            flashRGB(:,:,1) = reshape(round(ptb.maxdac .* ptb.invGamma(flashrgb(:,:,1), 1)), imgSize(1), imgSize(2));
            flashRGB(:,:,2) = reshape(round(ptb.maxdac .* ptb.invGamma(flashrgb(:,:,2), 2)), imgSize(1), imgSize(2));
            flashRGB(:,:,3) = reshape(round(ptb.maxdac .* ptb.invGamma(flashrgb(:,:,3), 3)), imgSize(1), imgSize(2));

            
            %find the index to the pixel with the largest contrast. at the end
            %of the trial, the RGB for this pixel will get dropped in the data
            %file.
            norms = sqrt(sum([flashImg .* flashImg], 3)); %in normalized gun space
            frameMaxInc = max(norms(:));
            if frameMaxInc > stim.flashMaxInc;
                [r, c] = find(norms == frameMaxInc, 1);
                stim.flashRGB(4) = flashRGB(r,c,1);  % DAC values b/w 0&maxdac
                stim.flashRGB(5) = flashRGB(r,c,2);
                stim.flashRGB(6) = flashRGB(r,c,3);
                stim.flashMaxInc = frameMaxInc;
            end


            %make a frame of the stimulus movie
            movie(a) = Screen('MakeTexture', ptb.w, TranslateToColourModeMex(flashRGB));
            DT_movie{length(DT_movie)+1} = flashRGB;
            %update the position of the gabor for the next frame refresh
            stim.gaborPhase = stim.gaborPhase + stim.flashDeltaPhase;
        end
    end

%
%   MAKE A GRATING MOVIE
%
%   On certain catch trials (controled by the user through the DTonline
%   GUI) a drifting sinusoidal grating is present instead of a gabor. This
%   allows the user to see if a neuron is still present despite a very low
%   baseline firing rate or responsiveness to a stimulus
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function movie = makeGratingMovie()
        pixpercycle = ptb.pixperdeg/(1./stim.prefLambda);
        sizeinpix = (stim.flashRect(3) - stim.flashRect(1));
        phase = stim.prefPhase;
        
        % Creating an aperture template
        [x,y] = meshgrid(linspace(-1,1,sizeinpix./2), linspace(-1,1,2*sizeinpix./2));
        % Only half the number of columns as rows since we're using colour mode
        aperture = (x.^2+y.^2) <= 1+2/sizeinpix;  % correction for binning
        xinc = (2*pi/pixpercycle)*cos((pi/2)-stim.prefTheta);
        yinc = (2*pi/pixpercycle)*sin((pi/2)-stim.prefTheta);
        [xramp, yramp] = meshgrid(xinc*([0:2:sizeinpix-1]), yinc*([0:sizeinpix-1]));
        
        %make the movie
        for i = 1:stim.flashNumFrames
            a = cos(xramp+yramp+phase);
            a = a.*aperture;
            im = zeros(sizeinpix, sizeinpix/2, 3);
            for plane = 1:3
                tmp = a.*(stim.rgbTrial(plane)-ptb.bkgndrgb(plane))+ptb.bkgndrgb(plane);
                tmp = round(tmp*ptb.maxdac)+1;
                tmp = ptb.invGamma(tmp, plane);
                tmp = round(tmp*ptb.maxdac);
                im(:,:,plane) = reshape(tmp, sizeinpix, sizeinpix/2);
            end

            phase = phase + stim.flashDeltaPhase;
            movie(i) = Screen('MakeTexture', ptb.w, TranslateToColourMode(im,1));
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
        stim.flashMaxInc = 0; %for determining differences b/w expected and presented cone contrasts
        stim.flashMovie = [];
        Screen('Close') %I think that this closes/deletes all offscreen textures
    end


%
%  SEND QUEST STATS TO REX
%
%   this fxn only gets called by REX druing quest experiments and sends
%   accross the modal value of the appropriate posterior distribution. Due
%   to the current architecture of Rex's state set (and how I update the
%   quest stats) this number actually reflects the modal value of the most
%   recent trial NOT the current trial. The current trials results won't be
%   known until 'setupTrial' gets called for the next trial.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function updateQuestStats()
        [m, idx] = max(q.Q{stim.tColor, stim.tSfs});
        pThreshMode = q.range{stim.tColor}(idx);
        
        %use send to rex to transmit the modal value of the posterior
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
        msg = zeros(1,9);
        nColors = length(q.range);
        for a = 1:nColors;
            msg((a.*3)-2) = min(q.range{a});
            msg((a.*3)-1) = max(q.range{a});
            msg(a.*3) = q.range{a}(2) - q.range{a}(1); %pretty kludgy
        end
        
        %now send the message using sendToRex
        sendToRex(udpCom, msg, 'double', 'sendQuestColorRanges();');
        
    end



%
%          MAKE A DRAWING RECTANGLE
%
% draws a square of side length "size" (in pix) arround a point on the
% monitor specified by centPos (in (x,y) cordinates) size is in units of
% tenths of degrees
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function rect = pos2rect(centPos, size)

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


%
%    GET/SET GABOR GLOBALS
%
% This functions loads the gabor/grating fits into the slave program. At
% the end of the grating or WN paradigm a global variable (gl) exists on
% the slave machine that contains all the information regarding receptive
% field parameters fit by these two programs. This function simply extracts
% the relvant RF data and converts to units suitable to the DT paradigm
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function loadGaborFits(fitMethod);
        %fitMethod (1) = Grating with prefIsolum colors, (2) = WN online
        %(3) = Grating with S- and L-M colors.

        global gl;

        %pull in the params and convert to units that REX (not the slave)
        %uses. The logic here is that REX feeds these parameters to the Mac
        %on every trial. So I'll send these to REX at the beginning of a
        %run, and REX will send them back b/4 each trial.
        if find([1,3] == fitMethod)
                if ~isfield(gl, 'grating');
                    return
                end

                if gl.grating.prefsf %prevent Infs when prefsf = 0
                    stim.prefLambda = (1./gl.grating.prefsf);               % deg/cyc
                else
                    stim.prefLambda = eps; %don't let the sfs = 0
                end
                stim.prefTheta = gl.grating.preforient;                     % radians
                stim.prefPhase = 0;
                stim.prefSigma = 4;                                         % hard coding the sigma !!!
                stim.prefGamma = 1;                                         % no other choice unless we modle the WN data
                stim.prefDiam = gl.grating.prefdiam;
                stim.RFPos = [gl.grating.x, gl.grating.y] .* 10;             % comes in dva. local standard is 10ths dva

                %if Grating.d has determined a prefered color compute the
                %appropriate color contrast ranges for an MoCS run.
                if isfield(gl.grating, 'prefIsolum')
                    switch fitMethod
                        case 1
                            exptColor = reshape(gl.grating.prefIsolum, 3, 2)';
                        case 3
                            exptColor = [1 -1 0;  0 0 1];
                    end
                    
                    load sCSF_Fits_Kali.mat
                    for a = 1:2;
                        [val, idx] = ismember(sign(exptColor(a,:)), sign(fp.colors), 'rows');
                        sensitivity = polyval(fp.polyCoeff(idx,:), gl.grating.prefsf);

                        %fp.colors contains unit vectors describing each color
                        %direction. I'll use these unit vectors to establish
                        %the LMS triplets for Rex to use. The scale factor for
                        %each of these color directions will be 0.25. This will
                        %be enforced on the Rex side. prefColors must be a %
                        %b/w 0 and 100 (not 0 & 1) so mulitiply the sensitivity
                        %by 100;
                        stim.exptColorsForRex(1, ((a-1)*3)+1 : a*3) = 3.2 .* (1./sensitivity) .* fp.colors(idx, :) .* 100;
                    end
                    
                    %make sure that the relative phases of the three cone types
                    %is consistent betweent the chosen cardinal and
                    %intermediate color
                    if (stim.exptColorsForRex(1:3) * stim.exptColorsForRex(4:6)') < 0;
                        stim.exptColorsForRex(1:3) = stim.exptColorsForRex(1:3) .* -1;
                    end
                    
                    %save the color that maximally excited the neuron. I'll
                    %use this color for grating catch trials.
                    stim.gratingColor = gl.grating.prefccdir;
                end

        elseif fitMethod == 2 %fits to white noise
                if ~isfield(gl, 'gabor');
                    return
                end
                stim.prefTheta = gl.gabor.theta;                           % radians
                stim.prefLambda = gl.gabor.lambda;                         % deg/cyc
                stim.prefPhase = 0;
                stim.prefSigma = gl.gabor.sigma .* 10;                     % DTdefault is 10ths of degrees but WN fits with degrees
                stim.prefGamma = gl.gabor.gamma;                           % 1 gives circles, else elipses

                %now determine the X,Y for the stimulus center. I'm going to send
                %these values to REX, so they need to be in tenths of degrees.
                x = [gl.stim.x + gl.gabor.xoffset]  .* 10; %these globals are in dva, so multiply by 10 to obtain
                y = [gl.stim.y + gl.gabor.yoffset]  .* 10;  %the local standard of 10ths of dva
                stim.RFPos = [x,y]; %needs to be in tenths of degrees for use locally
        end
    end


%HANDLE UDP COMMUNICATIONS
    function dealWithUdpMessages()

        msgSize = pnet(udpCom.sock, 'readpacket', 500, 'noblock');
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
            wtf %print out the error message
        end
    end

end %DTslave


%%
%******* NOTES AND QUESTIONS **************
% out or bounds errors occur when the LMS triplet is a NaN. The occurs b/c
% rex desires O CC conditions (log(0) is a NaN).
