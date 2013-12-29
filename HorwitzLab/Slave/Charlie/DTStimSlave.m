function DTStimSlave

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
% CAH 05/11 Striped all the functionality for staircasing and incorporated
%           functionality for optical stimulation. This update is a
%           significant departure from the regular version of DTslave

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
stim.gaborSigmaInTenthsDeg = [];
stim.gaborGamma = 0;            % a float; 1 give circles, else elipses
stim.gratingColor = [];         % the color for the grating catch trials (color of highest firing rate in grating.d)
stim.halfTrials = [];           % the number of trials per block per target i.e 0.5(ntrials per block)
stim.fpRect = [];               % drawing rectangle for the FP
stim.flashPos = [];
stim.nSigmas = [];
stim.RFPos = [0 0];             % this gets set by loadGaborFits
stim.rgbLevels = [];            % the rgb triplets for each stimulus type
stim.spPeriodLevels = [];       % a vector of possible spatial periods (specifed by rex as spt. freqs and converted below)
stim.T1Rect = [];               
stim.T2Rect = [];
stim.targRGB = [];              % gets set during monior setup.
stim.targSize = [];             % set by rex (in tenths dva)
stim.T1Dist = [];               %fractional distance to the stimulus from the FP
stim.T2Dist = [];

%TRIAL INFO
trial.drawFP = 0;
trial.drawFlash = 0;
trial.drawTarg = 0;
trial.sentTargMsg = 0;          % a flag to notify rex that the targs are now on.


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
        
        if(trial.drawTarg)
            Screen('FillRect', ptb.w, stim.targRGB, stim.T1Rect);
            Screen('FillRect', ptb.w, stim.targRGB, stim.T2Rect);
            if ~trial.sentTargMsg
                pnet(udpCom.sock, 'write', 'MACtargson>> >>');
                pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
                trial.sentTargMsg = 1;
            end
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
        
        if(trial.drawFP) %putting this last so that it is printed on top of everything else.
            Screen('FillRect', ptb.w, [0,0,0], stim.fpRect);
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
%  gaborLength      => in ms
%  gaborGamma       => the ratio of SDx:SDy
%  gaborSigma       => in tenths DVA 
%  gaborTheta       => in radians
%  gaborSpeed       => in cyc/sec
%  nSigmas          => +/- this many sigmas to truncate at
%  RFpos            => 2 vector (the x,y of the RF)
%  fixPos           => 2 vector (x,y of the FP)
%  fixPosSize       => in 1/10 DVA
%  targSize         => in 1/10 DVA
%  targDist1_2      => 2 vector ([T1, T2]), fractional dist to targ
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function setupExpConstants(LMSfromREX, ...
                                LMSlowScalar, ...
                                sptFreqfromREX, ...
                                nSptFreqs, ...
                                nContrasts, ...
                                nTrialsPerBlock, ...
                                nBlocks, ...
                                gaborLength, ...
                                gaborGamma, ...
                                gaborSigma, ...
                                gaborTheta, ...
                                gaborSpeed, ...
                                nSigmas, ...
                                RFpos, ...
                                fixPos, ...
                                fixPosSize, ...
                                targSize, ...
                                targDist1_2 ...
                                )
                            
        % initialize these parameters and experimental constants.
        stim.RFpos = RFpos;
        stim.targSize = targSize;
        stim.T1Dist = targDist1_2(1);
        stim.T2Dist = targDist1_2(2);
        stim.gaborTheta = gaborTheta;
        stim.gaborSigma = round(ptb.pixperdeg .* gaborSigma ./ 10);
        stim.gaborSigmaInTenthsDeg = gaborSigma;
        stim.gaborGamma = gaborGamma;
        stim.nSigmas = nSigmas;
        stim.flashDeltaPhase = gaborSpeed .* 2 .* pi .* (1 ./ ptb.refreshRate); %the amount to advance each frame
        stim.exptDone = 0;
        stim.counters = [];
        stim.blockCounter = nBlocks;
        stim.halfTrials = ceil(nTrialsPerBlock./2); %so that trials get allocated to each target evenly
        
        %deal with the time weighting function
        stim.flashNumFrames = ceil(ptb.refreshRate .* (gaborLength./1000));
        rampLength = ceil(stim.flashNumFrames ./ 4);
        ramp = linspace(0, 1, rampLength);  %ramp is 1/4th of the total duration on either side
        plateau = ones(1,stim.flashNumFrames - (rampLength .* 2));
        stim.flashTimeProf = [ramp, plateau, fliplr(ramp)];
        
        % a drawing rectangle for the fp
        stim.fpRect = pos2rect(fixPos, fixPosSize);

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
        noColor = LMSlowScalar == 0;
        if any(noColor)
            LMSfromREX(:,noColor) = []; %delete LMS triplets that are all zeros
            LMSlowScalar(noColor) = []; %don't forget to delete the scalars
        end
        numColors = size(LMSfromREX,2);
        
        
        %determine the extent of the monitor gamut for the specified color
        %directions
        norms = sqrt(diag(LMSfromREX' * LMSfromREX));
        units = LMSfromREX' ./ repmat(norms, 1, 3); %triplets go across rows!
        for c = 1:size(units, 1);
            if any(isnan(units(c,:)))
                maxContrast(c) = inf; %this happens only when LMSfromRex is [0 0 0]
            else
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
        end

        %setup the RGB triplets that define each color and contrast
        for a = 1:numColors;
            %make a scalingVector. Don't let the LMSlowScalar go
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
            if any(isnan(stim.rgbLevels{a})); keyboard; end
        end
                
        %set up the trial type matrix
        counters = fullfact([numColors, nSptFreqs, nContrasts, 2]); %2, one for each target
        counters(:,3) = counters(:,3) + 1; % reserve '1' for the zero contrast condition
        counters((end+1):(end+2),:) = [1,1,1,1;...
            1,1,1,2];  %take the zero contrast at each of the spatial locations
        stim.counters = [counters, ones(size(counters,1), 1).*stim.halfTrials]; %make a column of counters.
        
        % alert REX that contant setup is complete
        pnet(udpCom.sock, 'write', 'MACconstantSetupComplete>> >>');
        pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
    end

%
%   SETUPTRIAL
%
% This fxn accepts stimulus parameters from Rex and sets up the following
% trial according to the following inputs
%  lastTrialGood    => binary flag
%  optoStimOnly     => fixation trial: opto-stim only
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function setupTrial(lastTrialGood, laserCatchTrial)
        
        %pick a color direction, contrast, and spatial frequency to present.
        [stim.rgbTrial, stim.gaborLambda, stim.tColor, stim.tSfs, stim.tContrast, stim.tStimLoc] = MoCS(lastTrialGood, laserCatchTrial);
        
        %convert the rgb for each trial into DAC values
        if any(isnan(stim.rgbTrial))
            disp('End of Experiment');
            return
        end
        RGB = round(ptb.maxdac .* stim.rgbTrial) + 1;
        stim.flashRGB(1:3) = round(ptb.maxdac .* [ptb.invGamma(RGB(1), 1), ptb.invGamma(RGB(2), 2), ptb.invGamma(RGB(3), 3)]);
       
        
        %spatial location: On fixation trials stim.tStimLoc is initialized
        %to an empty vector but I still need to put it somewhere (present
        %the blank gabor in the T2 location)
        if stim.tStimLoc == 1
                stim.flashPos = stim.RFpos;
        elseif isempty(stim.tStimLoc) || (stim.tStimLoc == 2)
                stim.flashPos = -stim.RFpos;
        end
        
        %make some drawing rectangles        
        sizeInTenthsDeg = stim.gaborSigmaInTenthsDeg .* stim.nSigmas .* 2; %multiply by two so that +/-nSigmas mass is shown
        stim.flashRect = pos2rect(stim.flashPos, sizeInTenthsDeg);
        T1PosInDeg = round(stim.flashPos .* stim.T1Dist);
        stim.T1Rect = pos2rect(T1PosInDeg, stim.targSize);
        T2PosInDeg = round(-stim.flashPos .* stim.T2Dist);
        stim.T2Rect = pos2rect(T2PosInDeg, stim.targSize);
 
        
        %make the stimulus movie ahead of time.% zero this out the
        %flashFrameCounter here (and not in AllOff. Otherwise the number of
        %frames won't be droped with the trialParams.
        stim.flashFrameCounter = 0;
        halfSize = (stim.flashRect(3) - stim.flashRect(1)) ./ 2;  %base the size of the gabor off of the flash's texture window so that there are no mistakes
        row = -halfSize:2:halfSize-1; %subtract one so that you don't overlap the texture window;
        col = -halfSize:halfSize-1;
        [X, Y] = meshgrid(row, col);
        stim.flashMovie = makeGaborMovie(X,Y);
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
            case 1 % fixation trials (opto-stim only)
                rgbTrial = ptb.bkgndrgb; %in gun intensities
                lambda = stim.spPeriodLevels(1); %taking the first one by default...
                color = [];   %setting these to empty vectors should obviate book keeping issues when it comes time to decriment the "appropriate" trial type on the next trial
                sfs = [];
                contrast = [];
                targ = [];
        end
        
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
%  MONITOR SETUP
%
% Computes all the relevant monitor information including the gamma table,
% monitor spectra and M matrix.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function setupMonitor(mondist, screenwidth, calFilePath, targContrastModifier)

        %load in the calibration data and fill up the ptb structure
        load(calFilePath)
        calData = cals{end};
        smj = load('T_cones_smj.mat');
        fundamentals = smj.T_cones_smj;
        fundWavelengthSpacing = smj.S_cones_smj;
        ptb.bkgndRGB = round(255 .* calData.bgColor); %bkgnd voltages discritized b/w 0&255
        ptb.bkgndrgb = [calData.gammaTable(ptb.bkgndRGB(1)+1, 1), calData.gammaTable(ptb.bkgndRGB(2)+1, 2), ...
            calData.gammaTable(ptb.bkgndRGB(3)+1, 3)]; %add one to create and index b/w 1 & 256
        targContrastModifier = min(min(ptb.bkgndRGB), targContrastModifier);
        stim.targRGB = ptb.bkgndRGB - targContrastModifier;

        ptb.gammaTable = calData.gammaTable;
        ptb.monSpd = SplineSpd(calData.S_device, calData.P_device, fundWavelengthSpacing);
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
    function movie = makeGaborMovie(X,Y)        
        xprime = X .* cos(-stim.gaborTheta) + Y .* sin(-stim.gaborTheta);
        yprime = -X .* sin(-stim.gaborTheta) + Y .* cos(-stim.gaborTheta);
        rgb_increment = stim.rgbTrial - ptb.bkgndrgb;
        imgSize = size(xprime);
        flashImg = nan(size(xprime,1), size(xprime,2), 3);
        flashrgb = nan(size(xprime,1), size(xprime,2), 3);
        flashRGB = nan(size(xprime,1), size(xprime,2), 3);

        for a = 1:stim.flashNumFrames;

            %start by making the gabor, then multiplying by the increment from
            %background and by the temporal weighting fxn.
            gabor = exp(-(xprime.^2 + stim.gaborGamma.^2 .* yprime.^2) ./ (2.* stim.gaborSigma.^2)) .* cos(2 .* pi .* yprime ./ stim.gaborLambda + stim.gaborPhase);
            flashImg(:,:,1) = (gabor .* stim.flashTimeProf(a) .* rgb_increment(1)) + ptb.bkgndrgb(1);
            flashImg(:,:,2) = (gabor .* stim.flashTimeProf(a) .* rgb_increment(2)) + ptb.bkgndrgb(2);
            flashImg(:,:,3) = (gabor .* stim.flashTimeProf(a) .* rgb_increment(3)) + ptb.bkgndrgb(3);

            %now convert to DAC values
            flashrgb = round(ptb.maxdac .* flashImg) + 1; %intensities b/w 1 & maxdac +1
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
            %TranslateToColourModeMex(flashRGB);
            movie(a) = Screen('MakeTexture', ptb.w, TranslateToColourModeMex(flashRGB));
            %update the position of the gabor for the next frame refresh
            stim.gaborPhase = stim.gaborPhase + stim.flashDeltaPhase;
        end
    end


%**************************%
% ACCESSORY FUNCTIONS      %
%**************************%
%BLANK THE SCREEN
    function allOff()

        trial.drawFP = 0;
        trial.drawFlash = 0;
        trial.drawTarg = 0;

        %reset these gloabal parameters for the next trial
        stim.gaborPhase = 0;
        stim.flashMaxInc = 0; %for determining differences b/w expected and presented cone contrasts
        stim.flashMovie = [];
        trial.sentTargMsg = 0;
        Screen('Close') %I think that this closes/deletes all offscreen textures
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


%HANDLE UDP COMMUNICATIONS
    function dealWithUdpMessages()

        msgSize = pnet(udpCom.sock, 'readpacket', 500, 'noblock');
        if ~msgSize
            return
        end

        try
            %open the message
            message = pnet(udpCom.sock, 'read', msgSize, 'char');
            %disp(message);
            
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
