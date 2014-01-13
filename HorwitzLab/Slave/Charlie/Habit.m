function Habit
%
% SLAVE PROGRAM FOR CHROMATIC HABITUATION EXPERIMENT
%
% This tests detection thresholds with and without prior habituation.
% Habituation can only be along a single color direction. Up to 4 stimulus
% colors can be interleaved. Any number of spatial frequencies may be
% interleaved. Thresholds will be tested by the QUEST adaptive proceedure.
%
%
% CAH 12/2009  => Starting development.

%****** CONSTANTS AND PERTINENT VARIABLES ********
%UDP COMMUNICATION INFO
udp.sock = [];
udp.port = 6665;
udp.rexip = '192.168.1.120';
udp.plexip = '192.168.1.121';

%PTB CONSTANTS
ptb.bkgndlms = [];
ptb.bkgndRGB = [];
ptb.bkgndrgb = [];
ptb.exitNow = 0;
ptb.gammaTable = [];
ptb.invGamma = [];
ptb.maxdac = (2^16)-1;          %65536 dac values in colour mode
ptb.monCent = [];
ptb.monSpd = [];
ptb.monWidth = [];
ptb.monHeight = [];
ptb.M = [];
ptb.oldGamma = [];
ptb.pixperdeg = [];
ptb.refreshRate = [];
ptb.w = [];

% QUEST PARAMETERS
q.colorDirs = [];               % unit vectors for the colors tested
q.sdPrior = [];                 % stand dev common to all priors of thresholds
q.modBeta = [];                 % slope of psy fun common to all stim conditions
q.Q = {};                       % quest functions (one for each stim condition)
q.domain = {};                  % the domain over which each quest function is evaluated (in [start, end, dContrast] form... one for each color)
q.domainForRex = [];            % domain in a format that rex expects.

%TRIAL FLAGS
trial.exptDone = 0;
trial.drawFP = 0;
trial.drawTarg = 0;
trial.drawGabor = 0;
trial.drawHabit = 0;
trial.firstTrial = 1;
trial.lastTrialCorrect = [];    %for debugging

%GABOR PARAMETERS
gabor.blockCounter = [];
gabor.frameCounter = 0;
gabor.halfTrials = [];          %gets called each time a new block is initialized
gabor.gamma = [];               % ratio of gabor SDx to SDy
gabor.lambdas = [];
gabor.length = [];              % in milliseconds. stimulus duration
gabor.movie = [];
gabor.nSigmas = [];             %number of SDs to show (effectively the gabor size)
gabor.numFrames = [];           %number of frames the gabor should be on the screen
gabor.rect = [];
gabor.sigma = [];
gabor.speed = [];               % in cyc/sec
gabor.T1pos = [];               % a two vector (x,y). The "home" position of the gabor in pixels
gabor.theta = [];               % in radians
gabor.tTypes = [];              %a matrix [<color> <sfs> <targ> <counter>] of indicies to color and sfs types
gabor.timeProf = [];            % the temporal weighting function
gabor.tLoc = [];                % index to the current trial's location
gabor.trgb = [];                % the rgbs for the current trial
gabor.tRGB = [];                % a 6vector: [expectedRGB, actualRGB] of the current trial
gabor.tLambda = [];             % current trial's spatial period (in pix)
gabor.tColor = [];              % index of current trial's color
gabor.tSfs =[];                 % index of current trial's sfs

%TARGET PARAMETERS
targ.RGB = [];
targ.T1rect = [];
targ.T2rect = [];

%HABITUATION STIM PARAMETERS:
habit.centGrayRect = [];        %set but not currently used. 
habit.frameCounter = [];
habit.driftRate = [];
habit.length = [];              %length of presentation in ms
habit.movie = [];
habit.nFramesPerCycle = [];     %the number of movie frames in a single cycle of the drifting grating
habit.nFramesToPresent = [];    %the number of frames to present during a habituation period

%FIXATION PT PARAMETERS
fp.RGB = [0 0 0];
fp.rect = [];

%******* MAIN EXPERIMENTAL CODE (nested subfunctions to follow) ************%



% open a udp socket
[udp.sock, success] = pnetStart(udp.port);
if ~success;
    error('udp initialization failed');
end
pnet(udp.sock, 'setreadtimeout', 0);
pnet(udp.sock, 'setwritetimeout', 0);


%main loop
KbName('UnifyKeyNames');
escapeKey = KbName('ESCAPE'); %platform independent definition.
kbKey = KbName('-');
while(~ptb.exitNow)
    %check for udp messages
    dealWithUdpMessages();
    
    %draw a stimulus if need be (but only if a window is open)
    if ptb.w
        if(trial.drawTarg); drawTrialElement('targets'); end
        if(trial.drawGabor); drawTrialElement('gabor'); end
        if(trial.drawHabit); drawTrialElement('habit'); end
        if(trial.drawFP); drawTrialElement('fixationPoint'); end %at the end so that fp draws on top of everything else
        Screen('Flip', ptb.w);
        
        % send the 'last frame' message after the last flip. otherwise the
        % time of stimulus presentation will be misestimated by Rex.
        if trial.drawGabor && (gabor.frameCounter == gabor.numFrames)
            pnet(udp.sock, 'write', 'MACgaborLastFrame>> >>');
            pnet(udp.sock, 'writepacket', udp.rexip, udp.port);
            trial.drawGabor = 0;
        end
        if trial.drawHabit && (habit.frameCounter == habit.nFramesToPresent)
            pnet(udp.sock, 'write', 'MAChabitLastFrame>> >>');
            pnet(udp.sock, 'writepacket', udp.rexip, udp.port);
            trial.drawHabit = 0;
        end
    end
    
    %press ESC to bail out of the endless loop
    [dum1, dum2,keycode] = KbCheck();
    if (keycode(escapeKey)) %if you press escape
        if ptb.w
            pnet(udp.sock, 'close');
            Screen('LoadNormalizedGammaTable', ptb.w, ptb.oldGamma);
            Screen('CloseAll');
            ShowCursor;
        end
        return
    elseif (keycode(kbKey))
        sca
        keyboard
    end
end %main while loop

%******************************************************
%
%            NESTED SUPPORT FUNCTIONS
%
%******************************************************
    %
    %   TOGGLE STIMULUS ELEMENTS
    %
    % This function accepts commands from Rex and turns various trial
    % elements on and off. These binary flags are used in the main
    % experiemntal loop to determine which trial elements to draw to the
    % video buffer (which happens in drawTrialElement)
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function toggleTrialElement(element, onOff)
        switch lower(element)
            case 'fixationpoint'
                trial.drawFP = onOff;
            case 'gabor'
                trial.drawGabor = onOff;
            case 'targets'
                trial.drawTarg = onOff;
                trial.drawFP = 0;
            case 'habituate'
                trial.drawHabit = onOff;
            case 'alloff'
                trial.drawFP = 0;
                trial.drawGabor = 0;
                trial.drawTarg = 0;
                trial.drawHabit = 0;
        end     
    end

    %
    %  DRAW TRIAL ELEMENT
    %
    % A general function that gets called by the main while loop and
    % presents various parts of each trial. In particular, this function
    % just draws elements to the video buffer.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function drawTrialElement(trialElement)
        switch trialElement
            case 'fixationPoint'
                Screen('FillRect', ptb.w, fp.RGB, fp.rect);
            case 'targets'
                Screen('FillRect', ptb.w, targ.RGB, targ.T1rect);
                Screen('FillRect', ptb.w, targ.RGB, targ.T2rect);
            case 'gabor'
                gabor.frameCounter = gabor.frameCounter+1; %increment the counter.
                Screen('DrawTexture', ptb.w, gabor.movie(gabor.frameCounter), [], gabor.rect, [], 0);
                Screen('Close', gabor.movie(gabor.frameCounter));
                
                % send the '1st frame' msg before the flip. if you wait till
                % after the flip the time of the stimulus presentation will be
                % misrepresented
                if (gabor.frameCounter == 1);
                    pnet(udp.sock, 'write', 'MACgaborFirstFrame>> >>');
                    pnet(udp.sock, 'writepacket', udp.rexip, udp.port);
                end
            case 'habit'
                %setting habit time to zero is a pathological case. To
                %avoid errors, set the frame counter to the number of
                %frames to present and return from the subfunction.
                if habit.nFramesToPresent == 0;
                    habit.frameCounter = habit.nFramesToPresent;
                    return
                end
                
                %I loop over a single cycle of the habit stim, so determine
                %the actual frame number based on the global counter.
                habit.frameCounter = habit.frameCounter+1;
                frame = rem(habit.frameCounter, habit.nFramesPerCycle);
                if ~frame %rem == 0 is a pathological case
                    frame = habit.nFramesPerCycle;
                end
                Screen('DrawTexture', ptb.w, habit.movie(frame), [], [0, 0, ptb.monWidth, ptb.monHeight], [], 0);
                % send the '1st frame' msg before the flip. if you wait till
                % after the flip the time of the stimulus presentation will be
                % misrepresented
                if (habit.frameCounter == 1);
                    pnet(udp.sock, 'write', 'MAChabitFirstFrame>> >>');
                    pnet(udp.sock, 'writepacket', udp.rexip, udp.port);
                end
        end
    end

    %
    %         SETUP THE COLOR DIRECTIONS AND TRIAL COUNTERS
    %
    % Rex specifies all the color directions and spatial frequencies. Take
    % these parameters and convert to units used locally on the slave program.
    % Lastly, set up a counters matrix that keeps track of the trial types
    % presented. Input arguments have the following form:
    %  LMSfromREX       => 9 vector of LMS CC triplets (up to 3 dirs)
    %  alphaGuesses     => 3 vector, a threshold guess for each color direction
    %  sptFreqfromREX   => 2 vector, low and high sfs to test
    %  nSptFreqs        => scalar, how many sfs to test
    %  nTrialsPerBlock  => half this number are allocated to each targ location
    %  nBlocks          => number of block per experiment
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function setupGaborColors(LMSfromREX, alphaGuesses, sptFreqfromREX, nSptFreqs, nTrialsPerBlock, nBlocks)
        %Log space spatial frequencies between the upper and lower bound.
        sptFreqfromREX(~sptFreqfromREX) = eps; %no zeros allowed.
        sfs = logspace(log10(sptFreqfromREX(1)), log10(sptFreqfromREX(2)), nSptFreqs);
        degPerCyc = 1./sfs;
        gabor.lambdas = round(ptb.pixperdeg .* degPerCyc); %in pix
        
        %LMSfromREX should be a 9 element vector. Transform this so that
        %there are three rows and triplets go across rows
        exptColor = reshape(LMSfromREX, 3,3)';
        %determine the monitor gamut for each color
        norms = sqrt(sum(exptColor.^2,2));
        units = exptColor ./ repmat(norms, 1, 3);
        for c = 1:size(units, 1);
            if norms(c) == 0; %i.e. this color is disabled
                maxContrast(c) = 0;
                continue %skip to the next color
            end
            
            guess = 0;
            incSize = 0.1; %1 tenth of a %CC
            while(1)
                tmpLMS = units(c,:) .* guess;

                %peak of the gabor:
                peakrgb = ptb.M \((1+tmpLMS./100)' .* ptb.bkgndlms);
                peakRGB = round(ptb.maxdac .* peakrgb) + 1;
                %trough of the gabor:
                troughrgb = ptb.M \((1-tmpLMS./100)' .* ptb.bkgndlms);
                troughRGB = round(ptb.maxdac .* troughrgb) + 1;
                inGamut = all([troughRGB; peakRGB] < (ptb.maxdac+1)) && all([troughRGB; peakRGB] > 0);

                if inGamut
                    guess = guess + incSize;
                else
                    maxContrast(c) = guess - incSize;
                    break
                end
            end
        end

        
        %set up each quest prior distribution based on the maxContrast.
        %There's one for each color/sfs combo
        for c = 1:size(exptColor, 1);
            if maxContrast(c) == 0;
                domain{c} = zeros(1,50);
            else
                domain{c} = 0.01:0.05:maxContrast(c);
            end
            
            for s = 1:length(gabor.lambdas);
                prior = 1./(q.sdPrior.*sqrt(2.*pi)) .* exp(-(domain{c} - alphaGuesses(c)).^2 ./ (2.*q.sdPrior.^2));
                Q{c,s} = log(prior);
            end
        end
        
        % trim down the representation to reflect the non-zero color
        % directions
        validColor = find(sum(abs(exptColor),2) ~= 0);
        q.colorDirs = units(validColor, :); %delete LMS triplets that are all zeros
        q.domain = domain(validColor);
        q.Q = Q(validColor, :);
        
        % Rex (and the data file) expect domains for each color direction,
        % even if they're set to zero. Send rex what it expects
        q.domainForRex = domain;
        
        %set up the trial type matrix ([<color> <sfs> <targ> <trialCounter>])
        gabor.blockCounter = nBlocks;
        gabor.halfTrials = ceil(nTrialsPerBlock./2); %so that trials get allocated to each target evenly
        numColors = size(q.colorDirs,1);
        counters = fullfact([numColors, nSptFreqs, 2]); %2, one for each target
        gabor.tTypes = [counters, ones(size(counters,1), 1).*gabor.halfTrials]; %make a column of counters.
    end

    %
    %   HABITUATION STIMULUS SETUP
    %
    % Gets called by rex once per experiment, and initializes the
    % habituation movie. I make the movie once, and replay it when need be.
    % Relevant inputs:
    %
    % driftRate =>  in cyc/sec
    % sptFreq =>    in cyc/deg
    % orient =>     in radians
    % color =>      a three vector: [LMS] in CC between 0&100
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function setupHabituation(length, driftRate, sptFreq, orient, color)
        habit.length = length;
        pixpercycle = ptb.pixperdeg ./ sptFreq;
        xinc = (2*pi/pixpercycle) * cos((pi/2) - orient);
        yinc = (2*pi/pixpercycle) * sin((pi/2) - orient);
        [xramp, yramp] = meshgrid(xinc*(0:2:ptb.monWidth-1), yinc*(0:ptb.monHeight-1));
        habitlms = ((color./100)+1) .* ptb.bkgndlms';
        habitrgb = ptb.M \ habitlms(:); %same as inv(M)*v but matlab likes this better
        
        %determine the habit drift rate such that the number of frames per
        %cycle is an integer. throw an error if this drift rate is outside
        %of 0.5 Hz of the intended drift rate. I'm jumping through this
        %hoop b/c i just want to caluculate a single cycle worth of movie
        %and then loop it.
        numFrames = ptb.refreshRate / driftRate;
        if driftRate == 0
            habit.driftRate = 0;
            numFrames = 2; %I subtract one from this later!!!
        elseif numFrames == round(numFrames)
            habit.driftRate = driftRate;
        else
            habit.driftRate = ptb.refreshRate / round(numFrames);
            numFrames = round(numFrames);
            if abs(habit.driftRate - driftRate) > 0.5;
                error('habituation drift rate is incorrect')
            end
        end
        
        %make the movie
        phase = 0;
        habit.nFramesPerCycle = numFrames-1;
        deltaPhase = habit.driftRate * 2 * pi / ptb.refreshRate;
        for frame = 1:habit.nFramesPerCycle;
            amp = cos(xramp + yramp + phase);
            im = nan(ptb.monHeight, round(ptb.monWidth./2), 3);
            for gun = 1:3;
                tmp = amp .* (habitrgb(gun) - ptb.bkgndrgb(gun)) + ptb.bkgndrgb(gun);
                tmp = round(tmp .* ptb.maxdac) + 1; %index b/w 1&2^16
                tmp = ptb.invGamma(tmp, gun);
                tmp = round(tmp .* ptb.maxdac); %in DAC units
                im(:,:,gun) = reshape(tmp, ptb.monHeight, round(ptb.monWidth./2));
            end
            habit.movie(frame) = Screen('MakeTexture', ptb.w, TranslateToColourMode(im,1));
            phase = phase+deltaPhase;
            
            %display a waitbar thingy so the user doesn't go bonkers
            %wondering if the code is actually running
            if ~rem(frame, 3)
                angle = 360 .* frame./habit.nFramesPerCycle;
                arcRect = pos2rect([-88, 68], 4);
                Screen('FillArc', ptb.w, [], arcRect, 0, angle);
                Screen('Flip', ptb.w);
            end
        end
        
        %inform REX that the habituation setup is complete
        pnet(udp.sock, 'write', 'MAChabitSetupComplete>> >>');
        pnet(udp.sock, 'writepacket', udp.rexip, udp.port);
    end
    
    %
    %   TOPUP HABITUATION
    %
    % Gets called by Rex each time a habituator needs to be presented. This
    % basically just specifies how long the habit will be on the screen.
    % The center of the habit stim is painted gray. The size of this area
    % can change from trial to trial. eyewinSize in .1DVA and habitDuration
    % is in milliseconds.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function topupHabit(eyewinSize)
        %a hack for the massive habituation at the start of each run.
        if trial.firstTrial && habit.length
            habitDuration = 30000;
            habitDuration = habit.length;
        else
            habitDuration = habit.length;
        end
        habit.centGrayRect = pos2rect(ptb.monCent, 2*eyewinSize); %paint this area gray in drawTrialElement
        habit.nFramesToPresent = ptb.refreshRate * (habitDuration ./ 1000);
        habit.nFramesToPresent = ceil(habit.nFramesToPresent); %can't have a non integer!!
        habit.frameCounter = 0; %zero this here and not in allOff!!
    end


    %
    %   MONITOR SETUP
    %
    % Gets called by rex once per experiment. mondist and screenwidth are
    % specified in a rig-specific header file on each rex machine. Both
    % input arguments are in units of cm.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function setupMonitor(mondist, screenwidth, calPath)
        
        % compute relevant calibration data
        S = load('T_cones_smj.mat');
        load(calPath)
        cal = cals{end};
        fundamentals = S.T_cones_smj;
        ptb.gammaTable = cal.gammaTable;
        ptb.bkgndRGB = round(255 .* cal.bgColor); %bkgnd voltages discretized b/w 1&255
        ptb.bkgndrgb = [ptb.gammaTable(ptb.bkgndRGB(1)+1,1), ptb.gammaTable(ptb.bkgndRGB(2)+1,2), ptb.gammaTable(ptb.bkgndRGB(3)+1,3)];
        ptb.monSpd = SplineSpd(cal.S_device, cal.P_device, S.S_cones_smj);
        ptb.M = fundamentals * ptb.monSpd;
        ptb.invGamma = InvertGamma(cal, 1);
        ptb.bkgndlms = ptb.M * ptb.bkgndrgb';
        targ.RGB = ptb.bkgndRGB .* 0.70; %targs are 70% of bkgnd.
        
        %itialize PTB and compute all the relevant monitor dimensions
        if isempty(Screen('windows'));
            ptb.w = Screen('OpenWindow', 0, ptb.bkgndRGB);
        else
            winds = Screen('windows');
            if numel(winds) > 1
                for w = 2:numel(winds); %necessary in R2009b?
                    Screen('Close', winds(w)); %b/c the habit windows don't close themselves.
                end
            end
            ptb.w = Screen('windows');
            Screen('FillRect', ptb.w, ptb.bkgndRGB)
        end
        HideCursor;
        [ptb.monWidth, ptb.monHeight] = Screen('WindowSize', ptb.w);
        theta = atand((screenwidth./2) ./ mondist);
        ptb.pixperdeg = (ptb.monWidth ./ 2) ./ theta;
        ptb.oldGamma = Screen('LoadNormalizedGammaTable', ptb.w, repmat(linspace(0,1,256)', 1, 3));
        ptb.refreshRate = Screen('NominalFrameRate', ptb.w, 1);
        ptb.monCent = [ptb.monWidth, ptb.monHeight] ./ 2; %(x,y) at the center of the screen
        
                
        %debuging code:
        x = 0:255; %the normal domain of the gamma look up table
        xx = linspace(0, 255, 2^16); %the desired quantization of the gammaTable
        g1 = reshape(ptb.gammaTable, 256, 3);
        ptb.bigGamma = [spline(x, g1(:,1), xx)', spline(x, g1(:,2), xx)', spline(x, g1(:,3), xx)'];
        
        
        %inform REX that the monitor setup is complete
        pnet(udp.sock, 'write', 'MACmonitorSetupComplete>> >>');
        pnet(udp.sock, 'writepacket', udp.rexip, udp.port);
    end

    %
    %   SETUPTRIAL
    %
    %  lastTrialGood    => binary flag
    %  lastTrialCorrect => binary flag used in staircasing
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function setupTrial(lastTrialGood, lastTrialCorrect)
              
        %pick a color, contrast, and spatial frequency to present.
        [gabor.trgb, gabor.tLambda, gabor.tColor, gabor.tSfs, gabor.tLoc] = Quest(lastTrialGood, lastTrialCorrect);
        if any(isnan(gabor.trgb))
            disp('End of Experiment');
            return
        end
        
        %assign the stimulus location
        if gabor.tLoc == 1
                gabor.tPos = gabor.T1pos;
        elseif gabor.tLoc == 2
                gabor.tPos = -gabor.T1pos;
        end
        
        %recompute the gabor rectangle based on it's current location
        sizeInTenthsDeg = gabor.sigma .* gabor.nSigmas .* 2; %multiply by two so that +/-nSigmas mass is shown
        gabor.rect = pos2rect(gabor.tPos, sizeInTenthsDeg);
        
        %zero out the frame counter here (or else it won't get dropped in
        %the data file) and prepare the stimulus movie
        gabor.frameCounter = 0;
        gabor.movie = makeGaborMovie();
        
        %tell rex that the current trial is ready to go
        pnet(udp.sock, 'write', 'MACtrialSetupComplete>> >>');
        pnet(udp.sock, 'writepacket', udp.rexip, udp.port);
    end

    %
    %    QUEST TRIAL GENERATOR
    %
    % This function updates quest functions (stored seperately for each
    % color/sfs combo) on each trial. It also selects the next trial's
    % contrast according to the modal value of the quest function. 
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [rgbTrial, lambda, color, sfs, targ] = Quest(lastTrialGood, lastTrialCorrect)
        trial.lastTrialCorrect = lastTrialCorrect; %for debugging
        if lastTrialGood
            %update the counter from the previous trial
            if trial.firstTrial; trial.firstTrial = 0; end
            [dum, row] = ismember([gabor.tColor, gabor.tSfs, gabor.tLoc], gabor.tTypes(:,1:3), 'rows');
            gabor.tTypes(row, end) = gabor.tTypes(row, end) -1;
            
            %update the quest function from the previous trial
            tCC = norm(((ptb.M*gabor.trgb(:))-ptb.bkgndlms) ./ ptb.bkgndlms);
            tCC = tCC.*100; %convert to %CC 1to100;
            liklihood = 1-0.5.*exp(-(tCC./q.domain{gabor.tColor}).^q.modBeta);
            if ~lastTrialCorrect
                liklihood = 1-liklihood;
            end
            q.Q{gabor.tColor, gabor.tSfs} = q.Q{gabor.tColor, gabor.tSfs} + log(liklihood);
        end
        
        %obtain stimulus parameters for the upcomming trial
        [color, sfs, targ] = trialAllocator();
        if ~any([color, sfs, targ])
            [rgbTrial, lambda] = deal(NaN); % crappy way to deal with the end of an expt
            return
        end
        
        
        [dum, idx] = max(q.Q{color, sfs});
        coneContrast = q.domain{color}(idx);
        LMS = (q.colorDirs(color,:) .* coneContrast)./100; % in %CC (0to1)
        rgbTrial = [ptb.M \ ((1+LMS)' .* ptb.bkgndlms)]'; %make sure it comes out as a row vec
        lambda = gabor.lambdas(sfs);
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
    function [color, sfs, targ] = trialAllocator()
        %figure out which trial types are fair game. Don't advance to the
        %next block of trials until all trial types have been exhausted for
        %the current block
        typesAvailable = find(gabor.tTypes(:,end) ~= 0);
        %if no types are available than advance to the next block. If the
        %last block just finished than send the "all done" signal to REX.
        if isempty(typesAvailable)
            if (gabor.blockCounter == 1) %all done?
                pnet(udp.sock, 'write', 'MACexptDone>> >>');
                pnet(udp.sock, 'writepacket', udp.rexip, udp.port);
                trial.exptDone = 1;
                [color, sfs, targ] = deal(NaN); %return something so that the program doesn't bonk
                return
            else % decrement the block counter, reinitialize the tType counters
                gabor.blockCounter = gabor.blockCounter - 1;
                gabor.tTypes(:,end) = gabor.halfTrials;
                [color, sfs, targ] = trialAllocator(); %call recursivley for the new tType
                return
            end
        end

        %at long last, find a trial type for the next trial
        randomDraw = unidrnd(length(typesAvailable));
        tTypeInd = typesAvailable(randomDraw);
        color = gabor.tTypes(tTypeInd, 1);
        sfs = gabor.tTypes(tTypeInd, 2);
        targ = gabor.tTypes(tTypeInd, 3);
    end

    %
    %          MAKE THE GABOR MOVIE
    %
    % This function makes all frames of the stimulus so that the timing of
    % gabor presentation can be tightly controlled. It also keeps track of the
    % exptected and actual RGBs of the gabor 
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function movie = makeGaborMovie()
        %determine the gabor's RGB for the current trial
        RGB = round(ptb.maxdac .* gabor.trgb) + 1;
        gabor.tRGB = zeros(1,6);
        gabor.tRGB(1:3) = round(ptb.maxdac .* [ptb.invGamma(RGB(1), 1), ptb.invGamma(RGB(2), 2), ptb.invGamma(RGB(3), 3)]);
        
        %make the movie
        phase = 0;
        gaborMaxInc = 0;
        deltaPhase = gabor.speed .* 2 .* pi .* (1 ./ ptb.refreshRate); %the amount to advance each frame (in radians)
        halfSize = (gabor.rect(3) - gabor.rect(1)) ./ 2;  %base the size of the gabor off of the flash's texture window so that there are no mistakes
        row = -halfSize:2:halfSize-1; %subtract one so that you don't overlap the texture window;
        col = -halfSize:halfSize-1;
        [X, Y] = meshgrid(row, col);
        xprime = X .* cos(-gabor.theta) + Y .* sin(-gabor.theta);
        yprime = -X .* sin(-gabor.theta) + Y .* cos(-gabor.theta);
        rgb_increment = gabor.trgb - ptb.bkgndrgb;
        sigmaInPix = gabor.sigma ./ 10 .* ptb.pixperdeg;
        for a = 1:gabor.numFrames;
            %start by making the gabor, then multiplying by the increment from
            %background and by the temporal weighting fxn.
            template = exp(-(xprime.^2 + gabor.gamma.^2 .* yprime.^2) ./ (2.* sigmaInPix.^2)) .* cos(2 .* pi .* yprime ./ gabor.tLambda + phase);
            flashImg = repmat(template, [1, 1, 3]);
            flashImg(:,:,1) = (flashImg(:,:,1) .* gabor.timeProf(a) .* rgb_increment(1)) + ptb.bkgndrgb(1);
            flashImg(:,:,2) = (flashImg(:,:,2) .* gabor.timeProf(a) .* rgb_increment(2)) + ptb.bkgndrgb(2);
            flashImg(:,:,3) = (flashImg(:,:,3) .* gabor.timeProf(a) .* rgb_increment(3)) + ptb.bkgndrgb(3);
            

            %now convert to DAC values
            imgSize = size(template);
            flashrgb = round(ptb.maxdac .* flashImg) + 1; %intensities b/w 1 & maxdac +1
            flashRGB = ones(size(flashImg)); %preallocate so that indexing works
            flashRGB(:,:,1) = reshape(round(ptb.maxdac .* ptb.invGamma(flashrgb(:,:,1), 1)), imgSize(1), imgSize(2));
            flashRGB(:,:,2) = reshape(round(ptb.maxdac .* ptb.invGamma(flashrgb(:,:,2), 2)), imgSize(1), imgSize(2));
            flashRGB(:,:,3) = reshape(round(ptb.maxdac .* ptb.invGamma(flashrgb(:,:,3), 3)), imgSize(1), imgSize(2));
  
            %drop the RGB from the highest contrast pix into the data file
            bkgndImg = repmat(permute(ptb.bkgndrgb, [1,3,2]), size(flashImg,1), size(flashImg,2));
            diffFromBkgnd = flashImg - bkgndImg;
            norms = sqrt(sum(diffFromBkgnd.^2, 3)); %in normalized gun space
            frameMaxInc = max(norms(:));
            if frameMaxInc > gaborMaxInc;
                [r, c] = find(norms == frameMaxInc, 1);
                gabor.tRGB(1,4:6) = permute(flashRGB(r,c,:), [1,3,2]);  % DAC values b/w 0&maxdac
                gaborMaxInc = frameMaxInc;
            end

            %make a movie frame and update the spatial phase
            movie(a) = Screen('MakeTexture', ptb.w, TranslateToColourModeMex(flashRGB));
            phase = phase + deltaPhase;
        end
        

        actrgb = [ptb.bigGamma(gabor.tRGB(4)+1, 1), ptb.bigGamma(gabor.tRGB(5)+1, 2), ptb.bigGamma(gabor.tRGB(6)+1, 3)];
        actlms = ptb.M * actrgb(:);
        actLMS = (actlms - ptb.bkgndlms) ./ ptb.bkgndlms;
        reqlms = ptb.M * gabor.trgb(:);
        reqLMS = (reqlms - ptb.bkgndlms) ./ ptb.bkgndlms;
        prcntDiff = abs(norm(actLMS)-norm(reqLMS)) ./ norm(reqLMS);
        if prcntDiff > 0.1
                sca
                disp([[norm(actLMS), norm(reqLMS)].*100, trial.lastTrialCorrect]);
                error('actual and requested stimuli differ by more than 10 percent')
        end
    end

    %
    %  SET EXPT PARAMETERS
    %
    % Sets the local experimental parameters based on those specified by
    % REX. This function gets called before setupColorContrasts.
    %
    %  fpPos            => a 2 vector (x,y) in tenths of degrees
    %  fpSize           => in tenths of degrees
    %  targSize         => tenths of degrees
    %  T1 & T2Dist      => from fp in tenths of degrees
    %  RFpos            => a 2 vector, the (x,y) of the RF in tenths of degrees
    %  flashLength      => in milliseconds
    %  nSigmas          => number of SDs to present (effectively the gaborsize)
    %  gaborTheta       => in radians
    %  gaborSigma       => one sided standard dev in tenths of degrees
    %  gaborGamma       => a floating point number
    %  gaborSpeed       => temporal freq in cyc/sec
    %  priorSD          => The SD of the prior Quest distribution in CC
    %  beta             => slope of the models psyFun
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function setExptParams(fpPos, fpSize, targSize, T1dist, T2dist, RFpos, flashLength, nSigmas, gaborTheta, gaborSigma, gaborGamma, gaborSpeed, priorSD, beta)
        % fp params
        fp.rect = pos2rect(fpPos, fpSize);
        
        %targ params
        targ.T1rect = pos2rect(round(RFpos .* T1dist), targSize);
        targ.T2rect = pos2rect(round(-RFpos .* T2dist), targSize);
        
        % quest params
        q.sdPrior = priorSD;
        q.modBeta = beta;
        
        %gabor params
        gabor.T1pos = RFpos;
        gabor.length = flashLength;
        gabor.nSigmas = nSigmas;
        gabor.theta = gaborTheta;
        gabor.sigma = gaborSigma;
        gabor.gamma = gaborGamma;
        gabor.speed = gaborSpeed;
        gabor.numFrames = ceil(ptb.refreshRate .* (gabor.length./1000));

        %deal with the time weighting function here (b/c it doesn't change
        %between trials)
        rampLength = ceil(gabor.numFrames ./ 4);
        ramp = linspace(0, 1, rampLength);  %ramp is 1/4th of the total duration on either side
        plateau = ones(1,gabor.numFrames - (rampLength .* 2));
        gabor.timeProf = [ramp, plateau, fliplr(ramp)];
    end

%
%   SEND QUEST COLOR RANGES TO REX
%
%  Instead of sending the entire string of numbers I'll just send the
%  first, last, and stepsize. I'll send over the data for each color
%  direction tested as a vector of numbers. Each successive triplet is for
%  a new color. Since Rex is expecting a 9 vector i'll pad with zeros for
%  the cases where I use fewer than 3 colors.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function sendQuestColorRanges()
        msg = zeros(1,9);
        for a = 1:length(q.domainForRex);
            msg((a.*3)-2) = min(q.domainForRex{a});
            msg((a.*3)-1) = max(q.domainForRex{a});
            msg(a.*3) = q.domainForRex{a}(2) - q.domainForRex{a}(1); %pretty kludgy
        end
        
        %now send the message using sendToRex
        sendToRex(udp, msg, 'double', 'sendQuestColorRanges();');
        
    end
    
    %
    %          MAKE A DRAWING RECTANGLE
    %
    % draws a square of side length "size" (in .1DVA) arround a point on the
    % monitor specified by centPos ([x,y] in .1DVA). rect is the rectangle
    % in units of pixels
    %
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
    %   UDP COMMUNICATION
    %
    % Reads UDP string messages from Rex and 'evals' them as if these
    % statements were evaluated in the command window.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function dealWithUdpMessages()

        msgSize = pnet(udp.sock, 'readpacket', 260, 'noblock');
        if ~msgSize
            return
        end
        
        
        %open the message
        message = pnet(udp.sock, 'read', msgSize, 'char');
        %disp(message); %for debugging
        %see if we should return to MasterSlave
        if strcmpi('return', message(1:6));
            ptb.exitNow = 1;
        end
        
        try
            %try to execute the message
            eval(message);
            
        catch
            fprintf('Unknown message: %s\n', message);
            sca
            wtf %print out the error message
        end
    end

end %HABITslave



%%%%%%%%%%  QUESTIONS %%%%%%%%%%%%


% Can I make the while loop (while ~ptb.exitNow) ??

% does the floor/ceil thing in pos2rect matter??


%% TESTING PROCEEDURES

% I added a line to habit.m that dictated lastTrialCorrect when the
% stimulus was above 3% CC (line 513 in Quest subfun). This had the effect of simulating an observer
% with a hard threshold equal to 3%. In the offline analysis (habitUnpack)
% the modal values are messed up (which is expected b/c the
% lastTrialCorrect disagree between the slave and spot files). More
% importantly, the asympotic CC presented by the slave is exactly 3% CC.
% This demonstrates that the QUEST algorithm is working properly.