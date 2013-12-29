function DTslaveTraining

    % This slave program should accept commands from rex and present a flash of
    % color in a particular location in space for a particular amount of
    % time.

    %***************  SETUP GLOBALS AND CONSTANTS *************%
    global udpCom;
    udpCom.sock = [];
    udpCom.port = 6665;
    udpCom.rexip = '192.168.1.120';
    udpCom.plexip = '192.168.1.121';

    KbName('UnifyKeyNames');
    key_esc = KbName('escape');
    
    global ptb;
    ptb.w = [];
    ptb.monSize = [];       %(x, y, x, y) of TL and BR corners
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
    ptb.ccmode = 1;
    ptb.vpixx = ConnectedToViewPixx();
    ptb.exitNow = 0;        %a flag for integration with masterslave

    global stim;
    stim.flashFrameCounter = [];    %number of frames that flash has been on screen.
    stim.flashTimeProf = [];        %temporal profile of the stimulus
    stim.flashNumFrames = [];       %number of frames the flash is on the monitor
    stim.flashRect = [];            %the texture drawing rectangle
    stim.flashDeltaPhase = [];
    stim.gaborTheta = [];           %a float: the gabor's orientation in radians
    stim.gaborLambda = [];          %an int; wavelength in pixels
    stim.gaborPhase = [];           % a float; in radians
    stim.gaborSigma = [];           %a float; standard dev of gabor (comes in degs, but converted to pix)
    stim.gaborGamma = [];           %a float; 1 give circles, else elipses
    stim.gaborSpeed = [];           %a float; cycles per second
    stim.fpRGB = [];
    stim.fpRect = [];               %drawing rectangle for the FP
    stim.frameRectIpsi = [];        %drawing Rect for wire frame surrounding the flash
    stim.frameRectContra = [];
    stim.frameMode = 2;             %1=draw around flash only, 2=draw both
    stim.frameRGB = [0 0 0];        %frame is black by default
    stim.targCorrectRect = [];
    stim.targIncorrectRect = [];
    stim.targPosInDeg = [];
    stim.targMode = [];
    stim.targ_bad_rgb = [0 0 0];

    global trial;
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
                Screen('FillRect', ptb.w, stim.fpRGB, stim.fpRect);
            end

            if(trial.drawFrame)
                Screen('FrameRect', ptb.w, stim.frameRGB, stim.frameRectIpsi, 1);
                if (stim.frameMode == 2);
                    Screen('FrameRect', ptb.w, stim.frameRGB, stim.frameRectContra, 1);
                end
            end

            if(trial.drawTarg)
                Screen('FillRect', ptb.w, [0,0,0], stim.targCorrectRect);
                if (stim.targMode == 2);
                    Screen('FillRect', ptb.w, stim.targ_bad_rgb, stim.targIncorrectRect);
                end
            end

            if(trial.drawFlash)
                drawFlash();

                if(stim.flashFrameCounter == stim.flashNumFrames);
                    pnet(udpCom.sock, 'write', 'MACstimLastFrame>> >>');
                    pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
                    trial.drawFlash = 0;
                end

                if(stim.flashFrameCounter == 1)
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
        if (keyisdown && keycode(key_esc)) %if you press escape
            if(ptb.w > 0)
                pnet(udpCom.sock, 'close');
                sca();
            end
            return
        end
    end %main while loop
end %DTslave

%%
function fpToggle(onOff, fpRGB)
    global trial;
    global stim;

    trial.drawFP = onOff;
    stim.fpRGB = fpRGB/255;
end

%%
function wireFrameToggle(onOff, frameMode, frameOn)
    global trial;
    global stim;
    global ptb;

    trial.drawFrame = onOff;
    stim.frameMode = frameMode;

    if frameOn
        stim.frameRGB = [0 0 0];
    else
        stim.frameRGB = ptb.bkgndRGB/255;
    end
end

%%
function flashToggle(onOff)
    global trial;

    trial.drawFlash = onOff;

end

%%
function targToggle(onOff, targMode, bad_targ_rgb)
    global stim;
    global trial;
    global ptb;

    trial.drawTarg = onOff;
    trial.drawFP = 0;   %turn off the FP at the same instant you illuminate the targets.
    stim.targMode = targMode;
    bad_rgbs = round(255 * ((bad_targ_rgb./100) .* ptb.bkgndrgb)) + 1;
    stim.targ_bad_rgb = [ptb.smallInvGamma(bad_rgbs(1), 1), ptb.smallInvGamma(bad_rgbs(2), 2), ptb.smallInvGamma(bad_rgbs(3), 3)];
end

%%
function drawFlash()
    global ptb;
    global stim;

    %increment the counter
    stim.flashFrameCounter = stim.flashFrameCounter+1;
    Screen('Close');

    %start by making the gabor, then multiplying by the increment from
    %background and by the temporal weighting fxn.
    gabor = makeGabor(stim.gaborTheta, stim.gaborLambda, stim.gaborPhase, stim.gaborSigma, stim.gaborGamma);
    rgb_increment = stim.rgbTrial - ptb.bkgndrgb;
    flashImg = repmat(gabor, [1, 1, 3]);
    flashImg(:,:,1) = (flashImg(:,:,1) .* stim.flashTimeProf(stim.flashFrameCounter) .* rgb_increment(1)) + ptb.bkgndrgb(1);
    flashImg(:,:,2) = (flashImg(:,:,2) .* stim.flashTimeProf(stim.flashFrameCounter) .* rgb_increment(2)) + ptb.bkgndrgb(2);
    flashImg(:,:,3) = (flashImg(:,:,3) .* stim.flashTimeProf(stim.flashFrameCounter) .* rgb_increment(3)) + ptb.bkgndrgb(3);

    %for debugging descrepancies b/w reported and actual cone contrasts:
    stim.flashRGB(1:3) = rgb_increment;
    frameMaxInc = max(max(abs(flashImg(:,:,3) - ptb.bkgndrgb(3))));
    if frameMaxInc > stim.flashMaxInc;
        [r, c] = find(abs(flashImg(:,:,3) - ptb.bkgndrgb(3)) == frameMaxInc);
        r = r(1);
        c = c(1);
        stim.flashRGB(4) = flashImg(r,c,1) - ptb.bkgndrgb(1);
        stim.flashRGB(5) = flashImg(r,c,2) - ptb.bkgndrgb(2);
        stim.flashRGB(6) = flashImg(r,c,3) - ptb.bkgndrgb(3);
        stim.flashMaxInc = frameMaxInc;
    end


    % this is a hacknied way dealing with LMS triplets that are out of the
    % monitor gammut, but will work for now.
    try
        %now convert to DAC values
        imgSize = size(gabor);
        flashrgb = round(ptb.maxdac .* flashImg) + 1; %intensities b/w 1 & maxdac
        flashRGB = ones(size(flashImg)); %preallocate so that indexing works


        flashRGB(:,:,1) = reshape(ptb.invGamma(flashrgb(:,:,1), 1), imgSize(1), imgSize(2));
        flashRGB(:,:,2) = reshape(ptb.invGamma(flashrgb(:,:,2), 2), imgSize(1), imgSize(2));
        flashRGB(:,:,3) = reshape(ptb.invGamma(flashrgb(:,:,3), 3), imgSize(1), imgSize(2));
    catch
        %fprintf('Out of bounds error, could not display stimulus\n');
    end

    flashTexture = Screen('MakeTexture', ptb.w, TranslateToColourMode(flashRGB), [], [], 2);
    Screen('DrawTexture', ptb.w, flashTexture, [], stim.flashRect, [], 0);

    %update the position of the gabor for the next frame refresh
    stim.gaborPhase = stim.gaborPhase + stim.flashDeltaPhase;
end

%%
function allOff()
    global trial;

    trial.drawFP = 0;
    trial.drawFrame = 0;
    trial.drawFlash = 0;
    trial.drawTarg = 0;

end
%%
function rect = pos2rect(centPos, size)
    %draws a square of side length "size" (in pix) arround a point on the
    %monitor specified by centPos (in (x,y) cordinates) size is in
    %units of pixels

    global ptb;

    %centPos comes from rex in units of tenths of degrees so convert to
    %pixels first
    xyPosX = round(centPos(1) ./ 10 .* ptb.pixperdeg);
    xyPosY = round(centPos(2) ./ 10 .* ptb.pixperdeg);

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

%%
function setupTrial(flashLength, fpPos, fpSize, frameOffset, targSize, targDist, flashOn, gaborTheta, gaborLambda, gaborGamma, gaborSigma, gaborSpeed, flashSize, flashPos, lastTrialGood)
    %accepts parameters from rex then sets up the stimuli accordingly.
    %flashSize, fpSize and frameOffset should be in tenths of degrees.
    %stim.flashPos and fpPos should be with respect to the center of the
    %screen in tenths of degrees. flashLength should be in units of
    %milli-seconds. targDist should be b/w 0&1 signifying the distance b/w
    %the fp and the flash represented in thousandths.
    global stim;
    global ptb;

    % reset the frame number counter and delete all the textures that are
    % currently in the frame buffer. Reset the phase of the gabor
    stim.flashFrameCounter = 0;
    stim.gaborPhase = 0;
    Screen('Close') %I think that this closes/deletes all offscreen textures. redundant with above.


    %for debugging differences b/w expected and presented cone contrasts
    stim.flashMaxInc = 0;

    %set up the spatial location of the FP, wire frame, stimulus flash, and
    %saccade targets
    fpSize = round(fpSize ./ 10 .* ptb.pixperdeg);
    stim.fpRect = pos2rect(fpPos, fpSize);

    %deal with the flash
    flashSize = round(flashSize ./ 10 .* ptb.pixperdeg);
    stim.flashRect = pos2rect(flashPos, flashSize);

    frameOffset = round(frameOffset ./ 10 .* ptb.pixperdeg);
    stim.frameRectIpsi = pos2rect(flashPos, flashSize+frameOffset);
    stim.frameRectContra = pos2rect((flashPos.*-1), flashSize+frameOffset);

    %determine the location of the targets. Remeber that two right
    %triangles that have the same angles but different side lengths have
    %the same ratio of side lengths.
    targSize = round(targSize ./ 10 .* ptb.pixperdeg);
    stim.targPosInDeg = round(flashPos .* targDist);
    stim.targCorrectRect = pos2rect(stim.targPosInDeg, targSize);
    stim.targIncorrectRect = pos2rect((stim.targPosInDeg.*-1), targSize);

    %deal with the time weighting function
    stim.flashNumFrames = ceil(ptb.refreshRate .* (flashLength./1000));
    rampLength = ceil(stim.flashNumFrames ./ 6);
    ramp = linspace(0, 1, rampLength);  %ramp is 1/6th of the total duration on either side
    plateau = ones(1,stim.flashNumFrames - (rampLength .* 2));
    stim.flashTimeProf = [ramp, plateau, fliplr(ramp)];

    %set the flash to bkgnd for training purposes. Setting the time profile
    %to zero forces the flash to be at bkgnd b/c we'll add 0 .* the spatial
    %weighting fxn for each frame (i.e. never add anything to bkgnd)
    if ~flashOn
        stim.flashTimeProf = stim.flashTimeProf .* 0;
    end

    %setup the gabor
    stim.gaborTheta = gaborTheta;
    stim.gaborLambda = round(ptb.pixperdeg .* gaborLambda ./ 10);
    stim.gaborSigma = round(ptb.pixperdeg .* gaborSigma ./ 10);
    stim.gaborGamma = gaborGamma;
    stim.gaborSpeed = gaborSpeed;
    stim.flashDeltaPhase = gaborSpeed .* 2 .* pi .* (1 ./ ptb.refreshRate); %the amount to advance each frame

    %decrement the contrast counter from the last good trial
    if lastTrialGood
        stim.rgbCounters(stim.contrastType) = stim.rgbCounters(stim.contrastType) - 1;
    end

    %now pick the color of the stimulus at random from the available
    %possibilities, but only if the monkey got the last trial correct
    contrastsLeft = find(stim.rgbCounters ~= 0);
    if ~isempty(contrastsLeft)
        randomDraw = unidrnd(length(contrastsLeft));
        stim.contrastType = contrastsLeft(randomDraw);
        stim.rgbTrial = stim.rgbLevels(stim.contrastType, :);
        invGammaInds = round(ptb.maxdac .* stim.rgbTrial) + 1;
        %         stim.flashRGB = round(ptb.maxdac .* [ptb.invGamma(invGammaInds(1), 1), ptb.invGamma(invGammaInds(2), 2), ...
        %                                                ptb.invGamma(invGammaInds(3), 3)]);
    else
        pnet(udpCom.sock, 'write', 'MACexptdone>> >>');
        pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
    end
end


%%
function setupColorContrasts(LMSfromREX, nContrasts, nTrialsPerContrast)
    global stim ptb;

    %set up the trial types:
    if nContrasts > 1

        coneLevels(:, 2:nContrasts+1) = [logspace(log10(2), log10(LMSfromREX(1)), nContrasts);...
            logspace(log10(2), log10(LMSfromREX(2)), nContrasts);...
            logspace(log10(2), log10(LMSfromREX(3)), nContrasts)];
        nContrasts = nContrasts+1;
    else
        coneLevels = [logspace(0, log10(LMSfromREX(1)), nContrasts);...
            logspace(0, log10(LMSfromREX(2)), nContrasts);...
            logspace(0, log10(LMSfromREX(3)), nContrasts)];
    end

    coneLevels = ((coneLevels ./ 100) + 1) .* repmat(ptb.bkgndLMS, 1, nContrasts); % recall that the flashLMS get sent accross as percent

    stim.rgbLevels = [inv(ptb.M) * coneLevels]';  %transform to gun space

    stim.rgbCounters = ones(1,nContrasts) .* nTrialsPerContrast; %count down from the max num trials
end

%%
function setupMonitor(mondist, screenwidth, calFilePath)
    global ptb udpCom

    %load in the calibration data and fill up the ptb structure
    load('T_cones_smj.mat');
    load(calFilePath);
    cal = cals{end};
    fundamentals = T_cones_smj;
    ptb.bkgndRGB = round(255 * cal.bgColor);
    ptb.bkgndrgb = [cal.gammaTable(ptb.bkgndRGB(1)+1, 1), cal.gammaTable(ptb.bkgndRGB(2)+1, 2), ...
        cal.gammaTable(ptb.bkgndRGB(3)+1, 3)];
    
    ptb.gammaTable = cal.gammaTable;
    ptb.monSpd = SplineSpd(cal.S_device, cal.P_device, S_cones_smj);
    ptb.M = fundamentals * ptb.monSpd;
    ptb.invGamma = InvertGammaTable(cal.gammaInput, cal.gammaTable, ptb.maxdac+1);
    ptb.smallInvGamma = InvertGammaTable(cal.gammaInput, cal.gammaTable, 2^8);
    ptb.bkgndLMS = ptb.M * ptb.bkgndrgb';
    
    % start up the imaging pipeline
    if ~isempty(Screen('Windows'))
        ptb.w = max(Screen('Windows'));
        Screen('FillRect', ptb.w, ptb.bkgndRGB/255);
    else
        PsychImaging('PrepareConfiguration');
        if ptb.vpixx
            PsychImaging('AddTask', 'General', 'EnableDataPixxC48Output', ptb.ccmode); % in mode '1': every 2nd column of pixels is ignored
        else
            PsychImaging('AddTask', 'General', 'EnableBits++Color++Output', ptb.ccmode);
        end
        ptb.w = PsychImaging('OpenWindow', 0, ptb.bkgndRGB/255);
    end
    HideCursor();
    
    [ptb.monWidth ptb.monHeight] = Screen('WindowSize', ptb.w);
    theta = atand(screenwidth / 2 / mondist);
    ptb.pixperdeg = ptb.monWidth / 2 / theta;
    ptb.oldGamma = Screen('LoadNormalizedGammaTable', ptb.w, repmat(linspace(0,1,256)', 1, 3));
    ptb.refreshRate = Screen('NominalFrameRate', ptb.w, 1);
    ptb.monCent = [ptb.monWidth ptb.monHeight] / 2; %(x,y) at the center of the screen
    
    %inform REX that the monitor setup is complete
    pnet(udpCom.sock, 'write', 'MonitorSetupComplete>> >>');
    pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
end

%%
function gabor = makeGabor(theta, lambda, phi, sigma, gamma)
    global stim;

    halfSize = (stim.flashRect(3) - stim.flashRect(1)) ./ 2;  %base the size of the gabor off of the flash's texture window so that there are no mistakes
    interval = -halfSize:halfSize-1; %subtract one so that you don't overlap the texture window;
    [X, Y] = meshgrid(interval, interval);
    xprime = X .* cos(theta) + Y .* sin(theta);
    yprime = -X .* sin(theta) + Y .* cos(theta);
    gabor = exp(-(xprime.^2 + gamma.^2 .* yprime.^2) ./ (2.* sigma.^2)) .* cos(2 .* pi .* xprime ./ lambda + phi);

    %now average neighboring columns to deal with Bits++ in colour mode
    gabor = (gabor(:, 2:2:end) + gabor(:, 1:2:end)) ./ 2;
end

%%
function dealWithUdpMessages()
    global udpCom;
    global stim;
    global ptb;

    msgSize = pnet(udpCom.sock, 'readpacket', 200, 'no block');
    if ~msgSize
        return
    end

    try
        message = pnet(udpCom.sock, 'read', msgSize, 'char');

        %see if we should return to MasterSlave
        if strcmpi('return', message(1:6));
            ptb.exitNow = 1;
        end

        %disp(message);
        eval(message);
    catch ME
        fprintf('Trouble with message: "%s"\n', message);
        disp(getReport(ME));
    end
end

%%
%******* NOTES AND QUESTIONS **************
% 1) previous editions had a "blank texture" do I need one here?

%3) optional arg to framerate?

%9) fix the RGB rgb and gamma correction.

%11) check method of creating targRects (plus they seem to be offset)

%12) ptb.smallInvGamma is silly, remove from all fxns and figure out how to
    %modulate targ_bad intensity in a better fashion


