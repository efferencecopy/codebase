%% trying to make a movie of a trial of the DT paradigm

%% open a sample file
stro = nex2stro(findfile('S061009002'));

%% Define the trial parameters
DTindicies
trialNum = 224;  %high contrast = 325, low contrast = 224;, out of RF = 98
eyeOffsetX = 0.5;  %trial and error.... just using what looks best (in dva)
eyeOffsetY = 0.5;
bkgndrgb  = [stro.sum.exptParams.bkgnd_r, stro.sum.exptParams.bkgnd_g, stro.sum.exptParams.bkgnd_b];
bkgndrgb = [.5 .5 .5];

pixperdeg = stro.sum.exptParams.pixperdeg;
flashX = stro.trial(trialNum, flashXInd);       % in tenths of degrees
flashY = stro.trial(trialNum, flashYInd);       % in tenths of degrees
nSigmas = stro.sum.exptParams.flash_size;       % the number of sigmas in the flash rect
lambda = 50;                                    % spatial period in pix
sigma = stro.trial(trialNum, gaborSigmaInd);    % sd in tenths of degrees
theta = stro.trial(trialNum, gaborThetaInd);    % in radians
gamma = 1;                                      % sdX == sdY
driftRate = 3;                                  % in cyc/sec
devFrameRate = 30;                              % in frames per second for the resulting movie
oldFrameRate = stro.sum.exptParams.frame_rate;  % in frames per second for the dell 4 computer
epsize = 11;                                    % EP cursor size in pix
epHalfSize = floor(epsize/2);
fpHalfSize = 7;                                 % FP half size in pix
cometWidth = 9;                                 % comet width in pix
cometHalfWidth = floor(cometWidth/2);


%there are 5 distinct trial epochs. Start by determining how much time was
%spent in each of these epochs
preFlashFPTime = stro.trial(trialNum, flashOnInd) - stro.trial(trialNum, fpOnInd);
flashOnTime = stro.trial(trialNum, nFramesInd) .* (1./oldFrameRate);
postFlashDelay = (stro.trial(trialNum, targOnInd) - stro.trial(trialNum, flashOnInd)) - flashOnTime;
targsAndFPTime = 0.010; %just the way it is...
targsOnTime = stro.trial(trialNum, rewOnInd) - stro.trial(trialNum, targOnInd) - targsAndFPTime;


%now determine how many frames should be allocated to each epoch.
preFPFrames = 15;
preFlashFP_frames = ceil(preFlashFPTime .* devFrameRate);
flashOn_frames = ceil(flashOnTime .* devFrameRate);
postFlashDelay_frames = ceil(postFlashDelay .* devFrameRate);
targsAndFP_frames = ceil(targsAndFPTime .* devFrameRate);
targsOnly_frames = ceil(targsOnTime .* devFrameRate);
postChoiceFrames = 4; %so that we catch the saccade stabalizing
totNumFrames = sum([preFPFrames;...
                    preFlashFP_frames;...
                    flashOn_frames;...
                    postFlashDelay_frames;...
                    targsAndFP_frames;...
                    targsOnly_frames;...
                    postChoiceFrames]);

%% playing arround with making movie frames from a PTB texture frame
img = ones(768, 1024, 3) .* repmat(permute(bkgndrgb, [3,1,2]), [768, 1024,1]); %all pix at bkgnd
imCentX = round(size(img,2) ./ 2);
imCentY = round(size(img,1) ./ 2);

%compile the eyeposition data. I'm going to interpolate the given data
%based off the number of frames in the stimulus.
anlgStoreRate = 1./stro.sum.analog.storeRates{1}; % in seconds
elapsedTime = length(stro.ras{trialNum, 2}) .* anlgStoreRate;
anlgTimeVec = stro.ras{trialNum, 4}:anlgStoreRate:(stro.ras{trialNum, 4}+elapsedTime)-anlgStoreRate;
movieStartTime = stro.trial(trialNum, fpOnInd) - (preFPFrames ./ devFrameRate);
movieTimePts = [movieStartTime : 1./devFrameRate : movieStartTime + ((totNumFrames-1).*(1./devFrameRate))];
for a = 1:totNumFrames;
    errs = abs(anlgTimeVec - movieTimePts(a));
    ind = find(errs == min(errs));
    eyeX(a) = (stro.ras{trialNum, 2}(ind) .* 10) + eyeOffsetX; % in dva
    eyeY(a) = (stro.ras{trialNum, 3}(ind) .* 10) + eyeOffsetY; % in dva
end
eyeX = imCentX + round(eyeX .* pixperdeg);
eyeY = imCentY - round(eyeY .* pixperdeg);
eyeX = min(eyeX, 1014);     % don't let the eye position fall off the monintor
eyeX = max(eyeX, 10);       % and leave a bit of room for the eyePos marker
eyeY = min(eyeY, 758);
eyeY = max(eyeY, 10);
eyeImg = ones(epsize,epsize,3);
eyeImg = eyeImg .* repmat(permute([1,1,0], [3,1,2]), [size(eyeImg,1), size(eyeImg,2)]);
cometImg = ones(1,cometWidth,3) .* repmat(permute([1,1,0], [3,1,2]), 1, cometWidth);

% make the pre FP frames
clear M;
cometCounter = 1;
for a = 1:preFPFrames
    tmp_img = img;
    tmp_img(eyeY(a)-epHalfSize:eyeY(a)+epHalfSize, eyeX(a)-epHalfSize:eyeX(a)+epHalfSize, :) = eyeImg;
    if cometCounter>1
        dy = eyeY(cometCounter)-eyeY(cometCounter-1);
        dx = eyeX(cometCounter)-eyeX(cometCounter-1);
        if any([dy, dx])
            m = dy./dx;
            if ~isinf(m)
                b = eyeY(cometCounter) - m.*eyeX(cometCounter);
                lineX = min([eyeX(cometCounter-1),eyeX(cometCounter)]): max([eyeX(cometCounter-1),eyeX(cometCounter)]);
                for i = 1:length(lineX)
                    lineY = round(m.*lineX(i)+b);
                    tmp_img(lineY, lineX(i)-cometHalfWidth:lineX(i)+cometHalfWidth, :) = cometImg;
                end
            else
                lineY = min([eyeY(cometCounter-1),eyeY(cometCounter)]):max([eyeY(cometCounter-1),eyeY(cometCounter)]);
                for i = 1:length(lineY)
                    tmp_img(lineY(i), eyeX(cometCounter)-cometHalfWidth:eyeX(cometCounter)+cometHalfWidth, :) = cometImg;
                end
            end
        end
    end
    cometCounter = cometCounter+1;
    M(a) = im2frame(tmp_img);
end

% make the FP only frames
eyeInd = length(M);
for a = 1:preFlashFP_frames
    tmp_img = img;
    tmp_img([imCentY-fpHalfSize:imCentY+fpHalfSize], [imCentX-fpHalfSize:imCentX+fpHalfSize],:) = 0;
    tmp_img(eyeY(eyeInd+a)-epHalfSize:eyeY(eyeInd+a)+epHalfSize, eyeX(eyeInd+a)-epHalfSize:eyeX(eyeInd+a)+epHalfSize, :) = eyeImg;
    dy = eyeY(cometCounter)-eyeY(cometCounter-1);
    dx = eyeX(cometCounter)-eyeX(cometCounter-1);
    if any([dy, dx])
        m = dy./dx;
        if ~isinf(m)
                b = eyeY(cometCounter) - m.*eyeX(cometCounter);
                lineX = min([eyeX(cometCounter-1),eyeX(cometCounter)]): max([eyeX(cometCounter-1),eyeX(cometCounter)]);
                for i = 1:length(lineX)
                    lineY = round(m.*lineX(i)+b);
                    tmp_img(lineY, lineX(i)-cometHalfWidth:lineX(i)+cometHalfWidth, :) = cometImg;
                end
            else
                lineY = min([eyeY(cometCounter-1),eyeY(cometCounter)]):max([eyeY(cometCounter-1),eyeY(cometCounter)]);
                for i = 1:length(lineY)
                    tmp_img(lineY(i), eyeX(cometCounter)-cometHalfWidth:eyeX(cometCounter)+cometHalfWidth, :) = cometImg;
                end
            end
    end
    cometCounter = cometCounter+1;
    M(end+1) = im2frame(tmp_img);
end

% add the drifting gabor
phase = 0;
sigmaInPix = round((sigma./10) .* pixperdeg);
gbrX = round((flashX./10) .* pixperdeg);
gbrY = round((flashY./10) .* pixperdeg);
rampLength = ceil(flashOn_frames ./ 4);
ramp = linspace(0, 1, rampLength);  %ramp is 1/4th of the total duration on either side
plateau = ones(1,flashOn_frames - (rampLength .* 2));
timeProfile = [ramp, plateau, fliplr(ramp)];
flashSize = round((sigma ./ 10 .* nSigmas .* 2) .* pixperdeg); %convert from deg to pix.
halfSize = round(flashSize./2);
row = -halfSize:halfSize;
col = -halfSize:halfSize;
[X, Y] = meshgrid(row, col);
rgb_increment = [.09, -.09, 0]; %totally random for now
rgb_increment = permute(rgb_increment, [3,1,2]);
deltaPhase = driftRate .* 2 .* pi .* (1./devFrameRate);
eyeInd = length(M);
for a = 1:flashOn_frames;
    xprime = X .* cos(-theta) + Y .* sin(-theta);
    yprime = -X .* sin(-theta) + Y .* cos(-theta);
    gabor = exp(-(xprime.^2 + gamma.^2 .* yprime.^2) ./ (2.* sigmaInPix.^2)) .* cos(2 .* pi .* yprime ./ lambda + phase);
    gabor = repmat(gabor, [1,1,3]);
    gabor = gabor .* timeProfile(a) .* repmat(rgb_increment, [size(gabor,1), size(gabor,2),1]);
    gabor = gabor + repmat(permute(bkgndrgb, [3,1,2]), [size(gabor,1), size(gabor,2),1]);
    phase = phase+deltaPhase; %increment the phase
    
    %compile the image
    tmp_img = img;
    tmp_img([imCentY-fpHalfSize:imCentY+fpHalfSize], [imCentX-fpHalfSize:imCentX+fpHalfSize],:) = 0; %adding FP
    xRange = round(imCentX + [gbrX-size(gabor,2)./2 , gbrX+size(gabor,2)./2]);
    yRange = round(imCentY - [gbrY+size(gabor,1)./2 , gbrY-size(gabor,1)./2]);
    tmp_img(yRange(1):yRange(2)-1, xRange(1):xRange(2)-1, :) = gabor; % add the gabor
    tmp_img(eyeY(eyeInd+a)-epHalfSize:eyeY(eyeInd+a)+epHalfSize, eyeX(eyeInd+a)-epHalfSize:eyeX(eyeInd+a)+epHalfSize, :) = eyeImg; % add the eye pos
    dy = eyeY(cometCounter)-eyeY(cometCounter-1);
    dx = eyeX(cometCounter)-eyeX(cometCounter-1);
    if any([dy, dx])
        m = dy./dx;
        if ~isinf(m)
                b = eyeY(cometCounter) - m.*eyeX(cometCounter);
                lineX = min([eyeX(cometCounter-1),eyeX(cometCounter)]): max([eyeX(cometCounter-1),eyeX(cometCounter)]);
                for i = 1:length(lineX)
                    lineY = round(m.*lineX(i)+b);
                    tmp_img(lineY, lineX(i)-cometHalfWidth:lineX(i)+cometHalfWidth, :) = cometImg;
                end
            else
                lineY = min([eyeY(cometCounter-1),eyeY(cometCounter)]):max([eyeY(cometCounter-1),eyeY(cometCounter)]);
                for i = 1:length(lineY)
                    tmp_img(lineY(i), eyeX(cometCounter)-cometHalfWidth:eyeX(cometCounter)+cometHalfWidth, :) = cometImg;
                end
            end
    end
    
    cometCounter = cometCounter+1;
    M(end+1) = im2frame(tmp_img);
end

% pre targ wait frames. FP only
eyeInd = length(M);
for a = 1:postFlashDelay_frames+targsAndFP_frames;
    tmp_img = img;
    tmp_img([imCentY-fpHalfSize:imCentY+fpHalfSize], [imCentX-fpHalfSize:imCentX+fpHalfSize],:) = 0; %adding FP
    tmp_img(eyeY(eyeInd+a)-epHalfSize:eyeY(eyeInd+a)+epHalfSize, eyeX(eyeInd+a)-epHalfSize:eyeX(eyeInd+a)+epHalfSize, :) = eyeImg; % add the eye pos
    dy = eyeY(cometCounter)-eyeY(cometCounter-1);
    dx = eyeX(cometCounter)-eyeX(cometCounter-1);
    if any([dy, dx])
        m = dy./dx;
        if ~isinf(m)
            b = eyeY(cometCounter) - m.*eyeX(cometCounter);
            lineX = min([eyeX(cometCounter-1),eyeX(cometCounter)]): max([eyeX(cometCounter-1),eyeX(cometCounter)]);
            for i = 1:length(lineX)
                lineY = round(m.*lineX(i)+b);
                tmp_img(lineY, lineX(i)-cometHalfWidth:lineX(i)+cometHalfWidth, :) = cometImg;
            end
        else
            lineY = min([eyeY(cometCounter-1),eyeY(cometCounter)]):max([eyeY(cometCounter-1),eyeY(cometCounter)]);
            for i = 1:length(lineY)
                tmp_img(lineY(i), eyeX(cometCounter)-cometHalfWidth:eyeX(cometCounter)+cometHalfWidth, :) = cometImg;
            end
        end
    end
    cometCounter = cometCounter+1;
    M(end+1) = im2frame(tmp_img);
end

%make the target only frames.
targX = flashX./2; %in tenths dva
targY = flashY./2; %in tenths dva
targX = round((targX./10) .* pixperdeg);
targY = round((targY./10) .* pixperdeg);
T1loc = round([imCentY-targY, imCentX+targX]);
T2loc = round([imCentY+targY, imCentX-targX]);
eyeInd = length(M);
for a = 1:targsOnly_frames;
    tmp_img = img;
    tmp_img([T1loc(1)-fpHalfSize:T1loc(1)+fpHalfSize], [T1loc(2)-fpHalfSize:T1loc(2)+fpHalfSize],:) = 0;
    tmp_img([T2loc(1)-fpHalfSize:T2loc(1)+fpHalfSize], [T2loc(2)-fpHalfSize:T2loc(2)+fpHalfSize],:) = 0;
    tmp_img(eyeY(eyeInd+a)-epHalfSize:eyeY(eyeInd+a)+epHalfSize, eyeX(eyeInd+a)-epHalfSize:eyeX(eyeInd+a)+epHalfSize, :) = eyeImg; % add the eye pos
    dy = eyeY(cometCounter)-eyeY(cometCounter-1);
    dx = eyeX(cometCounter)-eyeX(cometCounter-1);
    if any([dy, dx])
        m = dy./dx;
        if ~isinf(m)
                b = eyeY(cometCounter) - m.*eyeX(cometCounter);
                lineX = min([eyeX(cometCounter-1),eyeX(cometCounter)]): max([eyeX(cometCounter-1),eyeX(cometCounter)]);
                for i = 1:length(lineX)
                    lineY = round(m.*lineX(i)+b);
                    tmp_img(lineY, lineX(i)-cometHalfWidth:lineX(i)+cometHalfWidth, :) = cometImg;
                end
            else
                lineY = min([eyeY(cometCounter-1),eyeY(cometCounter)]):max([eyeY(cometCounter-1),eyeY(cometCounter)]);
                for i = 1:length(lineY)
                    tmp_img(lineY(i), eyeX(cometCounter)-cometHalfWidth:eyeX(cometCounter)+cometHalfWidth, :) = cometImg;
                end
            end
    end
    cometCounter = cometCounter+1;
    M(end+1) = im2frame(tmp_img);
end

%make the post trial frames
eyeInd = length(M);
for a = 1:postChoiceFrames;
    tmp_img = img;
    tmp_img(eyeY(eyeInd+a)-epHalfSize:eyeY(eyeInd+a)+epHalfSize, eyeX(eyeInd+a)-epHalfSize:eyeX(eyeInd+a)+epHalfSize, :) = eyeImg; % add the eye pos
    dy = eyeY(cometCounter)-eyeY(cometCounter-1);
    dx = eyeX(cometCounter)-eyeX(cometCounter-1);
    if any([dy, dx])
        m = dy./dx;
        if ~isinf(m)
                b = eyeY(cometCounter) - m.*eyeX(cometCounter);
                lineX = min([eyeX(cometCounter-1),eyeX(cometCounter)]): max([eyeX(cometCounter-1),eyeX(cometCounter)]);
                for i = 1:length(lineX)
                    lineY = round(m.*lineX(i)+b);
                    tmp_img(lineY, lineX(i)-cometHalfWidth:lineX(i)+cometHalfWidth, :) = cometImg;
                end
            else
                lineY = min([eyeY(cometCounter-1),eyeY(cometCounter)]):max([eyeY(cometCounter-1),eyeY(cometCounter)]);
                for i = 1:length(lineY)
                    tmp_img(lineY(i), eyeX(cometCounter)-cometHalfWidth:eyeX(cometCounter)+cometHalfWidth, :) = cometImg;
                end
            end
    end
    cometCounter = cometCounter+1;
    M(end+1) = im2frame(tmp_img);
end




%% COMPILE THE MOVIE WITH MPGWRITE
repeat = 1;     %default = 1
pSearch = 2;    %default = 0
bSearch = 2;    %default = 1
reference = 1;  %default = 0
pixRange = 3;  %default = 10
iFrame = 8;     %default = 8
pFrame = 10;    %default = 10
bFrame = 25;    %default = 25

options = [repeat, pSearch, bSearch, reference, pixRange, iFrame, pFrame, bFrame];
mpgwrite(M, gray, 'testmpg.mpg', options)

%% Trying to make the sound file.
Fs = 42000; %for compatability with iMovie?

%define what spikes should sound like...
sampsPerSpike = 10;
%soundClip = exp((1-[1:sampsPerSpike])./3);
soundClip = [0.99, 0.99, 0.99, 0.99, 0.99, -0.99, -0.99, -0.99, -0.99, -0.99];

%create a reward beep
beepDuration = .200;                % in seconds
beepSamples = Fs .* beepDuration;
beepFreq = 900;                     % in cyc/sec (Hz)...
cycPerBeep = beepFreq .* beepDuration;
beep = linspace(0, cycPerBeep.*2.*pi, beepSamples);
beep = sin(beep) .* 0.99;

%compile the sound
movieDuration = totNumFrames .* (1./devFrameRate);   %in seconds
nSamples = round(Fs .* movieDuration);
spikeTimes = stro.ras{trialNum, 1} - movieStartTime;
soundSamples = zeros(1, nSamples);
for a = 1:length(spikeTimes);
    idx = round(Fs .* spikeTimes(a));
    soundSamples(idx:idx+sampsPerSpike-1) = soundClip;
end
rewTime = stro.trial(trialNum, rewOnInd) - movieStartTime;
idx = round(Fs .* rewTime);
soundSamples(idx:idx+beepSamples-1) = beep;

%preview the sound and write it to a .wav file
sound(soundSamples, Fs)
wavwrite(soundSamples, Fs, 'testsound')



