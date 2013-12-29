% Playing aorund with gratings
%orientation = pi/4;
%phase = 0;
%period = 10;
%a = mkgrating(20,20,orientation,period,phase);
% Averaging neighboring columns
%b = (a(:,1:2:end)+a(:,2:2:end))/2;
%b = round(b*65535/2+65535/2);
%c = TranslatetoColourmode(cat(3,b,b,b));
%image(c/256)

%%
%
% Playing around with stimuli where the cone contrast
% either stays fixed across a background change or does not.
%
% GDLH 10/26/08

try
MONDIST = 100; % cm
gl.screenWidthcm = 40; %  from REX
nstim = 2;
nreps = 10;
stimsizeinpix = 10;

load ('T_cones_smj.mat');
fundamentals = T_cones_smj;
load('/Monitor Calibration/Monitor data/Dell 4/Dell4BitsCal');
cal = cals{end};
P_device = SplineSpd(SToWls(cal.S_device), cal.P_device, SToWls([380 5 81]));
M = fundamentals*P_device;
invM = inv(M);

gl.bkgndRGB = round(255*cal.bgColor)';
gl.cal.gammaTable = cal.gammaTable;
gl.cal.invgammaTable = InvertGamma(gl.cal.gammaTable, 1);
gl.cal.monSpd = cal.P_device;
gl.mondistcm = MONDIST;

bkgndrgb = [gl.cal.gammaTable(gl.bkgndRGB(1)+1,1); ...
               gl.cal.gammaTable(gl.bkgndRGB(2)+1,2); ...
               gl.cal.gammaTable(gl.bkgndRGB(3)+1,3)];

bkgndlms = M*bkgndrgb;
   
% SCREEN stuff
gl.windowPtr = Screen('OpenWindow',0, 255*cal.bgColor);
[screenwidthpix, screenheightpix]  = Screen('WindowSize', gl.windowPtr);
gl.screenCenterXpix = screenwidthpix/2;
gl.screenCenterYpix = screenheightpix/2;
pixpercm = screenwidthpix/gl.screenWidthcm;
cmperdeg = gl.screenWidthcm/(2*atan2(gl.screenWidthcm/2, gl.mondistcm)*180/pi);
gl.pixperdeg = pixpercm*cmperdeg;  % Small pixels
gl.framerate = Screen('FrameRate', gl.windowPtr, 1);

clut = repmat(linspace(0,1,256),3,1)';
Screen('LoadNormalizedGammaTable', 0, clut);

% Creating drawing windows
for j = 1:nstim
    x = 5*cos((j-1)/(nstim)*2*pi);
    y = 5*sin((j-1)/(nstim)*2*pi);
    drawrect = round([gl.screenCenterXpix+(x*gl.pixperdeg)-stimsizeinpix ...
                gl.screenCenterYpix-(y*gl.pixperdeg)-stimsizeinpix ...
                gl.screenCenterXpix+(x*gl.pixperdeg)+stimsizeinpix...
                gl.screenCenterYpix-(y*gl.pixperdeg)+stimsizeinpix]);

    if(rem(drawrect(1), 2)) %if the rectangle starts on an odd pixel
        drawrect(1) = drawrect(1) - 1;
        drawrect(3) = drawrect(3) - 1;
    end
    drawrects(j,:) = drawrect;
end

for i = 1:nreps  % repeats
% picking colors for stimuli
stimcc = []; stimlms = [];
for j = 1:nstim
    a = normrnd(0,.02);
    b = normrnd(0,.3);
    stimcc(j,:) = [a+.1 -a+.1 b+.1];
    stimlms(j,:) = bkgndlms'.*(1+stimcc(j,:));
end
stimrgb = (inv(M)*stimlms')';

% Figuring out change in background
a = normrnd(0,.015);
b = normrnd(0,.1);
bkgnddeltacontrast = [a -a b];
newbkgndlms = bkgndlms'.*(1+bkgnddeltacontrast);
newbkgndrgb = inv(M)*newbkgndlms';

% Finding new colors for stimuli
alpha = 0;
newstimlms = repmat(newbkgndlms,nstim,1).*(1+stimcc);   % Changing cone exc. ratios
whichstim = unidrnd(nstim);
whichstim = 1;
newstimlms(whichstim,:) = alpha.*newstimlms(whichstim,:)+(1-alpha).*stimlms(whichstim,:);
newstimrgb = (inv(M)*newstimlms')';

% Making background screen 1
im = zeros(screenheightpix, screenwidthpix/2, 3);
for plane = 1:3
    tmp = bkgndrgb(plane);
    tmp = round(tmp*size(gl.cal.invgammaTable,1)-1)+1;
    tmp = gl.cal.invgammaTable(tmp, plane);
    tmp = round(tmp*(size(gl.cal.invgammaTable,1)-1));
    im(:,:,plane) = repmat(tmp, screenheightpix, screenwidthpix/2);
end
bkgndtex=Screen('MakeTexture', gl.windowPtr, TranslateToColourMode(im,1));

% Making stimulus textures
for j = 1:nstim
    im = zeros(2*stimsizeinpix, stimsizeinpix, 3);
    for plane = 1:3
        tmp = stimrgb(j,plane);
        tmp = round(tmp*size(gl.cal.invgammaTable,1)-1)+1;
        tmp = gl.cal.invgammaTable(tmp, plane);
        tmp = round(tmp*(size(gl.cal.invgammaTable,1)-1));
        im(:,:,plane) = repmat(tmp, 2*stimsizeinpix, stimsizeinpix);
    end
    tex(j)=Screen('MakeTexture', gl.windowPtr, TranslateToColourMode(im,1));
end

% Making screen 2
im = zeros(screenheightpix, screenwidthpix/2, 3);
for plane = 1:3
    tmp = newbkgndrgb(plane);
    tmp = round(tmp*size(gl.cal.invgammaTable,1)-1)+1;
    tmp = gl.cal.invgammaTable(tmp, plane);
    tmp = round(tmp*(size(gl.cal.invgammaTable,1)-1));
    im(:,:,plane) = repmat(tmp, screenheightpix, screenwidthpix/2);
end
newbkgndtex=Screen('MakeTexture', gl.windowPtr, TranslateToColourMode(im,1));

for j = 1:nstim
    im = zeros(2*stimsizeinpix, stimsizeinpix, 3);
    for plane = 1:3
        tmp = newstimrgb(j,plane);
        tmp = round(tmp*size(gl.cal.invgammaTable,1)-1)+1;
        tmp = gl.cal.invgammaTable(tmp, plane);
        tmp = round(tmp*(size(gl.cal.invgammaTable,1)-1));
        im(:,:,plane) = repmat(tmp, 2*stimsizeinpix, stimsizeinpix);
    end
    newtex(j)=Screen('MakeTexture', gl.windowPtr, TranslateToColourMode(im,1));
end
Screen('DrawTexture',gl.windowPtr,newbkgndtex,[],[],[],0);
Screen('Flip',gl.windowPtr,0,0);
pause(1);
Screen('DrawTexture',gl.windowPtr,newbkgndtex,[],[],[],0);
for j = 1:nstim
    Screen('DrawTexture',gl.windowPtr,newtex(j),[],drawrects(j,:),[],0);
    Screen('Close',newtex(j));
end
Screen('Close',newbkgndtex);
Screen('Flip',gl.windowPtr,0,0);
pause(2);
Screen('DrawTexture',gl.windowPtr,bkgndtex,[],[],[],0);
for j = 1:nstim
    Screen('DrawTexture',gl.windowPtr,tex(j),[],drawrects(j,:),[],0);
    Screen('Close',tex(j));
end
Screen('Flip',gl.windowPtr,0,0);
pause(1);
Screen('DrawTexture',gl.windowPtr,bkgndtex,[],[],[],0); 
Screen('Close',bkgndtex);
Screen('Flip',gl.windowPtr,0,0);
pause(1);
end
Screen('CloseAll');
catch
    lasterror
    sca;
end

%%
% Match to sample task with a background that can have a spatial gradient
%
% GDLH 10/27/08


try
MONDIST = 100; % cm
gl.screenWidthcm = 40; %  from REX
stimx = 5;
stimy = 0;
ngrads = 5;
gradpos = 4;
stimsizeinpix = 10;
samplecc = [.2 .2 .2];
samplergb = inv(M)*(bkgndlms'.*(1+samplecc))';

% Parameters that define a trial
% 1) side with target
% 2) gradient type
% 3) distractor type

trialtypes = fullfact([2,ngrads,ngrads]);

gradlimscc = [.01 -.01 0; 0 0 0];
gradmagnitudes = linspace(-1,1,ngrads);
gradccs = gradmagnitudes'*gradlimscc(1,:);
gradlmss = repmat(bkgndlms',ngrads,1).*(1+gradccs);
distractorlmss = gradlmss.*(1+repmat(samplecc,ngrads,1));  % Cone excitations for the distractors
distractorrgb = inv(M)*distractorlmss(whichdist,:)';

load ('T_cones_smj.mat');
fundamentals = T_cones_smj;
load('/Monitor Calibration/Monitor data/Dell 4/Dell4BitsCal');
cal = cals{end};
P_device = SplineSpd(SToWls(cal.S_device), cal.P_device, SToWls([380 5 81]));
M = fundamentals*P_device;
invM = inv(M);

gl.bkgndRGB = round(255*cal.bgColor)';
gl.cal.gammaTable = cal.gammaTable;
gl.cal.invgammaTable = InvertGamma(gl.cal.gammaTable, 1);
gl.cal.monSpd = cal.P_device;
gl.mondistcm = MONDIST;

bkgndrgb = [gl.cal.gammaTable(gl.bkgndRGB(1)+1,1); ...
               gl.cal.gammaTable(gl.bkgndRGB(2)+1,2); ...
               gl.cal.gammaTable(gl.bkgndRGB(3)+1,3)];

bkgndlms = M*bkgndrgb;

% SCREEN stuff
gl.windowPtr = Screen('OpenWindow',0, 255*cal.bgColor);
[screenwidthpix, screenheightpix]  = Screen('WindowSize', gl.windowPtr);
gl.screenCenterXpix = screenwidthpix/2;
gl.screenCenterYpix = screenheightpix/2;
pixpercm = screenwidthpix/gl.screenWidthcm;
cmperdeg = gl.screenWidthcm/(2*atan2(gl.screenWidthcm/2, gl.mondistcm)*180/pi);
gl.pixperdeg = pixpercm*cmperdeg;  % Small pixels

clut = repmat(linspace(0,1,256),3,1)';
Screen('LoadNormalizedGammaTable', 0, clut);

% Creating drawing widows for the targets
targpos = [0 0; stimx stimy; -stimx stimy];
for j = 1:3
    drawrect = round([gl.screenCenterXpix+(targpos(j,1)*gl.pixperdeg)-stimsizeinpix ...
                      gl.screenCenterYpix-(targpos(j,2)*gl.pixperdeg)-stimsizeinpix ...
                      gl.screenCenterXpix+(targpos(j,1)*gl.pixperdeg)+stimsizeinpix...
                      gl.screenCenterYpix-(targpos(j,2)*gl.pixperdeg)+stimsizeinpix]);

    if(rem(drawrect(1), 2)) %if the rectangle starts on an odd pixel
           drawrect(1) = drawrect(1) - 1;
           drawrect(3) = drawrect(3) - 1;
    end
    drawrects(j,:) = drawrect;
end

for i = 1:10
% Things that vary trial by trial
whichgrad = unidrnd(ngrads);
whichdist = unidrnd(ngrads);
whichside = unidrnd(2);

% Making background screen without the gradient
im = zeros(screenheightpix, screenwidthpix/2, 3);
for plane = 1:3
    tmp = bkgndrgb(plane);
    tmp = round(tmp*size(gl.cal.invgammaTable,1)-1)+1;
    tmp = gl.cal.invgammaTable(tmp, plane);
    tmp = round(tmp*(size(gl.cal.invgammaTable,1)-1));
    im(:,:,plane) = repmat(tmp, screenheightpix, screenwidthpix/2);
end
bkgndtex=Screen('MakeTexture', gl.windowPtr, TranslateToColourMode(im,1));

% Making background screen with the gradient
im = zeros(screenheightpix, screenwidthpix/2, 3);
gradientlms = [];
gradendpointsinpix = round([gl.screenCenterXpix-(gradpos*gl.pixperdeg),...
                            gl.screenCenterXpix+(gradpos*gl.pixperdeg)]./2);
for cone = 1:3
    tmp = [repmat(gradccs(whichgrad,cone),1,gradendpointsinpix(1)),...
           linspace(gradccs(whichgrad,cone),gradlimscc(2,cone),gradendpointsinpix(2)-gradendpointsinpix(1)),...
           repmat(0,1,screenwidthpix/2-gradendpointsinpix(2))];
    gradientlms(cone,:) = bkgndlms(cone)*(1+tmp);
end
if (whichside == 2);
    gradientlms = fliplr(gradientlms);
end
gradientrgb = inv(M)*gradientlms;
for plane = 1:3
    tmp = gradientrgb(plane,:);
    tmp = round(tmp*size(gl.cal.invgammaTable,1)-1)+1;
    tmp = gl.cal.invgammaTable(tmp, plane);
    tmp = round(tmp*(size(gl.cal.invgammaTable,1)-1));
    im(:,:,plane) = repmat(tmp', screenheightpix, 1);
end
bkgndgradtex=Screen('MakeTexture', gl.windowPtr, TranslateToColourMode(im,1));

% Making the sample
im = zeros(2*stimsizeinpix, stimsizeinpix, 3);
for plane = 1:3
    tmp = samplergb(plane);
    tmp = round(tmp*size(gl.cal.invgammaTable,1)-1)+1;
    tmp = gl.cal.invgammaTable(tmp, plane);
    tmp = round(tmp*(size(gl.cal.invgammaTable,1)-1));
    im(:,:,plane) = repmat(tmp, 2*stimsizeinpix, stimsizeinpix);
end
sampletex=Screen('MakeTexture', gl.windowPtr, TranslateToColourMode(im,1));

% Making the distractor
im = zeros(2*stimsizeinpix, stimsizeinpix, 3);
for plane = 1:3
    tmp = distractorrgb(plane);
    tmp = round(tmp*size(gl.cal.invgammaTable,1)-1)+1;
    tmp = gl.cal.invgammaTable(tmp, plane);
    tmp = round(tmp*(size(gl.cal.invgammaTable,1)-1));
    im(:,:,plane) = repmat(tmp, 2*stimsizeinpix, stimsizeinpix);
end
distractortex=Screen('MakeTexture', gl.windowPtr, TranslateToColourMode(im,1));

% drawing
Screen('DrawTexture',gl.windowPtr,bkgndtex,[],[],[],0);

Screen('Flip',gl.windowPtr,0,0);
pause(1);
Screen('DrawTexture',gl.windowPtr,bkgndtex,[],[],[],0);
Screen('DrawTexture',gl.windowPtr,sampletex,[],drawrects(1,:),[],0);
Screen('Close',bkgndtex);
Screen('Flip',gl.windowPtr,0,0);
pause(1);
Screen('DrawTexture',gl.windowPtr,bkgndgradtex,[],[],[],0);
if (whichside == 1)  % distractor on left
    Screen('DrawTexture',gl.windowPtr,sampletex,[],drawrects(2,:),[],0);
    Screen('DrawTexture',gl.windowPtr,distractortex,[],drawrects(3,:),[],0);
else
    Screen('DrawTexture',gl.windowPtr,sampletex,[],drawrects(3,:),[],0);
    Screen('DrawTexture',gl.windowPtr,distractortex,[],drawrects(2,:),[],0);
end
Screen('Close',bkgndgradtex);
Screen('Flip',gl.windowPtr,0,0);
pause(1);
end
Screen('CloseAll');
catch
    lasterror
    sca;
end
%%
% Here's a place where Patrick can work on creating a smoothed bar that
% drifts across the RF of a neuron.
% GDLH 9/20/11
USECOLOURMODE = 0;
MONDIST = 100; % cm
STIMSIZEINDEG = 4;
gl.screenWidthcm = 40; % cm
stimuluscc = [-.07 .07 0]; % cone contrast units
stimpos = [-3 0]; % in degrees of visual angle relative to the center of the screen

load ('T_cones_smj10.mat');
fundamentals = T_cones_smj10;
load('/Monitor Calibration/Monitor data/Dell 4/Dell4BitsCal');
cal = cals{end};
P_device = SplineSpd(SToWls(cal.S_device), cal.P_device, SToWls([380 5 81]));
M = fundamentals*P_device;

gl.bkgndRGB = round(255*cal.bgColor)';
gl.cal.gammaTable = cal.gammaTable;
gl.cal.invgammaTable = InvertGamma(gl.cal.gammaTable, USECOLOURMODE);
gl.cal.monSpd = cal.P_device;
gl.mondistcm = MONDIST;

bkgndrgb = [gl.cal.gammaTable(gl.bkgndRGB(1)+1,1); ...
               gl.cal.gammaTable(gl.bkgndRGB(2)+1,2); ...
               gl.cal.gammaTable(gl.bkgndRGB(3)+1,3)];
bkgndlms = M*bkgndrgb;




try
    % ---- This ia all boilerplate stuff that sets up for displaying ---
    % Note: "SCREEN" is a multi-purpose function that is at the heart of
    % the Psychophysics Toolbox. Type "help screen" at the command prompt for
    % more information.
    gl.windowPtr = Screen('OpenWindow',0, 255*cal.bgColor);
    [screenwidthpix, screenheightpix]  = Screen('WindowSize', gl.windowPtr);
    gl.screenCenterXpix = screenwidthpix/2;
    gl.screenCenterYpix = screenheightpix/2;
    pixpercm = screenwidthpix/gl.screenWidthcm;
    cmperdeg = gl.screenWidthcm/(2*atan2(gl.screenWidthcm/2, gl.mondistcm)*180/pi);
    gl.pixperdeg = pixpercm*cmperdeg;  % Number os pixels per degree of visual angle

    stimsizeinpix = round(gl.pixperdeg*STIMSIZEINDEG);

    % Creating the spatial structure of the stimulus
    %stimspatialprofile = normpdf(linspace(-3,3,stimsizeinpix),0,1)'*normpdf(linspace(-3,3,stimsizeinpix),0,1);
    %stimspatialprofile = stimspatialprofile./max(stimspatialprofile(:));
    
    % Added by JPW to attempt creation of a gabor
    gaborSigma = STIMSIZEINDEG*10; %Std in 1/10ths degree VA - seems suspicious... making sigma = the whole stim size * 10 works pretty well... -JPW
    gaborGamma = 1; %Ratio between x and y std
    gaborLambda = stimsizeinpix; %Spatial period (pixels/cycle)
    gaborTheta = pi/4; %Orientation in radians
    flashDeltaPhase = 0; 
    gaborPhase = 0;
    
    if isinteger(stimsizeinpix)
        halfSize = stimsizeinpix/2;
        row = -halfSize:halfSize-1; %subtract one so that you don't overlap the texture window;
        col = -halfSize:halfSize-1;
    else
        halfsize = (stimsizeinpix-1)/2;
        row = -halfsize:halfsize;
        col = -halfsize:halfsize;
    end
        
    [X, Y] = meshgrid(row, col);
    xprime = X .* cos(-gaborTheta) + Y .* sin(-gaborTheta);
    yprime = -X .* sin(-gaborTheta) + Y .* cos(-gaborTheta);
    
    gaussian = exp(-(xprime.^2 + gaborGamma.^2 .* yprime.^2)./ (2.* gaborSigma.^2));
    
    % Drawing loop
    %for i=1:300
    
    gabor = cos(2 .* pi .* yprime ./ gaborLambda + gaborPhase);
    stimspatialprofile = gabor;
    
    % Loading a linear normalized gamma table
    % We don't want the graphics card to deal with gamma compensation, we
    % have to do it ourselves if we want to take advantage of the Bits++
    clut = repmat(linspace(0,1,256),3,1)';
    Screen('LoadNormalizedGammaTable', 0, clut);

    % Defining a drawing rectangle where the stimulus will be shown
    % [left, bottom, right, top]
    drawrect = round([gl.screenCenterXpix+(stimpos(1)*gl.pixperdeg)-stimsizeinpix/2 ...
        gl.screenCenterYpix-(stimpos(2)*gl.pixperdeg)-stimsizeinpix/2 ...
        gl.screenCenterXpix+(stimpos(1)*gl.pixperdeg)+stimsizeinpix/2 ...
        gl.screenCenterYpix-(stimpos(2)*gl.pixperdeg)+stimsizeinpix/2]);
    % Shifting the drawing rectangle if it happens to start on an
    % odd-numbered pixel. This is only important when using the Bits++
    if(rem(drawrect(1), 2))
        drawrect(1) = drawrect(1) - 1;
        drawrect(3) = drawrect(3) - 1;
    end

    % Creating an NxMx3 image to be displayed.
    im = zeros(stimsizeinpix, stimsizeinpix, 3); % preallocating space
    stimulusrgb = inv(M)*(bkgndlms'.*(1+stimuluscc))'; % computing rgbs at the peak
    for gun = 1:3
        % rescaling stimspatialprofile to go between bkgndrgb and stimulusrgb (instead of between 0 and 1) 
        tmp = (stimulusrgb(gun)-bkgndrgb(gun))*stimspatialprofile+bkgndrgb(gun);
        % Rounding the intensity values at each pixel so that we can use
        % these values as indices into the inverse gamma table
        tmp = round(tmp*size(gl.cal.invgammaTable,1))+1;
        % Checking for out of gamut errors
        if (any(tmp(:) > size(gl.cal.invgammaTable,1)) | any(tmp(:) < 1))
            [gun min(tmp(:)) max(tmp(:))]
            error('Out of gamut error');
        end
        % converting intensities to voltages
        tmp = gl.cal.invgammaTable(tmp, gun);
        tmp = round(tmp*(size(gl.cal.invgammaTable,1)-1));
        im(:,:,gun) = reshape(tmp, stimsizeinpix, stimsizeinpix);
    end
    stimtex=Screen('MakeTexture', gl.windowPtr, im);

    % drawing
    for i = 1:300
        Screen('DrawTexture',gl.windowPtr,stimtex,[],drawrect+[i 0 i 0],[],0);
        %gaborPhase = gaborPhase + flashDeltaPhase;
        Screen('DrawTexture',gl.windowPtr,stimtex,[],drawrect,[],0);
        Screen('Flip',gl.windowPtr);
    end
    Screen('Close',stimtex);
    Screen('CloseAll');
catch
    sca;
    rethrow(lasterror);
end

%%
% Trying to make a random dot kinematogram 4/6/13 GDLH

try
    MONDIST = 100; % cm
    gl.screenWidthcm = 40; %  from REX
    load('/Monitor Calibration/Monitor data/Dell 4/Dell4BitsCal');
    cal = cals{end};
    gl.bkgndRGB = round(255*cal.bgColor)';
    gl.cal.gammaTable = cal.gammaTable;
    gl.cal.invgammaTable = InvertGamma(gl.cal.gammaTable, 1);
    gl.cal.monSpd = cal.P_device;
    gl.mondistcm = MONDIST;
    
    bkgndrgb = [gl.cal.gammaTable(gl.bkgndRGB(1)+1,1); ...
        gl.cal.gammaTable(gl.bkgndRGB(2)+1,2); ...
        gl.cal.gammaTable(gl.bkgndRGB(3)+1,3)];
    
    % SCREEN stuff
%    gl.windowPtr = Screen('OpenWindow',0, 255*cal.bgColor);
    % New way of setting up a screen
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask','General','UseFastOffscreenWindows');
    [gl.windowPtr, gl.winRect] = PsychImaging('OpenWindow', 0, 255*cal.bgColor, [], [], 2); % 2 buffers
    [screenwidthpix, screenheightpix]  = Screen('WindowSize', gl.windowPtr);
    gl.screenCenterXpix = screenwidthpix/2;
    gl.screenCenterYpix = screenheightpix/2;
    pixpercm = screenwidthpix/gl.screenWidthcm;
    cmperdeg = gl.screenWidthcm/(2*atan2(gl.screenWidthcm/2, gl.mondistcm)*180/pi);
    gl.pixperdeg = pixpercm*cmperdeg;  % Small pixels
    gl.framerate = Screen('FrameRate', gl.windowPtr, 1);
    
    % Taken from DotDemo.m
%    Screen('BlendFunction',gl.windowPtr, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    fps = Screen('FrameRate',gl.windowPtr);
%    ifi = Screen('GetFlipInterval',gl.windowPtr);
%    if fps==0
%        fps=1/ifi;
%    end
    vbl = Screen('Flip', gl.windowPtr);
    
    % Setting up the stimulus
    MAXNFRAMES = 200;
    gl.stim.diam = 3;
    gl.stim.x = 1;
    gl.stim.y = -2;
    gl.stim.speed = 2; % DVA/sec
    gl.stim.dotdiam = .05;
    gl.stim.ndots = 100;

    stimsizeinpix = round(gl.pixperdeg*gl.stim.diam/2);  % /2 counting in double-width pixels
    clut = repmat(linspace(0,1,256),3,1)';
    Screen('LoadNormalizedGammaTable', 0, clut);
    disp('Got here 1')    
%     drawrect = round([gl.screenCenterXpix+(gl.stim.x*gl.pixperdeg)-stimsizeinpix ...
%         gl.screenCenterYpix-(gl.stim.y*gl.pixperdeg)-stimsizeinpix ...
%         gl.screenCenterXpix+(gl.stim.x*gl.pixperdeg)+stimsizeinpix...
%         gl.screenCenterYpix-(gl.stim.y*gl.pixperdeg)+stimsizeinpix]);
%     
%     if(rem(drawrect(1), 2)) %if the rectangle starts on an odd pixel
%         drawrect(1) = drawrect(1) - 1;
%         drawrect(3) = drawrect(3) - 1;
%     end
%     
%     im = ones(stimsizeinpix*2, stimsizeinpix,3);
%     
%     [x,y] = meshgrid(linspace(-1,1,stimsizeinpix),linspace(-1,1,stimsizeinpix*2));
%     mask = x.^2. + y.^2 > 1;
%     for i = 1: 3
%         im(:,:,i) = im(:,:,i)*cal.bgColor(i);
%         im(:,:,i) = im(:,:,i) .* mask;
%     end
%     tex=Screen('MakeTexture', gl.windowPtr, TranslateToColourMode(im*2^16,1));
% 
%     Screen('DrawTexture',gl.windowPtr,tex,[],drawrect,[],0);
%     disp('got here 4')
%     Screen('Flip',gl.windowPtr,0,0);

    % now drawing some dots
    pfs = gl.stim.speed * gl.pixperdeg / fps;  % dot speed in pixels per frame
    s = gl.stim.dotdiam * gl.pixperdeg;    % Size of dots in pixels
    [dots_xy] = unidrnd(stimsizeinpix*2,2,gl.stim.ndots,4);
 %   txy = [gl.screenCenterXpix+(gl.stim.x*gl.pixperdeg) gl.screenCenterYpix-(gl.stim.y*gl.pixperdeg);...
 %       gl.screenCenterXpix-(gl.stim.x*gl.pixperdeg) gl.screenCenterYpix+(gl.stim.y*gl.pixperdeg);...
 %       gl.screenCenterXpix+(gl.stim.y*gl.pixperdeg) gl.screenCenterYpix+(gl.stim.x*gl.pixperdeg);...
 %       gl.screenCenterXpix-(gl.stim.y*gl.pixperdeg) gl.screenCenterYpix-(gl.stim.x*gl.pixperdeg)];
 
    % Sign flip on y values below to make the transformation
    % from DVA ('=' = up) to screen pixels ('+' = down)
    rotmat = [0 -1; 1 0]; % 90 degree rotation
    txy = [];
    for i = 1:4 % rotating by 90 degrees each time
        txy(i,:) = [gl.screenCenterXpix gl.screenCenterYpix]+gl.pixperdeg*[gl.stim.x -gl.stim.y]*rotmat^(i-1)';
    end
    
    distractorpositions = [1:4];
    signalposition = unidrnd(4);
    distractorpositions(distractorpositions == signalposition) = [];
    
    signaldir = rand<.5; % My convention: down is '0', up is '1'
    dirsign = -(2*signaldir-1); % dirsign: '-1' means signal up, '1' means signal down
    % The way the drawing works, the top of the screen are small y
    % values.
    for i = 1:MAXNFRAMES
        for j = 1:4
            L = sum((dots_xy(:,:,j)-stimsizeinpix).^2) < stimsizeinpix.^2;
            Screen('DrawDots',gl.windowPtr,dots_xy(:,L,j),s,[],txy(j,:))
        end
        Screen('Flip',gl.windowPtr,0,0);
        dots_xy(:,:,signalposition) = dots_xy(:,:,signalposition) + repmat([0;dirsign*pfs],[1 size(dots_xy,2) 1]); % first plane are signal dots
        dots_xy(:,:,distractorpositions) = dots_xy(:,:,distractorpositions) + repmat([0;-dirsign*pfs],[1 size(dots_xy,2) 3]);
        dots_xy(dots_xy(:) > stimsizeinpix*2) = 0;
        dots_xy(dots_xy(:) < 0) = stimsizeinpix*2;    
    end
    % Screen('Close',tex);
    sca;
catch
    disp('errored out');
    sca;
end

signalposition
signaldir