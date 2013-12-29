% ConeIsoTestBits.m
% Script for making cone isolating static stimuli to be shown to a color
% blind observer to test out the monitor calibration quality.
%
% Adding stuff so that it works with the Bits++ in Colour mode
%
% GDLH 6/8/07

% Flag for making measurements with the PR-650
MEASURE = 1;

% 1 = Gaussian blob, 0 = square
GAUSSIAN = 0;

% Loading calibration data
CALPATH = '/Monitor Calibration/Monitor data/Dell 4';
CALFILE = 'Dell4BitsCal';
CALIDX = 3;
load ([CALPATH,'/',CALFILE]);
cal = cals{CALIDX};
MAXDAC = 65535;

% Display parameters
NPIXPERSIDE = 200;  % # of pixels on each side of the image
DISPLAYTIME = 3;    % % of seconds to display each image
LMS = [0 0 0; 1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1; 0 0 0];  % color direction of each stimulus

if (GAUSSIAN)
    imagevect = normpdf(linspace(-4,4,NPIXPERSIDE),0,1);
    imagemat = imagevect'*imagevect;
    imagemat = imagemat./max(imagemat(:));
else
    imagemat = ones(NPIXPERSIDE,NPIXPERSIDE);
end
imagemat = imagemat(:,1:2:end);  
% getting rid of every other column in the image because the Bits++ in
% colour mode doubles the width of every pixel.

% Preparing the M matrix and the inverse gamma table
load ('T_cones_smj.mat');
fundamentals = T_cones_smj;
cal = cals{CALIDX};
P_device = SplineSpd(SToWls(cal.S_device), cal.P_device, SToWls([380 5 81]));
M = fundamentals*P_device;
rgbAmbient = FindModelWeights(cal.P_ambient, cal.P_device);
invGammaTable = InvertGamma(cal, 1);

if (MEASURE)
    disp('Turn on Photometer and hit <return>');
    disp('Then focus on square and hit <return> again');
    pause();
    PR650init(1,0);
    spds = [];
end

% Remember, cal.bgColor contains the *voltages* that we use to get the right
% background color, not the intensities.  For now, I'm just going to 
% set the background to the color that it was during the calibration, 
% that is, represented with only 8 bit precision.
bkgndrgb = [cal.gammaTable(round(256*cal.bgColor(1)),1);...
            cal.gammaTable(round(256*cal.bgColor(2)),2);...
            cal.gammaTable(round(256*cal.bgColor(3)),3)];

% Figuring out what the rgbs should be for each stimulus
rgbs = [];
for i = 1:size(LMS,1)
   rgb = inv(M)*LMS(i,:)';
   undershoot = rgb./-bkgndrgb;
   overshoot = rgb./(1-bkgndrgb);
   factor = 1/(max([overshoot; undershoot]));
   rgb = rgb.*factor;   % Going out as far as we can in the color direction
                        % given by rgb, staying within the monitor gamut
   if (all(isnan(rgb)))
       rgb = [0 0 0]';
   end
   % No need to worry about the ambient; cal.bgColor already contains the
   % ambient built into it.
   rgbs(i,:) = rgb;
end

% Opening the window
w = Screen('OpenWindow',0,round(256*cal.bgColor));
HideCursor;
[winWidth, winHeight] = Screen('WindowSize',w);
winCenterX = winWidth/2;
winCenterY = winHeight/2;
drawrect = [-floor(NPIXPERSIDE/2) -floor(NPIXPERSIDE/2) ceil(NPIXPERSIDE/2) ceil(NPIXPERSIDE/2)];
drawrect = drawrect + [winCenterX winCenterY winCenterX winCenterY];
if (drawrect(1)/2 ~= floor(drawrect(1)/2))  % ensure that the image starts on an even pixel
    drawrect(1) = drawrect(1) + 1;
    drawrect(3) = drawrect(3) + 1;
end

% Don't use the gamma compensation stuff on the graphics card with the
% Bits++.  Doing so would presumably decrease the available color resolution (the
% digital output of the graphics card is only 8 bits/channel and changing
% the "NormalizedGammaTable" to anything but a straight line would
% presumably use these 8 bits inefficiently - a uniformly distributed input
% would not lead to a uniformly distributed output).
Screen('LoadNormalizedGammaTable',0,repmat(linspace(0,1,256)',1,3));

% Present a square for the PR650 to aim at.  Wait for a <return>.
if (MEASURE)
   Screen('FillRect',w,round(255*[cal.bgColor+unifrnd(-.2,.2,3,1)]), drawrect);
   Screen('Flip',w);
   pause();
end

% Main drawing loop
for i = 1:size(rgbs,1)
    rgbim = cat(3,rgbs(i,1)*imagemat,rgbs(i,2)*imagemat,rgbs(i,3)*imagemat);
    rgbim = rgbim+repmat(permute(bkgndrgb,[3 2 1]),NPIXPERSIDE, NPIXPERSIDE/2);
    RGBim = round(MAXDAC*rgb2RGB(rgbim, invGammaTable));   
    RGBim = TranslateToColourMode(RGBim,1);
    tex=Screen('MakeTexture', w, RGBim);
    Screen('DrawTexture',w,tex,[],drawrect,[],0);
    Screen('Flip',w);
    Screen('Close',tex);
    pause(DISPLAYTIME);
    if (MEASURE)
        spds(i,:) = PR650measspd;
    end
end
Screen('Close',w);
ShowCursor;

% Calculating the cone contrasts from each of the measured spds.
% Ideally, the contrast to the isolated cone will be large and the contrast
% to the other cones will be negligible.  I make a measurement of the
% background color at the begining and at the end of the procedure.  The first
% background measurement provides the baseline against which contrast is
% defined.  The second one shows how much contamination you get (how much
% deviation from perfect cone isolation you get) with a repeated
% measurement.
LMS = [];
if (MEASURE)
    whiteidx = 1;
    LMSbkgnd = fundamentals*spds(whiteidx,:)';
    for i = 1:size(rgbs,1)
        LMS(i,:) = fundamentals*spds(i,:)';
    end
    LMSbkgndmat = repmat(LMSbkgnd',size(LMS,1),1);
    conecontrasts = (LMS-LMSbkgndmat)./LMSbkgndmat;
    figure;
    bar(conecontrasts);
end