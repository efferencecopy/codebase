% ConeIsoTest.m
% Script for making cone isolating static stimuli to be shown to a color
% blind observer to test out the monitor calibration quality.
%
% See "ConeIsoTestBits.m" for a veriosn that works with the Bits++ in
% Colour mode.
%
% GDLH 6/8/07

% Flag for making measurements with the PR-650
MEASURE = 0;

% Loading calibration data
CALPATH = '/Monitor Calibration/Monitor data/Dell 4';
CALFILE = 'Dell4Cal';
CALIDX = 3;
load ([CALPATH,'/',CALFILE]);
cal = cals{CALIDX};

% Display parameters
NPIXPERSIDE = 50;
DISPLAYTIME = 5;
LMS = [0 0 0; 1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 2; 0 0 -2; 0 0 0];

% Setting up the NormalizedGammaTable
bkgndrgb = cal.bgColor;
oldGammaTable = Screen('ReadNormalizedGammaTable',0);
newGammaTable = InvertGamma(cal, 0);

% Preparing the M matrix
load ('T_cones_sp.mat');
fundamentals = T_cones_smj;
cal = cals{CALIDX};
P_device = SplineSpd(SToWls(cal.S_device), cal.P_device, SToWls([380 5 81]));
M = fundamentals*P_device;
rgbAmbient = FindModelWeights(cal.P_ambient, cal.P_device);

% Computing RGBs
RGBs = [];
for i = 1:size(LMS,1)
   rgb = inv(M)*LMS(i,:)';
   undershoot = rgb./-cal.bgColor;
   overshoot = rgb./(1-cal.bgColor);
   factor = 1/(max([overshoot; undershoot]));
   rgb = rgb.*factor;
   if (all(isnan(rgb)))
       rgb = [0 0 0]';
   end
   % No need to worry about the ambient; cal.bgColor already contains the
   % ambient built into it.
   rgb = rgb+cal.bgColor;
   RGBs(i,:) = rgb';
end

if (MEASURE)
    disp('Turn on Photometer and hit <return>');
    disp('Then focus on square and hit <return> again');
    pause();
    PR650init(1,0);
    spds = [];
end

% Opening the window
w = Screen('OpenWindow',0,round(255*cal.bgColor));
HideCursor;
[winWidth, winHeight] = Screen('WindowSize',w);

winCenterX = winWidth/2+1;
winCenterY = winHeight/2+1;
drawrect = [-floor(NPIXPERSIDE/2) -floor(NPIXPERSIDE/2) ceil(NPIXPERSIDE/2) ceil(NPIXPERSIDE/2)];
drawrect = drawrect + [winCenterX winCenterY winCenterX winCenterY];
Screen('LoadNormalizedGammaTable',0,newGammaTable);

if (MEASURE)
   Screen('FillRect',w,round(255*[cal.bgColor+unifrnd(-.2,.2,3,1)]), drawrect);
   Screen('Flip',w);
end
pause();  % initial pause - hit return to get out

% Main drawing loop
for i = 1:size(RGBs,1)
    rgb = RGBs(i,:);
    image = 255*repmat(permute(rgb,[3 1 2]),NPIXPERSIDE,NPIXPERSIDE);
    tex=Screen('MakeTexture', w, image);
    Screen ('DrawTexture',w,tex,[],drawrect,[],1);
    Screen('Flip',w);
    Screen('Close',tex);
    pause(DISPLAYTIME);
    if (MEASURE)
        spds(i,:) = PR650measspd;
    end
end
Screen('Close',w);
ShowCursor;
Screen('LoadNormalizedGammaTable',0,oldGammaTable);

LMS = [];
if (MEASURE)
    whiteidx = find(all(RGBs==repmat(cal.bgColor',size(RGBs,1),1),2), 1);
    LMSbkgnd = fundamentals*spds(whiteidx,:)';
    for i = 1:size(RGBs,1)
        LMS(i,:) = fundamentals*spds(i,:)';
    end
    LMSbkgndmat = repmat(LMSbkgnd',size(LMS,1),1);
    conecontrasts = (LMS-LMSbkgndmat)./LMSbkgndmat;
    bar(conecontrasts);
end