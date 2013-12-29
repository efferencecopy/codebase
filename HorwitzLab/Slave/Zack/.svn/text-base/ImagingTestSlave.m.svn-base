%#ok<*DEFNU>
function ImagingTestSlave()
% Communication variables
global udpCom gl
udpCom.port = 6665;
udpCom.sock = nan;
udpCom.rexip = '192.168.1.120';

% Variables that are destructively modified by subfunctions Because we are
% guaranteed to be running this code after having already setup the global
% gl structure (via the white noise paradigm) there's no sense in
% reassigning everything.
gl.vpixx = ConnectedToViewPixx();

gl.fliprequest = false;
gl.timetoleave = false;
gl.framerate = 0;
gl.framecounter = 0;
gl.framecountermax = 0;

gl.fp.on = false;
gl.fp.x = 0;
gl.fp.y = 0;
gl.fp.size = 0;
gl.fp.rgb = [0 0 0];
gl.fp.drawrect = [0 0 0 0];

gl.windowPtr = 0;

[udpCom.sock, success] = pnetStart(udpCom.port);
if ~success, return, end

pnet(udpCom.sock, 'setreadtimeout', 0);
pnet(udpCom.sock, 'setwritetimeout', 0);

messageIsAvailable = 0;
while ~messageIsAvailable  % endless loop.
    messageIsAvailable = pnet(udpCom.sock, 'readpacket', 1000, 'noblock');
    
    if messageIsAvailable
        DealWithMessage(messageIsAvailable);
        messageIsAvailable = 0;
    end
    
    if gl.timetoleave, return, end
    
    if gl.windowPtr > 0
        if gl.stim.on
            DrawStim();
        end
        if gl.fp.on
            DrawFP();
        end
        if gl.fliprequest
            DoFlip();
        end
    end
    
    [keyisdown,nil,keycode] = KbCheck();
    if keyisdown
        parse_key_press(keycode)
    end
end
end

function parse_key_press(keycode)
global gl udpCom
switch find(keycode,1)
    case KbName('escape')
        ShowCursor();
        pnet(udpCom.sock, 'close');
        sca();
        gl.timetoleave = true;
end
end

function DoFlip()
global gl udpCom

if gl.stim.on    
    Screen('Flip', gl.windowPtr);
    if ~gl.framecounter
       pnet(udpCom.sock, 'write', 'MACSTIMON>> >>');
       pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
    end
    if gl.framecounter == gl.framecountermax - 1
        pnet(udpCom.sock, 'write', 'MACSTIMOFF>> >>');
        pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
        Screen('Close', gl.tex);
        gl.stim.on = false;
    end
    gl.framecounter = gl.framecounter + 1;
end
gl.fliprequest = false;
end

% Varargin{1} and {2} are expected to be the calibration and fundamentals
% filename (or full paths).
function InitDisplay(mondist, screenwidth, varargin)
global gl

assert(nargin > 2 && ~isempty(varargin{1}), ...
    'You need to provide a calibration file name!');

calfilename = varargin{1};
load(calfilename);
cal = cals{end}; %#ok<USENS>
gl.mondistcm = mondist;
gl.screenWidthcm = screenwidth;
gl.bkgndRGB = round(255*cal.bgColor)';
gl.cal.gammaTable = cal.gammaTable;
gl.cal.monSpd = cal.P_device;
gl.cal.invgammaTable = InvertGammaTable(cal.gammaInput, gl.cal.gammaTable, 2^16);

if nargin > 3 && ~isempty(varargin{2})
    fundfilename = varargin{2};
    s = load(fundfilename);
    fns = fieldnames(s);
    gl.cal.fundamentals = s.(fns{1})';
    P_device = SplineSpd(SToWls(cal.S_device), cal.P_device, SToWls([380 5 81]));
    gl.cal.M = gl.cal.fundamentals'*P_device;
    gl.cal.invM = inv(gl.cal.M);
end

% A flag to control which stimulus drawing method to use (for comparison)
if nargin < 5 || isempty(varargin{3})
    gl.new_method = 0;
else
    gl.new_method = varargin{3};
end

gl.bkgndrgb = [gl.cal.gammaTable(gl.bkgndRGB(1)+1,1);...
    gl.cal.gammaTable(gl.bkgndRGB(2)+1,2);...
    gl.cal.gammaTable(gl.bkgndRGB(3)+1,3)];

% enable colour mode on either device
% from BitsPlusPlus.m: "Both the graphics hardwares gamma table and the
% Bits++ internal clut are set to an identity mapping while this mode is
% active." i.e., we don't need to invoke 'LoadNormalizedGammaTable'.
if gl.new_method
    gl.ccmode = 1; % in mode '1': every 2nd column of pixels is ignored
    % start up the imaging pipeline
    if ~isempty(Screen('Windows'))
        gl.windowPtr = max(Screen('Windows'));
        Screen('FillRect', gl.windowPtr, gl.bkgndRGB/255);
    else
        PsychImaging('PrepareConfiguration');
        if gl.vpixx
            PsychImaging('AddTask', 'General', 'EnableDataPixxC48Output',   gl.ccmode);
        else
            PsychImaging('AddTask', 'General', 'EnableBits++Color++Output', gl.ccmode);
        end
        gl.windowPtr = PsychImaging('OpenWindow', 0, gl.bkgndRGB/255);
    end
else
    Screen('LoadNormalizedGammaTable', 0, repmat(linspace(0,1,256),3,1)');
    gl.windowPtr = Screen('OpenWindow', 0, gl.bkgndRGB);
end

gl.framerate = Screen('NominalFrameRate', gl.windowPtr, 1);
[screenwidthpix, screenheightpix] = Screen('WindowSize', gl.windowPtr);
gl.screenWidthpix = screenwidthpix;
gl.screenHeightpix = screenheightpix;
gl.screenCenterXpix = screenwidthpix/2;
gl.screenCenterYpix = screenheightpix/2;

pixpercm = gl.screenWidthpix/gl.screenWidthcm;
theta = atan2(gl.screenWidthcm/2, gl.mondistcm)*180/pi;
cmperdeg = gl.screenWidthcm/2/theta;
gl.pixperdeg = pixpercm*cmperdeg;

gl.bkgndrgb = [gl.cal.gammaTable(gl.bkgndRGB(1)+1,1);...
    gl.cal.gammaTable(gl.bkgndRGB(2)+1,2);...
    gl.cal.gammaTable(gl.bkgndRGB(3)+1,3)];

gl.fp.on = false;
gl.stim.on = false;

HideCursor();
end

% "Show" functions are called from REX.
% They set up the appropriate fields in the "gl" structure (and set the
% "...on" toggle field to true).
function ShowFP(x, y, size, fpr, fpg, fpb)
global gl

gl.fp.x = x/10;
gl.fp.y = y/10;
gl.fp.size = size/10;
gl.fp.rgb = [fpr fpg fpb]/255;

xx = gl.fp.x*gl.pixperdeg;
yy = gl.fp.y*gl.pixperdeg;

fpsizeinpix = round(gl.pixperdeg*gl.fp.size);

gl.fp.drawrect = [gl.screenCenterXpix + xx - floor(fpsizeinpix/2) + 1, ...
    gl.screenCenterYpix - yy - floor(fpsizeinpix/2) + 1, ...
    gl.screenCenterXpix + xx + ceil(fpsizeinpix/2), ...
    gl.screenCenterYpix - yy + ceil(fpsizeinpix/2)];

gl.fp.on = true;
end

% "Draw" functions get called on every screen refresh so long as the
% corresponding graphical object is to be displayed.  These are the
% functions that actually draw the object to the screen (and they make a
% fliprequest so that the drawn object will actually appear).
function DrawFP()
global gl

Screen('FillRect', gl.windowPtr, gl.fp.rgb, gl.fp.drawrect);
gl.fliprequest = true;
end

function HideFP()
global gl

gl.fp.on = false;
gl.fliprequest = true;
end

% testing mock up:
% gl.pixperdeg=65.87; gl.framerate=120; nframesplateau=40; nframesramp=20; sigma=.4; nsigmas=3; tf=3; ggamma=1; sf=3; theta=0; phi=0;
function PrepareGabor(nframesplateau, nframesramp, theta, sf, phi, sigma, ggamma, tf, nsigmas)
global gl
sf = .5;
tf = 60;
lambda = 1 / sf;
ramp = linspace(0, 1, nframesramp);
plateau = ones(1, nframesplateau);
temporalprofile = [ramp plateau fliplr(ramp)];
nframes = length(temporalprofile);
stimsizeindeg = sigma * nsigmas; % This is half stim size (we we go out +/- nsigmas)
stimsizeinpix = round(stimsizeindeg * gl.pixperdeg);  % full stim size in doublewide pixels
[x,y] = meshgrid(stimsizeindeg * linspace(-1, 1, stimsizeinpix), ...
    stimsizeindeg * linspace(-1, 1, 2 * stimsizeinpix));
% x and y are in dva
X = x * cos(-theta) + y * sin(-theta);
Y =-x * sin(-theta) + y * cos(-theta);

deltaphase = tf * 2 * pi / gl.framerate;
phases = phi + (0:nframes-1)*deltaphase;
phases = reshape(phases, [1 1 nframes]);
temporalprofile = reshape(temporalprofile, [1 1 nframes]);
expterm = bsxfun(@times, exp(-(X.^2 + ggamma^2 * Y.^2) / 2 / sigma^2), temporalprofile);
costerm = cos(bsxfun(@plus, 2 * pi * Y / lambda, phases));
gl.stim.template = expterm .* costerm;
end

% test mockup:
% stimx = 8; stimy = -20; lcont = .6; mcont = .6; scont = .6;
% stimx = 8; stimy = -20; lcont = .09; mcont = -.09; scont = 0;
function ShowStim(stimx, stimy, stimconecontrast)
global gl

stimx = stimx / 10;
stimy = stimy / 10;

% Creating the drawing window
x = stimx * gl.pixperdeg;
y = stimy * gl.pixperdeg;
stimsizeinpix = size(gl.stim.template, 2); % in doublewide pix

gl.drawrect = round([gl.screenCenterXpix + x - stimsizeinpix ...
    gl.screenCenterYpix - y - stimsizeinpix ...
    gl.screenCenterXpix + x + stimsizeinpix ...
    gl.screenCenterYpix - y + stimsizeinpix]);

if rem(gl.drawrect(1), 2) %if the rectangle starts on an odd pixel
    gl.drawrect(1) = gl.drawrect(1) - 1;
    gl.drawrect(3) = gl.drawrect(3) - 1;
end

assert(gl.drawrect(3) - gl.drawrect(1) == 2 * stimsizeinpix ...
    && gl.drawrect(4) - gl.drawrect(2) == 2 * stimsizeinpix, ...
    'Incorrect draw window size!');

% Getting RGBs
bkgndlms = gl.cal.M * gl.bkgndrgb;
gl.stim.rgb = gl.cal.invM * (bkgndlms .* (1 + stimconecontrast'));

% Checking for out of gamut errors
% Right now just squeezing it back into the gamut - need to send a
% message to REX that this was done.
lims = min([gl.bkgndrgb'; 1 - gl.bkgndrgb']);
if any(abs(gl.stim.rgb - gl.bkgndrgb) > lims')
    scalefactor = max(abs(gl.stim.rgb - gl.bkgndrgb) ./ lims');  % +(1/NGAMMASTEPS);
    gl.stim.rgb = (gl.stim.rgb - gl.bkgndrgb) ./ scalefactor + gl.bkgndrgb;
end

gl.framecountermax = size(gl.stim.template, 3);

gl.tex = zeros(1,gl.framecountermax);

NGAMMASTEPS = size(gl.cal.invgammaTable, 1);
for i=1:gl.framecountermax
    frame = gl.stim.template(:,:,i);
    img = zeros(size(gl.stim.template, 1), size(gl.stim.template, 2), 3);
    
    for plane = 1:3
        tmp = frame * (gl.stim.rgb(plane) - gl.bkgndrgb(plane)) + gl.bkgndrgb(plane);
        tmp = round(tmp * (NGAMMASTEPS - 1)) + 1;
        tmp = gl.cal.invgammaTable(tmp,plane);
        img(:,:,plane) = reshape(tmp, size(img(:,:,1)));
    end
    
    if ~gl.new_method
        img = round(img * (NGAMMASTEPS - 1));
        gl.tex(i) = Screen('MakeTexture', gl.windowPtr, TranslateToColourMode(img));
    else
        gl.tex(i) = Screen('MakeTexture', gl.windowPtr, TranslateToColourMode(img, [], gl.ccmode), [], [], 2);
    end
    %Screen('Close', gl.tex(i));
end

gl.framecounter = 0;

gl.stim.on = true;
end

function DrawStim()
global gl

Screen('DrawTexture', gl.windowPtr, gl.tex(gl.framecounter+1), [], gl.drawrect, [], 0);
% Screen('Close', gl.tex(gl.framecounter+1));
gl.fliprequest = true;
end

function HideStim()
global gl

gl.stim.on = false;
gl.fliprequest = true;
end

function AllOff()
global gl

gl.stim.on = false;
gl.fp.on = false;
gl.fliprequest = true;
end

function DealWithMessage(msgSize)
global udpCom gl

message = pnet(udpCom.sock, 'read', msgSize, 'char');
if strncmp(message, 'return', 6)
    stk = dbstack();  % Check whether called from another function or from command line
    if ~strcmp(stk(end).name, mfilename)
        gl.timetoleave = true;
    end
end

try
    eval(message);
catch ME
    fprintf('Trouble with message: "%s"\n', message);
    disp(getReport(ME));
end
end