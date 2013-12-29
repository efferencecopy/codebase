%#ok<*DEFNU>
function IsoSampSlave()
    % Communication variables
    global udpCom KEY gl
    udpCom.port = 6665;
    udpCom.sock = nan;
    udpCom.rexip = '192.168.1.120';

    % Keypress constants
    KbName('UnifyKeyNames');
    KEY.ESC = KbName('ESCAPE');
    
    gl.fliprequest = 0;
    gl.timetoleave = 0;
    gl.framecounter = 0;
    gl.framecountermax = 0;

    gl.fp.on = 0;
    gl.fp.x = 0;
    gl.fp.y = 0;
    gl.fp.size = 0;
    gl.fp.rgb = [0 0 0];
    gl.fp.drawrect = [0 0 0 0];

    gl.stim.on = 0;
    gl.stim.template = [];
    
    gl.vpixx = ConnectedToViewPixx();
    gl.ccmode = 1; % in mode '1': every 2nd column of pixels is ignored

    [udpCom.sock, okay] = pnetStart(udpCom.port);
    if ~okay, return; end

    pnet(udpCom.sock, 'setreadtimeout', 0);
    pnet(udpCom.sock, 'setwritetimeout', 0);
    disp('In IsoSampSlave');

    messageIsAvailable = 0;
    while ~messageIsAvailable  % endless loop.
        messageIsAvailable = pnet(udpCom.sock, 'readpacket', ...
                                  1000, 'noblock');

        if messageIsAvailable
            DealWithMessage(messageIsAvailable);
            messageIsAvailable = 0;
        end

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

        if gl.timetoleave, return; end

        [keyisdown,junk,keycode] = KbCheck(); %#ok<ASGLU>
        if keyisdown && keycode(KEY.ESC)
            ShowCursor();
            pnet(udpCom.sock, 'close');
            sca();
            return
        end
    end
end

function DoFlip()
    global gl udpCom

    Screen('Flip', gl.windowPtr);
    if gl.stim.on
        if ~gl.framecounter
           pnet(udpCom.sock, 'write', 'MACSTIMON>> >>');
           pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
        end
        if gl.framecounter == gl.framecountermax - 1
            pnet(udpCom.sock, 'write', 'MACSTIMOFF>> >>');
            pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
            gl.stim.on = 0;
        end
        gl.framecounter = gl.framecounter + 1;
    end
    gl.fliprequest = 0;
end

function InitDisplay(mondist, screenwidth, varargin)
    global gl

    calfilename = varargin{1};
    load(calfilename);
    cal = cals{end}; %#ok<USENS> it comes from loading the above file
    gl.mondistcm = mondist;
    gl.screenWidthcm = screenwidth;
    gl.bkgndRGB = round(255 * cal.bgColor)';
    gl.cal.gammaTable = cal.gammaTable;
    gl.cal.monSpd = cal.P_device;
    gl.cal.invgammaTable = InvertGammaTable(cal.gammaInput, cal.gammaTable, 2^16);
   
    if nargin > 3 && ~isempty(varargin{2})
        fundfilename = varargin{2};
        s = load(fundfilename);
        fns = fieldnames(s);
        gl.cal.fundamentals = s.(fns{1})';
        wavelength_spacing = s.(fns{2});
        P_device = SplineSpd(SToWls(cal.S_device), cal.P_device, SToWls(wavelength_spacing));
        gl.cal.M = gl.cal.fundamentals'*P_device;
        gl.cal.invM = inv(gl.cal.M);
    end
    
    % start up the imaging pipeline
    if ~isempty(Screen('Windows'))
        gl.windowPtr = max(Screen('Windows'));
        Screen('FillRect', gl.windowPtr, cal.bgColor);
    else
        PsychImaging('PrepareConfiguration');
        if gl.vpixx
            PsychImaging('AddTask', 'General', 'EnableDataPixxC48Output', gl.ccmode);
        else
            PsychImaging('AddTask', 'General', 'EnableBits++Color++Output', gl.ccmode);
        end
        gl.windowPtr = PsychImaging('OpenWindow', 0, cal.bgColor);
    end

    gl.framerate = Screen('FrameRate', gl.windowPtr, 1);
    [screenwidthpix, screenheightpix]  = Screen('WindowSize', gl.windowPtr);
    gl.screenWidthpix = screenwidthpix;
    gl.screenHeightpix = screenheightpix;
    gl.screenCenterXpix = screenwidthpix / 2;
    gl.screenCenterYpix = screenheightpix / 2;

    pixpercm = gl.screenWidthpix / gl.screenWidthcm;
    theta = atan2(gl.screenWidthcm / 2, gl.mondistcm) * 180 / pi;
    cmperdeg = gl.screenWidthcm / 2 / theta;
    gl.pixperdeg = pixpercm * cmperdeg;

    gl.bkgndrgb = [gl.cal.gammaTable(gl.bkgndRGB(1) + 1,1)
                   gl.cal.gammaTable(gl.bkgndRGB(2) + 1,2)
                   gl.cal.gammaTable(gl.bkgndRGB(3) + 1,3)];

    gl.stim.on = 0;
    gl.fp.on = 0;
    HideCursor();
end

function PrepareGabor(flash_time, theta, sf, phi, sigma, ggamma, tf, nsigmas)
global gl

lambda = 1 / sf;
nframes = ceil(gl.framerate * flash_time / 1000);
ramplength = ceil(nframes / 4);
ramp = linspace(0, 1, ramplength);
plateau = ones(1, nframes - 2 * ramplength);
temporalprofile = [ramp plateau fliplr(ramp)];
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

function ShowFP(x, y, size, fpr, fpg, fpb)
    global gl

    gl.fp.x = x / 10;
    gl.fp.y = y / 10;
    gl.fp.size = size / 10;
    gl.fp.rgb = [fpr fpg fpb]/255;

    xx = gl.fp.x * gl.pixperdeg;
    yy = gl.fp.y * gl.pixperdeg;

    fpsizeinpix = round(gl.pixperdeg * gl.fp.size);

    gl.fp.drawrect = [gl.screenCenterXpix + xx - floor(fpsizeinpix / 2) + 1, ...
                      gl.screenCenterYpix - yy - floor(fpsizeinpix / 2) + 1, ...
                      gl.screenCenterXpix + xx + ceil(fpsizeinpix / 2), ...
                      gl.screenCenterYpix - yy + ceil(fpsizeinpix / 2)];
    gl.fp.on = 1;
end

function DrawFP()
    global gl
    Screen('FillRect', gl.windowPtr, gl.fp.rgb, gl.fp.drawrect);
    gl.fliprequest = 1;
end

function HideFP()
    global gl
    gl.fp.on = 0;
    gl.fliprequest = 1;
end

function ShowStim(stimx, stimy, l, m, s)
    global gl

    stimx = stimx / 10;
    stimy = stimy / 10;
    stimconecontrast = [l m s] / 100; % make stim dim

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
        fprintf('Squeezing back into gamut by a factor of %f!\n', scalefactor);
        gl.stim.rgb = (gl.stim.rgb - gl.bkgndrgb) ./ scalefactor + gl.bkgndrgb;
    end

    gl.framecounter = 0;
    gl.framecountermax = size(gl.stim.template, 3);
    gl.stim.on = 1;
 end

function DrawStim()
global gl

NGAMMASTEPS = size(gl.cal.invgammaTable, 1);
frame = gl.stim.template(:,:,gl.framecounter + 1);
img = zeros(size(gl.stim.template, 1), size(gl.stim.template, 2), 3);

for plane = 1:3
    tmp = frame * (gl.stim.rgb(plane) - gl.bkgndrgb(plane)) + gl.bkgndrgb(plane);
    tmp = round(tmp * (NGAMMASTEPS - 1)) + 1;
    tmp = gl.cal.invgammaTable(tmp,plane);
    img(:,:,plane) = reshape(tmp, size(img(:,:,1)));
end

tex = Screen('MakeTexture', gl.windowPtr, TranslateToColourMode(img, [], gl.ccmode), [], [], 2);
Screen('DrawTexture', gl.windowPtr, tex, [], gl.drawrect, [], 0);
Screen('Close', tex);
gl.fliprequest = 1;
end

function HideStim()
    global gl
    gl.stim.on = 0;
    gl.fliprequest = 1;
end

function AllOff()
    global gl
    gl.stim.on = 0;
    gl.fp.on = 0;
    gl.fliprequest = 1;
end

function DealWithMessage(msgSize)
    global udpCom gl
    message = pnet(udpCom.sock, 'read', msgSize, 'char');

    if strncmp(message, 'return', 6)
        stk = dbstack();
        if ~strcmp(stk(end).name, mfilename)
            gl.timetoleave = 1;
        end
    end

    try
        eval(message);
    catch ME
        fprintf('Trouble with message: "%s"\n', message);
        disp(getReport(ME));
    end
end