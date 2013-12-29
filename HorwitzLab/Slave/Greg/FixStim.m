%#ok<*DEFNU>
function FixStim
    % Communication variables
    global udpCom
    udpCom.port = 6665;
    udpCom.sock = nan;
    udpCom.rexip = '192.168.1.120';

    % Keypress constants
    global KEY
    KbName('UnifyKeyNames');
    KEY.ESC = KbName('escape');
    
    global gl
    % Variables that are destructively modified by subfunctions
    % Because we are guaranteed to be running this code after having
    % already
    % setup the global gl structure (via the white noise paradigm)
    % there's no sense in reassigning everything.
    gl.usecolourmode = 1;  % Are we using the Bits++
    gl.vpixx = ConnectedToViewPixx();
    gl.ccmode = 1;
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

    gl.targ.on = 0;
    gl.targ.x = 0;  % Position of RF x
    gl.targ.y = 0;  % Position of RF y
    gl.targ.size = 0;
    gl.targ.rgb = [0 0 0];
    gl.targ.drawrect = [];
    
    gl.dots.on = 0;
    gl.dots.apsizeinpix = 0;
    gl.dots.dotsizeinpix = 0;
    gl.dots.ppf = 0; % pixels per frame
    gl.dots.ppfdist = 0; % for training 
    gl.dots.positions = []; % By convention the first element is signal position
    gl.dots.dir = 0;
    
    [udpCom.sock, Success] = pnetStart(udpCom.port);
    if ~Success, return, end
    
    pnet(udpCom.sock, 'setreadtimeout', 0);
    pnet(udpCom.sock, 'setwritetimeout', 0);
    disp('In FixStim');

    messageIsAvailable = 0;
    while ~messageIsAvailable  % endless loop.
        messageIsAvailable = pnet(udpCom.sock, 'readpacket', 1000, 'noblock');
        
        if messageIsAvailable
            DealWithMessage(messageIsAvailable);
            messageIsAvailable = 0;
        end
        
        if gl.dots.on
            DrawDots();
        end
        if gl.targ.on
            DrawTarg();
        end
        if gl.fp.on
            DrawFP();
        end
        if gl.fliprequest
            DoFlip();
        end

        if gl.timetoleave, return, end
        
        [keyisdown,junk,keycode] = KbCheck(); %#ok<ASGLU>
        if keyisdown && keycode(KEY.ESC)
            ShowCursor();
            pnet(udpCom.sock, 'close');
            Screen('CloseAll');
            return;
        end
    end
end

% DoFlip() taken from WhiteNoise
function DoFlip() % tidy this up to work with targ
    global gl udpCom
    
    Screen('Flip', gl.windowPtr, 0, 0);

    if(gl.targ.on == 1 && gl.framecounter == 0)
           pnet(udpCom.sock, 'write', 'MACSTIMON>> >>');
           pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
    end
    if (gl.framecounter == gl.framecountermax)
        pnet(udpCom.sock, 'write', 'MACSTIMOFF>> >>');
        pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
        gl.fliprequest = 0;
        gl.framecounter = 0;
        gl.framecountermax  = 0;
    end
    gl.framecounter = gl.framecounter + 1;
end

% DOES NOT STAND ALONE  Must inherit some parameters (e.g m fundamentals)
% from another function.
function InitDisplay(mondist, screenwidth, calfilename, fundfilename)
    global gl
    
    load(calfilename);
    cal = cals{end}; %#ok<USENS>
    gl.mondistcm = mondist;
    gl.screenWidthcm = screenwidth;
    gl.bkgndRGB = round(255*cal.bgColor)';
    gl.cal.gammaTable = cal.gammaTable;
    gl.cal.monSpd = cal.P_device;
    gl.cal.invgammaTable = InvertGammaTable(cal.gammaInput,cal.gammaTable,2^16);

    if (nargin > 3) % For backwards compatibility GDLH 2/17/13
        s = load(fundfilename);
        fns = fieldnames(s);
        P_device = SplineSpd(SToWls(cal.S_device), cal.P_device, SToWls([380 5 81]));
        gl.cal.fundamentals = eval(['s.',fns{1}])';
        gl.cal.M = gl.cal.fundamentals'*P_device;
        gl.cal.invM = inv(gl.cal.M);
    end
    
    if ~isempty(Screen('Windows'))
        gl.windowPtr = max(Screen('Windows'));
        Screen('FillRect', gl.windowPtr, cal.bgColor);
    else
        PsychImaging('PrepareConfiguration');
        if gl.vpixx
            PsychImaging('AddTask', 'General', 'EnableDataPixxC48Output', gl.ccmode); % in mode '1': every 2nd column of pixels is ignored
        else
            PsychImaging('AddTask', 'General', 'EnableBits++Color++Output', gl.ccmode);
        end
        gl.windowPtr = PsychImaging('OpenWindow', 0, cal.bgColor);
    end
    
    % Gamma correction is done in software so just set hardware gamma
    % lookup to the unity line.
    clut = repmat(linspace(0,1,256),3,1)';
    Screen('LoadNormalizedGammaTable', gl.windowPtr, clut);

    gl.framerate = Screen('FrameRate', gl.windowPtr, 1);
    [screenwidthpix, screenheightpix] = Screen('WindowSize', gl.windowPtr);
    gl.screenWidthpix = screenwidthpix;
    gl.screenHeightpix = screenheightpix;
    gl.screenCenterXpix = screenwidthpix/2;
    gl.screenCenterYpix = screenheightpix/2;

    pixpercm = gl.screenWidthpix/gl.screenWidthcm;
    theta = atan2(gl.screenWidthcm/2, gl.mondistcm)*180/pi;
    cmperdeg = gl.screenWidthcm/(2*theta);
    gl.pixperdeg = pixpercm*cmperdeg;

    gl.bkgndrgb = [gl.cal.gammaTable(gl.bkgndRGB(1)+1,1);...
        gl.cal.gammaTable(gl.bkgndRGB(2)+1,2);...
        gl.cal.gammaTable(gl.bkgndRGB(3)+1,3)];

    gl.fp.on = 0;
    gl.targ.on = 0;

    HideCursor();
end

% "Show" functions are called from REX.  They set up the appropriate fields
% in the "gl" structure (and set the "...on" toggle field to 1).
function ShowFP(x, y, fpsize, fpr, fpg, fpb)
    global gl

    gl.fp.x = x/10;
    gl.fp.y = y/10;
    gl.fp.size = fpsize/10;
    gl.fp.rgb = [fpr fpg fpb];
    
    xx = gl.fp.x*gl.pixperdeg;
    yy = gl.fp.y*gl.pixperdeg;

    fpsizeinpix = round(gl.pixperdeg*gl.fp.size);

    gl.fp.drawrect = [gl.screenCenterXpix + xx - floor(fpsizeinpix/2) + 1, ...
        gl.screenCenterYpix - yy - floor(fpsizeinpix/2) + 1, ...
        gl.screenCenterXpix + xx + ceil(fpsizeinpix/2)    , ...
        gl.screenCenterYpix - yy + ceil(fpsizeinpix/2)];

    gl.fp.on = 1;
end

% "Draw" functions get called on every screen refresh so long as the
% corresponding graphical object is to be displayed.  These are the
% functions that actually draw the object to the screen (and they make a
% fliprequest so that the drawn object will actually appear).
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

function ShowTarg(x, y, targsize, fpr, fpg, fpb)
    global gl

    gl.targ.x = x/10;
    gl.targ.y = y/10;
    gl.targ.size = targsize/10;
    gl.targ.rgb = [fpr fpg fpb];
    
    xx = gl.targ.x*gl.pixperdeg;
    yy = gl.targ.y*gl.pixperdeg;

    targsizeinpix = round(gl.pixperdeg*gl.targ.size);

    new_rect = [gl.screenCenterXpix + xx - floor(targsizeinpix/2) + 1
        gl.screenCenterYpix - yy - floor(targsizeinpix/2) + 1
        gl.screenCenterXpix + xx + ceil(targsizeinpix/2)
        gl.screenCenterYpix - yy + ceil(targsizeinpix/2)];

    gl.targ.drawrect = [gl.targ.drawrect new_rect];
    
    if size(gl.targ.drawrect,2) > 1
        % put the smaller (vis stim) in the first column
        targ_sizes = gl.targ.drawrect(3,:) - gl.targ.drawrect(1,:);
        [nil,sort_idx] = sort(targ_sizes);
        gl.targ.drawrect(:,sort_idx) = gl.targ.drawrect;
    end
    
    gl.targ.on = 1;
    gl.framecounter = 0;
end

function DrawTarg()
    global gl

    if ~isempty(gl.targ.drawrect)
        Screen('FillRect', gl.windowPtr, gl.targ.rgb, gl.targ.drawrect);
    end
    gl.fliprequest = 1;
end

function HideTarg()
    global gl

    if ~isempty(gl.targ.drawrect)
        % remove the first column which is the vis stim first
        % i.e., enforces that the vis stim is _always_ extinguished first
        gl.targ.drawrect(:,1) = [];
    else
        gl.targ.on = 0;
    end
    gl.fliprequest = 1;
end

function ShowLMScc(x,y,height,width,lcc,mcc,scc)
    global gl
    
    %[x y height width lcc mcc scc]
    
    stimcc = [lcc mcc scc];
    bkgndlms = gl.cal.M*gl.bkgndrgb;
    stimlms = bkgndlms.*(1+stimcc');
    gl.targ.rgb = gl.cal.invM*stimlms;
    if (any(gl.targ.rgb > 1) || any(gl.targ.rgb < 0))
        gl.targ.rgb = gl.bkgndrgb;
    end
    % Need to convert to DAC values
    % Not using the Bits++
    for plane = 1:3
        tmp = round(gl.targ.rgb(plane)*size(gl.cal.invgammaTable,1)-1)+1;
        tmp = gl.cal.invgammaTable(tmp, plane);
        gl.targ.rgb(plane) = round(tmp); 
    end
    
    gl.targ.x = x/10;
    gl.targ.y = y/10;
    gl.targ.size = [height width]./10;
    
    xx = gl.targ.x*gl.pixperdeg;
    yy = gl.targ.y*gl.pixperdeg;

    targheightinpix = round(gl.pixperdeg*gl.targ.size(1));
    targwidthinpix = round(gl.pixperdeg*gl.targ.size(2));

    gl.targ.drawrect = ...
       [gl.screenCenterXpix + xx - floor(targwidthinpix/2) + 1 ...
        gl.screenCenterYpix - yy - floor(targheightinpix/2) + 1 ...
        gl.screenCenterXpix + xx + ceil(targwidthinpix/2)     ...
        gl.screenCenterYpix - yy + ceil(targheightinpix/2)];

    gl.targ.on = 1;
end

function HideLMScc()
    global gl
    
    gl.targ.on = 0;
    gl.fliprequest = 1;
end

function ShowDots(trainingflag, x,y,speed,ndots,dir,apsize,dotsize,framecountermax, speedmult)
    global gl;

    if nargin < 10
        speedmult = 1;
    end
    if (trainingflag)
        napertures = 1;
    else
        napertures = 4;
    end
    apsize = apsize/10; % tenths of degrees
    speed = speed/10; % tenths of degrees/sec
    dotsize = dotsize/100; % hundredths of degrees
    x = x/10; % tenths of degrees
    y = y/10; % tenths of degrees
    
    gl.dots.apsizeinpix = round(gl.pixperdeg*apsize/2);  % /2 counting in double-width pixels
    gl.dots.ppf = speed * gl.pixperdeg / gl.framerate;  % dot speed in pixels per frame
    gl.dots.ppfdist = gl.dots.ppf*speedmult;
    gl.dots.dotsizeinpix = dotsize * gl.pixperdeg;    % Size of dots in pixels
%    gl.dots.positions = unidrnd(gl.dots.apsizeinpix*2,2,ndots,4);  % individual dot initial positions
    gl.dots.positions = unifrnd(-gl.dots.apsizeinpix, gl.dots.apsizeinpix,2,ndots,4);  % individual dot initial positions
    
    gl.dots.dir = dir;
    gl.framecountermax = framecountermax;
    
    % Sign flip on y values below to make the transformation
    % from DVA ('=' = up) to screen pixels ('+' = down)
    rotmat = [0 -1; 1 0]; % 90 degree rotation
    gl.dots.aperturepositions = [];
    for i = 1:napertures % rotating by 90 degrees each time. First is signal location.
        gl.dots.aperturepositions(i,:) = [gl.screenCenterXpix gl.screenCenterYpix]+gl.pixperdeg*[x -y]*rotmat^(i-1)';
    end
    gl.dots.on = 1;
end

function DrawDots()
    global gl
    
    dirsign = -(2*gl.dots.dir-1); % dirsign: '-1' means signal up, '1' means signal down
    % The way the drawing works, the top of the screen are small y
    % values.
try
    for j = 1:size(gl.dots.aperturepositions,1)
        L = sum(gl.dots.positions(:,:,j).^2) < gl.dots.apsizeinpix.^2;
        Screen('DrawDots',gl.windowPtr,gl.dots.positions(:,L,j),gl.dots.dotsizeinpix,[],gl.dots.aperturepositions(j,:))
    end
catch
    sca
    keyboard
end
    gl.dots.positions(:,:,1) = gl.dots.positions(:,:,1) + repmat([0;dirsign*gl.dots.ppf],...
        [1 size(gl.dots.positions,2) 1]); % first plane are signal dots
    gl.dots.positions(:,:,2:4) = gl.dots.positions(:,:,2:4) + repmat([0;-dirsign*gl.dots.ppfdist],...
        [1 size(gl.dots.positions,2) 3]);
    gl.dots.positions(gl.dots.positions(:) > gl.dots.apsizeinpix) = -gl.dots.apsizeinpix;
    gl.dots.positions(gl.dots.positions(:) < -gl.dots.apsizeinpix) = gl.dots.apsizeinpix;
    
    gl.fliprequest = 1;
end

function AllOff()
    global gl

    gl.targ.drawrect = [];
    gl.targ.on = 0;
    gl.dots.on = 0;
    gl.fp.on = 0;
    gl.fliprequest = 1;
    gl.framecounter = 0;
end

function DealWithMessage(msgSize)
    global udpCom gl
    message = pnet(udpCom.sock, 'read', msgSize, 'char');
    
    if strncmp(message, 'return', 6)
        a = dbstack;  % Check whether called from another function or from command line
        if ~strcmp(a(end).name, mfilename)
            gl.timetoleave = 1;
        end
    end
    
    try
        eval(message);
    catch ME
        fprintf('Ignoring uninterpretable message: "%s"\n', message);
        disp(getReport(ME));
    end
end