%#ok<*DEFNU>
function memsacc()
    % Communication variables
    global udpCom KEY gl
    udpCom.port = 6665;
    udpCom.sock = nan;
    udpCom.rexip = '192.168.1.120';

    % Keypress constants
    KEY.ESC = 41;

    % Variables that are destructively modified by subfunctions
    % Because we are guaranteed to be running this code after having
    % already
    % setup the global gl structure (via the white noise paradigm)
    % there's no sense in reassigning everything.
    gl.usecolourmode = 1;  % Are we using the Bits++

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
    gl.stim.x = 0;
    gl.stim.y = 0;
    gl.stim.size = 0;
    gl.stim.rgb = [0 0 0];
    gl.stim.drawrect = [0 0 0 0];

    [udpCom.sock, Success] = pnetStart(udpCom.port);
    if ~Success, return, end
    
    pnet(udpCom.sock, 'setreadtimeout', 0);
    pnet(udpCom.sock, 'setwritetimeout', 0);
    disp('In memsacc');

    messageIsAvailable = 0;
    while ~messageIsAvailable  % endless loop.
        messageIsAvailable = pnet(udpCom.sock, 'readpacket', 1000, 'noblock');
        
        if messageIsAvailable
            DealWithMessage(messageIsAvailable);
            messageIsAvailable = 0;
        end

        if gl.fp.on
            DrawFP();
        end
        if gl.stim.on
            DrawStim();
        end
        if gl.fliprequest
            DoFlip();
        end

        if gl.timetoleave, return, end
        
        [keyisdown,junk,keycode] = KbCheck(); %#ok<ASGLU>
        if keyisdown && keycode(KEY.ESC)
            ShowCursor();
            pnet(udpCom.sock, 'close');
            sca();
            return;
        end
    end
end

function DoFlip()
    global gl
    
    Screen('Flip', gl.windowPtr, 0, 0);
    gl.fliprequest = 0;
end

% DOES NOT STAND ALONE  Must inherit some parameters (e.g m fundamentals)
% from another function.
function InitDisplay(mondist, screenwidth, calfilename)
    global gl
    
    load(calfilename);
    cal = cals{end}; %#ok<USENS>
    gl.mondistcm = mondist;
    gl.screenWidthcm = screenwidth;
    gl.bkgndRGB = round(255*cal.bgColor)';
    gl.cal.gammaTable = cal.gammaTable;
    gl.cal.monSpd = cal.P_device;
    gl.cal.invgammaTable = InvertGamma(gl.cal.gammaTable, 1);

    if isempty(Screen('windows'))
        gl.windowPtr = Screen('OpenWindow', 0, gl.bkgndRGB);
    else
        gl.windowPtr = Screen('windows');
        gl.windowPtr = gl.windowPtr(1);
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
function ShowFP(x, y, size, fpr, fpg, fpb)
    global gl

    gl.fp.x = x/10;
    gl.fp.y = y/10;
    gl.fp.size = size/10;
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

function ShowStim(x, y, size, stimr, stimg, stimb)
    global gl
    
    gl.stim.x = x/10;
    gl.stim.y = y/10;
    gl.stim.size = size/10;
    gl.stim.rgb = [stimr stimg stimb];
    
    xx = gl.stim.x*gl.pixperdeg;
    yy = gl.stim.y*gl.pixperdeg;

    fpsizeinpix = round(gl.pixperdeg*gl.stim.size);

    gl.stim.drawrect = [gl.screenCenterXpix + xx - floor(fpsizeinpix/2) + 1, ...
        gl.screenCenterYpix - yy - floor(fpsizeinpix/2) + 1, ...
        gl.screenCenterXpix + xx + ceil(fpsizeinpix/2)    , ...
        gl.screenCenterYpix - yy + ceil(fpsizeinpix/2)];

    gl.stim.on = 1;
end

% "Draw" functions get called on every screen refresh so long as the
% corresponding graphical object is to be displayed.  These are the
% functions that actually draw the object to the screen (and they make a
% fliprequest so that the drawn object will actually appear).
function DrawFP()
    global gl

    Screen('Fillrect', gl.windowPtr, gl.fp.rgb, gl.fp.drawrect);
    gl.fliprequest = 1;
end

function DrawStim()
    global gl

    Screen('Fillrect', gl.windowPtr, gl.stim.rgb, gl.stim.drawrect);
    gl.fliprequest = 1;
end

function HideFP()
    global gl

    gl.fp.on = 0;
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
        stk = dbstack();  % Check whether called from another function or from command line
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