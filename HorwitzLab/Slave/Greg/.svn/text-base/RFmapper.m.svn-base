function RFmapper
    % RFmapper.m
    %
    % This will be a program that allows the user to map a RF by moving a bar of
    % light around on the screen.  It will be kind of like "play mode" in
    % Cortex.  It will present the fixation point and extinguish it in response
    % to commands from Rex.
    %
    % GDLH 3/6/07
    %
    % Porting over to work with pnet.
    %
    % GDLH 1/22/06
    %
    % Keyboard commands as of 1/25/08
    % Q: rotate counter clockwise
    % W: rotate clockwise
    % A: decrease width
    % S: increase width
    % Z: decrease height
    % X: increase height
    % 7 (on numeric pad): decrease red
    % 8 (on numeric pad): increase red
    % 4 (on numeric pad): decrease green
    % 5 (on numeric pad): increase green
    % 1 (on numeric pad): decrease blue
    % 2 (on numeric pad): increase blue
    %
    % If you click the button, the stimulus blinks for one frame and the (x, y)
    % coordinates are put into gl.bar.xy.  This gets inherited by WhiteNoise.m
    % (and all other Matlab functions because it's global) and sent to REX when
    % this information is requested.
    %
    % GDLH 2/4/08
    
    % Communication variables
    global udpCom
    udpCom.port = 6665;
    udpCom.sock = nan;
    udpCom.rexip = '192.168.1.120';
    
    % Keypress constants
    global KEY
    KbName('UnifyKeyNames');
    KEY.ESC = KbName('ESCAPE');
    KEY.Q = KbName('q');
    KEY.W = KbName('w');
    KEY.A = KbName('a');
    KEY.S = KbName('s');
    KEY.Z = KbName('z');
    KEY.X = KbName('x');
    KEY.NUMSEVEN = KbName('7');
    KEY.NUMEIGHT = KbName('8');
    KEY.NUMFOUR = KbName('4');
    KEY.NUMFIVE = KbName('5');
    KEY.NUMONE = KbName('1');
    KEY.NUMTWO = KbName('2');
    KEY.BLACKBKGND = KbName('-_');
    KEY.REGBKGND = KbName('=+');
    
    % Variables that are destructively modified by subfunctions
    global gl
    gl.mondist = 0;
    gl.screenwidth = 0;
    gl.pixperdeg = 0;
    gl.bkgndrgb = [0 0 0];
    gl.windowPtr = 0;
    gl.fliprequest = 0;
    gl.timetoleave = 0;
    
    gl.vpixx = ConnectedToViewPixx();
    gl.ccmode = 1;
    
    gl.fp.on = 0;
    gl.fp.x = 0;
    gl.fp.y = 0;
    gl.fp.size = 0;
    gl.fp.rgb = [0 0 0];
    gl.fp.drawrect = [];
    
    gl.bar.on = 1;
    gl.bar.blink = 0;
    gl.bar.height = 1;
    gl.bar.width = 1;
    gl.bar.rgb = [1 1 1];
    gl.bar.theta = 0;
    if (~isfield(gl.bar,'xy'))
        gl.bar.xy = [0 0];
    end
    
    [udpCom.sock, Success] = pnetStart(udpCom.port);
    if ~(Success)
        return
    end
    pnet(udpCom.sock, 'setreadtimeout', 0);
    pnet(udpCom.sock, 'setwritetimeout', 0);
    disp('In RFmapper');
    
    messageIsAvailable = 0;
    while ~messageIsAvailable  % endless loop.
        messageIsAvailable = pnet(udpCom.sock, 'readpacket', 1000, 'noblock');
        if (messageIsAvailable)
            DealWithMessage(messageIsAvailable);
            messageIsAvailable = 0;
        end
        if (gl.windowPtr > 0)
            [keyisdown,secs,keycodes] = KbCheck();
            if (keyisdown)
                keycode = find(keycodes);
                if keycode == KEY.ESC
                    pnet(udpCom.sock, 'close');
                    sca();
                    return
                else
                    InterpretKey(keycode);
                end
            end
            if (gl.bar.blink)
                Blinkbar();
            end
            if (gl.bar.on)
                Drawbar();
            end
            if (gl.fp.on)
                DrawFP();
            end
            if (gl.fliprequest)
                Screen('Flip',gl.windowPtr);
                gl.fliprequest = 0;
            end
        end
        if (gl.timetoleave)
            gl.stim.x = gl.bar.xy(1);
            gl.stim.y = gl.bar.xy(2);
            return
        end
    end
end

% InitDisplay below was taken verbatim from WhiteNoiseGH.
% Should these be a single function somewhere?

function InitDisplay(mondist, screenwidth, varargin) % varargin{1} and {2} are defined to be the calibration and fundamentals
    global gl
    
    calfilename = varargin{1};
    load(calfilename);
    cal = cals{end}; %#ok<USENS>
    gl.mondistcm = mondist;
    gl.screenWidthcm = screenwidth;
    gl.bkgndrgb = cal.bgColor';
    gl.cal.gammaTable = cal.gammaTable;
    gl.cal.monSpd = cal.P_device;
    gl.cal.invGamma = InvertGammaTable(cal.gammaInput, cal.gammaTable, 2^16);
    
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
            PsychImaging('AddTask', 'General', 'EnableDataPixxC48Output', gl.ccmode); % in mode '1': every 2nd column of pixels is ignored
        else
            PsychImaging('AddTask', 'General', 'EnableBits++Color++Output', gl.ccmode);
        end
        gl.windowPtr = PsychImaging('OpenWindow', 0, cal.bgColor);
    end
    
    gl.framerate = Screen('NominalFrameRate', gl.windowPtr, 1);
    [screenwidthpix, screenheightpix]  = Screen('WindowSize', gl.windowPtr);
    gl.screenWidthpix = screenwidthpix; % using Bits++ in Colour mode each pixel has a 1x2 aspect ratio
    gl.screenHeightpix = screenheightpix;  % but Psychophysicstoolbox doesn't (need to) know about this
    gl.screenCenterXpix = screenwidthpix/2;
    gl.screenCenterYpix = screenheightpix/2;
    
    pixpercm = gl.screenWidthpix/gl.screenWidthcm;
    theta = atan2(gl.screenWidthcm/2, gl.mondistcm)*180/pi;
    cmperdeg = gl.screenWidthcm/(2*theta);
    gl.pixperdeg = pixpercm*cmperdeg;
    
    gl.stim.on = 0;
    gl.fp.on = 0;
    
    HideCursor();
end

%%
function ShowFP(x, y, size, fpr, fpg, fpb)
    global gl
    
    gl.fp.x = x/10;
    gl.fp.y = y/10;
    gl.fp.size = size/10;
    gl.fp.rgb = [fpr fpg fpb]/255;
    
    fpsizeinpix = round(gl.pixperdeg*gl.fp.size);
    
    gl.fp.drawrect = [gl.screenCenterXpix+(gl.fp.x*gl.pixperdeg)-floor(fpsizeinpix/2)...
        gl.screenCenterYpix-(gl.fp.y*gl.pixperdeg)-floor(fpsizeinpix/2)...
        gl.screenCenterXpix+(gl.fp.x*gl.pixperdeg)+ceil(fpsizeinpix/2)...
        gl.screenCenterYpix-(gl.fp.y*gl.pixperdeg)+ceil(fpsizeinpix/2)];
    
    gl.fp.on = 1;
end

%%
function DrawFP()
    global gl
    
    Screen('FillOval', gl.windowPtr, gl.fp.rgb, gl.fp.drawrect);
    gl.fliprequest = 1;
end

%%
function HideFP()
    global gl
    
    gl.fp.on = 0;
end

%%
function DealWithMessage(msgSize)
    global udpCom gl
    
    message = pnet(udpCom.sock, 'read', msgSize, 'char');
    if (strncmp(message,'return',6))
        a = dbstack;  % Check whether called from another function or from command line
        if (~strcmp(a(end).name, mfilename))
            gl.timetoleave = 1;
        end
    end
    try
        eval(message);
    catch ME
        fprintf('Ignoring uninterpretable message: "%s"\n',message);
        disp(getReport(ME));
    end
end

%%
function Drawbar()
    global gl
    
    barheightinpix = round(gl.bar.height*gl.pixperdeg);
    barwidthinpix = round(gl.bar.width*gl.pixperdeg);
    
    [mouseposX, mouseposY, buttons] = GetMouse();
    if (buttons(1))
        gl.bar.xy = [(mouseposX-gl.screenCenterXpix)/gl.pixperdeg ...
            -(mouseposY-gl.screenCenterYpix)/gl.pixperdeg];  % y convention is opposite for mouse and screen
        gl.bar.blink = 1;
    end
    
    drawrect = [mouseposX-floor(barwidthinpix/2)...
        mouseposY-floor(barheightinpix/2)...
        mouseposX+ceil(barwidthinpix/2)...
        mouseposY+ceil(barheightinpix/2)];
    
    img = bsxfun(@plus, zeros(barwidthinpix,barheightinpix,3), reshape(gl.bar.rgb, [1 1 3]));
    tex = Screen('MakeTexture', gl.windowPtr, img, [], [], 2);
    
    Screen ('DrawTexture', gl.windowPtr, tex, [], drawrect, round(gl.bar.theta), 1);
    Screen('Close', tex);
    
    gl.fliprequest = 1;
end

%%
function Blinkbar()
    global gl
    
    if (gl.bar.on)
        gl.bar.on = 0;
    else
        gl.bar.on = 1;
        gl.bar.blink = 0;
    end
    gl.fliprequest = 1;
end

%%
% function SendRFPos(x,y)
%     global gl
%     
%     w = gl.windowPtr;
%     barheightinpix = round(gl.bar.height*gl.pixperdeg);
%     barwidthinpix = round(gl.bar.width*gl.pixperdeg);
%     [screenWidth, screenHeight] = Screen('WindowSize',w);
%     
%     % x and y to the nearest 10th of a degree
%     x = round((x-screenWidth/2)*10/gl.pixperdeg);
%     y = round((y-screenHeight/2)*10/gl.pixperdeg);
%     messageIsAvailable = 0;
%     
%     % matlabUDP('send', sprintf('RF(%d,%d)',x,-y));
%     % Reverse the sign of y because of differences in PTB and REX
%     % conventions.
%     
% end

%%
function InterpretKey(keycode)
    if (length(keycode) > 1)
        keycode = keycode(1);
    end
    
    UpdateBar(keycode);
    UpdateBkgnd(keycode);
end

%%
function UpdateBkgnd(keycode)
    global gl KEY
    
    if keycode == KEY.BLACKBKGND
        Screen('FillRect', gl.windowPtr, [0 0 0]);
    elseif keycode == KEY.REGBKGND
        Screen('FillRect', gl.windowPtr, gl.bkgndrgb);
    end
end

%%
function UpdateBar(keycode)
    global gl KEY
    
    % gives about a 2% change in gun intensity per step
    if gl.vpixx % 10-bit color depth
        GUNSTEP = 5 * 2^2 / (2^10-1);
    else % bits++: 14-bit color depth
        GUNSTEP = 5 * 2^6 / (2^14-1);
    end
    
    barheightinpix = gl.bar.height*gl.pixperdeg;
    barwidthinpix = gl.bar.width*gl.pixperdeg;
    
    switch (keycode)
        case KEY.Q
            gl.bar.theta = mod(gl.bar.theta-1,180);
        case KEY.W
            gl.bar.theta = mod(gl.bar.theta+1,180);
        case KEY.A
            gl.bar.width = max([1 barwidthinpix-1])/gl.pixperdeg;
        case KEY.S
            gl.bar.width = min([gl.screenWidthpix barwidthinpix+1])/gl.pixperdeg;
        case KEY.Z
            gl.bar.height = max([1 barheightinpix-1])/gl.pixperdeg;
        case KEY.X
            gl.bar.height = min([gl.screenHeightpix barheightinpix+1])/gl.pixperdeg;
        case KEY.NUMSEVEN
            gl.bar.rgb(1) = max([0 gl.bar.rgb(1)-GUNSTEP]);
        case KEY.NUMEIGHT
            gl.bar.rgb(1) = min([1 gl.bar.rgb(1)+GUNSTEP]);
        case KEY.NUMFOUR
            gl.bar.rgb(2) = max([0 gl.bar.rgb(2)-GUNSTEP]);
        case KEY.NUMFIVE
            gl.bar.rgb(2) = min([1 gl.bar.rgb(2)+GUNSTEP]);
        case KEY.NUMONE
            gl.bar.rgb(3) = max([0 gl.bar.rgb(3)-GUNSTEP]);
        case KEY.NUMTWO
            gl.bar.rgb(3) = min([1 gl.bar.rgb(3)+GUNSTEP]);
        case KEY.BLACKBKGND
            gl.bar.rgb(1:3) = 0;
    end
end
