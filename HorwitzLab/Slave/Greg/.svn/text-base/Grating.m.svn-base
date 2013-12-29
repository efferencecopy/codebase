function Grating
% Grating.m
%
%    Slave problem for displaying drifing gratings.
%
% GDLH 3/1/08

    % Communication variables
    global udpCom;
    udpCom.port = 6665;
    udpCom.sock = nan;
    udpCom.rexip = '192.168.1.120';

    % Keypress constants
    global KEY
    KbName('UnifyKeyNames');
    KEY.ESC = KbName('escape');
    
    % Variables that are destructively modified by subfunctions
    global gl;
    gl.mondistcm = 0;
    gl.screenWidthcm = 0;
    gl.screenWidthpix = 0;
    gl.screenHeightpix = 0;
    gl.screenCenterXpix = 0;
    gl.screenCenterYpix = 0;
    gl.pixperdeg = 0;
    gl.bkgndRGB = [0 0 0];
    gl.bkgndrgb = [0 0 0];
    gl.windowPtr = 0;
    gl.fliprequest = 0;
    gl.framerate = 0;
    gl.framecounter = 0;
    gl.framecountermax = 0;
    gl.timetoleave = 0;
    gl.vpixx = ConnectedToViewPixx();
    gl.ccmode = 1;
    
    gl.ep.x = 0;
    gl.ep.y = 0;
    
    gl.fp.on = 0;
    gl.fp.x = 0;
    gl.fp.y = 0;
    gl.fp.size = 0;
    gl.fp.rgb = [0 0 0];
    
    gl.grating.on = 0;
    gl.grating.drawrect = [0 0 0 0];
    gl.grating.diam = 1;
    gl.grating.sf = 1;
    gl.grating.tf = 1;
    gl.grating.conecontrast = [0 0 0];
    gl.grating.rgb = [0 0 0];
    gl.grating.phase = 0;
    gl.grating.orient = 0;
    gl.grating.aperture = [];
    
    gl.cal.gammaTable = [];
    gl.cal.invgammaTable = [];
    gl.cal.monSpd = [];
    gl.cal.fundamentals = [];
    gl.cal.M = zeros(3);
    gl.cal.invM = zeros(3);

    [udpCom.sock, Success] = pnetStart(udpCom.port);
    if ~(Success)
        return
    end
    pnet(udpCom.sock, 'setreadtimeout', 0);
    pnet(udpCom.sock, 'setwritetimeout', 0);
    disp('In Grating');

    messageIsAvailable = 0;
    while ~messageIsAvailable  % endless loop.
        messageIsAvailable = pnet(udpCom.sock, 'readpacket', 1000, 'no block');
        if (messageIsAvailable)
            DealWithMessage(messageIsAvailable);
            messageIsAvailable = 0;   
        end     
        if (gl.windowPtr > 0)
            [keyisdown,secs,keycode] = KbCheck();
            if (gl.grating.on)
                Drawgrating();
            end
            if (gl.fp.on)
                DrawFP();
            end
            if (gl.fliprequest)
                DoFlip();
            end
        end
        if (gl.timetoleave)
            return;
        end
        [keyisdown,secs,keycode] = KbCheck();
        if (keyisdown && keycode(KEY.ESC))
           ShowCursor;
           pnet(udpCom.sock, 'close');
           Screen('CloseAll');
           return;
        end
    end
end

%%
% DoFlip() taken from WhiteNoise
% Should these be a single function somewhere? 
function DoFlip()
    global gl;
    global udpCom;
    
    Screen('Flip',gl.windowPtr,0,0);
    if (gl.grating.on)
        if(gl.framecounter == 0)
           pnet(udpCom.sock, 'write', 'MACSTIMON>> >>');
           pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
        end
        gl.framecounter = gl.framecounter + 1;
    end
    gl.fliprequest = 0;    
end

%%
% InitDisplay below was taken verbatim from WhiteNoiseGH.
% Should these be a single function somewhere?

function InitDisplay(mondist, screenwidth, calfilename, fundfilename)
    global gl;
    
    load(calfilename);
    cal = cals{end};
    gl.mondistcm = mondist;
    gl.screenWidthcm = screenwidth;
    gl.bkgndRGB = round(255*cal.bgColor)';
    gl.cal.gammaTable = cal.gammaTable;
    gl.cal.invgammaTable = InvertGammaTable(cal.gammaInput, gl.cal.gammaTable, 2^16);
    gl.cal.monSpd = cal.P_device;
    gl.stim.on = 0;
    gl.fp.on = 0;

    s = load(fundfilename);
    fns = fieldnames(s);
    P_device = SplineSpd(SToWls(cal.S_device), cal.P_device, SToWls([380 5 81]));
    gl.cal.fundamentals = eval(['s.',fns{1}])';
    gl.cal.M = gl.cal.fundamentals'*P_device;
    gl.cal.invM = inv(gl.cal.M);
    
    gl.bkgndrgb = [gl.cal.gammaTable(gl.bkgndRGB(1)+1,1),...
                   gl.cal.gammaTable(gl.bkgndRGB(2)+1,2),...
                   gl.cal.gammaTable(gl.bkgndRGB(3)+1,3)]';

    
    % Gamma correction is done in software so just set hardware gamma
    % lookup to the unity line.
    clut = repmat(linspace(0,1,256),3,1)';
    Screen('LoadNormalizedGammaTable', gl.windowPtr, clut);
        
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
    gl.framerate = Screen('FrameRate', gl.windowPtr, 1);
    [screenwidthpix, screenheightpix]  = Screen('WindowSize', gl.windowPtr);
    gl.screenWidthpix = screenwidthpix; % using Bits++ in Colour mode each pixel has a 1x2 aspect ratio
    gl.screenHeightpix = screenheightpix;  % but Psychophysicstoolbox doesn't (need to) know about this
    gl.screenCenterXpix = screenwidthpix/2;
    gl.screenCenterYpix = screenheightpix/2;
    
    pixpercm = gl.screenWidthpix/gl.screenWidthcm;
    theta = atan2(gl.screenWidthcm/2, gl.mondistcm)*180/pi;
    cmperdeg = gl.screenWidthcm/(2*theta);
    gl.pixperdeg = pixpercm*cmperdeg;
    
    HideCursor;
end

%%
% "Show" functions are called from REX.  They set up the appropriate fields 
% in the "gl" structure (and set the "...on" toggle field to 1). 
function ShowFP(x, y, size, fpr, fpg, fpb)
    global gl;

    gl.fp.x = x/10;
    gl.fp.y = y/10;
    gl.fp.size = size/10;
    gl.fp.rgb = [fpr, fpg, fpb]/255;
    
    fpsizeinpix = round(gl.pixperdeg*gl.fp.size);

    gl.fp.drawrect = [gl.screenCenterXpix+(gl.fp.x*gl.pixperdeg)-floor(fpsizeinpix/2)+1 ...
                gl.screenCenterYpix-(gl.fp.y*gl.pixperdeg)-floor(fpsizeinpix/2)+1 ...
                gl.screenCenterXpix+(gl.fp.x*gl.pixperdeg)+ceil(fpsizeinpix/2)...
                gl.screenCenterYpix-(gl.fp.y*gl.pixperdeg)+ceil(fpsizeinpix/2)];
    
    gl.fp.on = 1;
    
end

%%
% "Draw" functions get called on every screen refresh so long as the
% corresponding graphical object is to be displayed.  These are the
% functions that actually draw the object to the screen (and they make a
% fliprequest so that the drawn object will actually appear).
function DrawFP()
    global gl;

    Screen('Fillrect', gl.windowPtr, gl.fp.rgb, gl.fp.drawrect);
    
    gl.fliprequest = 1;
end

%%
function HideFP()
    global gl;
   
    gl.fp.on = 0;
    gl.fliprequest = 1;

end

%%
function ShowGrating(rfx, rfy, diam, sf, tf,lcont, mcont, scont,  orient, phase, ncycles)

    global gl;
    
    gl.grating.x = rfx/10;
    gl.grating.y = rfy/10;
    gl.grating.diam = diam;
    gl.grating.sf = sf;
    gl.grating.tf = tf;
    gl.grating.conecontrast = [lcont mcont scont];
    gl.grating.orient = orient;
    gl.grating.phase = phase;
    if (tf == 0)  % If TF is zero, intrepret ncycles as number of frames
        gl.framecountermax = ncycles;
    else
        gl.framecountermax = floor(ncycles/tf*gl.framerate);
    end
    
    if (gl.grating.diam == 0)  % Why would we ever need this?
        return;
    end
    
    % Calculating the drawing rectangle
    stimsizeinpix = round(gl.pixperdeg*gl.grating.diam/2);  % /2 counting in double-width pixels 
    x = (gl.grating.x+gl.ep.x)*gl.pixperdeg;
    y = (gl.grating.y+gl.ep.y)*gl.pixperdeg;
    gl.grating.drawrect = round([gl.screenCenterXpix+x-stimsizeinpix ...
                gl.screenCenterYpix-y-stimsizeinpix ...
                gl.screenCenterXpix+x+stimsizeinpix ...
                gl.screenCenterYpix-y+stimsizeinpix]);
    if(rem(gl.grating.drawrect(1), 2)) %if the rectangle starts on an odd pixel
        gl.grating.drawrect(1) = gl.grating.drawrect(1) - 1;
        gl.grating.drawrect(3) = gl.grating.drawrect(3) - 1;
    end
    gl.ep.x = 0; % resetting
    gl.ep.y = 0;
    
    % Calculating grating rgbs     
    bkgndlms = gl.cal.M*gl.bkgndrgb;
    gratinglms = bkgndlms.*[1+gl.grating.conecontrast'];
    gl.grating.rgb = gl.cal.invM*gratinglms;
    if (any(gl.grating.rgb > 1) | any(gl.grating.rgb < 0))
        gl.grating.rgb = gl.bkgndrgb;
    end
    
    % Creating an aperture template
    [x,y] = meshgrid(linspace(-1,1,stimsizeinpix), linspace(-1,1,2*stimsizeinpix));
    % Only half the number of columns as rows since we're using colour mode
    gl.grating.aperture = (x.^2+y.^2) <= 1+2/stimsizeinpix;  % correction for binning
    
    gl.framecounter = 0;
    gl.grating.on = 1;

end
%%
function Drawgrating()
    global gl;
    global udpCom
    
    w = gl.windowPtr;
    
    pixpercycle = gl.pixperdeg/gl.grating.sf;
    sizeinpix = round(gl.pixperdeg*gl.grating.diam/2)*2;
    
    xinc = (2*pi/pixpercycle)*cos((pi/2)-gl.grating.orient);
    yinc = (2*pi/pixpercycle)*sin((pi/2)-gl.grating.orient);
    [xramp, yramp] = meshgrid(xinc*([0:2:sizeinpix-1]), yinc*([0:sizeinpix-1]));
    a = cos(xramp+yramp+gl.grating.phase);
    a = a.*gl.grating.aperture;
    im = zeros(sizeinpix, sizeinpix/2, 3);
    for plane = 1:3
        tmp = a.*(gl.grating.rgb(plane)-gl.bkgndrgb(plane))+gl.bkgndrgb(plane);
        tmp = round(tmp*size(gl.cal.invgammaTable,1)-1)+1;
        tmp = gl.cal.invgammaTable(tmp, plane);
       %tmp = round(tmp*(size(gl.cal.invgammaTable,1)-1));
        im(:,:,plane) = reshape(tmp, sizeinpix, sizeinpix/2);
    end
    
    if (gl.framecounter == gl.framecountermax)
        pnet(udpCom.sock, 'write', 'MACSTIMOFF>> >>');
        pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
        gl.grating.on = 0;
    end
    phaseinc = gl.grating.tf*2*pi/gl.framerate;
    gl.grating.phase = mod(gl.grating.phase+phaseinc, 2*pi);
  
    tex=Screen('MakeTexture', gl.windowPtr, TranslateToColourMode(im,[],gl.ccmode),[],[],2);
    Screen('DrawTexture',gl.windowPtr,tex, [], gl.grating.drawrect,[],0);
    Screen('Close',tex);
    
    gl.fliprequest = 1;
end

%%  
function AllOff()
    global gl;
   
    gl.grating.on = 0;
    gl.fp.on = 0;
    gl.fliprequest = 1;
    
end

%%
function eyepos(x, y)
    global gl;

    gl.ep.x = x/40;
    gl.ep.y = y/40;
end

%%
function DealWithMessage(msgSize)
    global udpCom;
    global gl;  % needs to be here for functions evaled by this one
    message = pnet(udpCom.sock, 'read', msgSize, 'char');
    if (strncmp(message,'return',6))
        a = dbstack;  % Check whether called from another function or from command line 
        if (~strcmp(a(end).name, mfilename))
            gl.timetoleave = 1;
        end
    end
    try 
        eval(message)
    catch
        fprintf('Ignoring uninterpretable message: "%s"\n',message);
        error = lasterror;
        disp(error.message);
        disp(error.identifier);
        disp(error.stack);
    end
end


