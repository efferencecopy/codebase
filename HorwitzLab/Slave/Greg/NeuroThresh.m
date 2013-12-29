function NeuroThresh
% NeuroThresh.m
%
%    Slave program for displaying Gabor patches.  This will be used
% for adaptively finding the neuron's "threshold".
%
% GDLH 6/9/09


    % Communication variables
    global udpCom;
    udpCom.port = 6665;
    udpCom.sock = nan;
    udpCom.rexip = '192.168.1.120';

    % Keypress constants
    global KEY;
    KbName('UnifyKeyNames');
    KEY.ESC = KbName('escape');
    
    global gl;
   % gl.mondistcm = 0;
   % gl.screenWidthcm = 0;
   % gl.screenWidthpix = 0;
   % gl.screenHeightpix = 0;
   % gl.screenCenterXpix = 0;
   % gl.screenCenterYpix = 0;
  %  gl.pixperdeg = 0;
  %  gl.bkgndrgb = [0 0 0];
   % gl.windowPtr = 0;
   % gl.cal.gammaTable = []; 
   % gl.cal.monSpd = [];
   % gl.cal.fundamentals = [];
   % gl.cal.M = zeros(3);
   % gl.cal.invM = zeros(3);
    
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
    
    gl.vpixx = ConnectedToViewPixx();
    gl.ccmode = 1;
   
    gl.stim.on = 0;
    gl.stim.template = [];
    
    [udpCom.sock, Success] = pnetStart(udpCom.port);
    if ~(Success)
        return
    end
    pnet(udpCom.sock, 'setreadtimeout', 0);
    pnet(udpCom.sock, 'setwritetimeout', 0);
    disp('In NeuroThresh');

    messageIsAvailable = 0;
    while ~messageIsAvailable  % endless loop.
        messageIsAvailable = pnet(udpCom.sock, 'readpacket', 1000, 'no block');
        if (messageIsAvailable)
            DealWithMessage(messageIsAvailable);
            messageIsAvailable = 0;   
        end     
        if (gl.windowPtr > 0)
            if (gl.stim.on)
                DrawStim();
            end
            if (gl.fp.on)
                DrawFP();
            end
            if (gl.fliprequest)
                DoFlip();
            end
        end
        if (gl.timetoleave)
            gl.stim.template = [];
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
function DoFlip()
    global gl;
    global udpCom;
    
    Screen('Flip',gl.windowPtr,0,0);
    if (gl.stim.on)
        if(gl.framecounter == 0)
           pnet(udpCom.sock, 'write', 'MACSTIMON>> >>');
           pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
        end
        if (gl.framecounter == gl.framecountermax-1)
            pnet(udpCom.sock, 'write', 'MACSTIMOFF>> >>');
            pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
            gl.stim.on = 0;
        end
        gl.framecounter = gl.framecounter + 1;
    end
    gl.fliprequest = 0;    
end

%%

function InitDisplay(mondist, screenwidth, calfilename, fundfilename)

    global gl;
    
    load(calfilename);
    cal = cals{end};
    gl.mondistcm = mondist;
    gl.screenWidthcm = screenwidth;
    gl.bkgndRGB = round(255*cal.bgColor)';
    gl.cal.gammaTable = cal.gammaTable;
    gl.cal.invgammaTable = InvertGammaTable(cal.gammaInput, ...
        gl.cal.gammaTable, 2^16);
    gl.cal.monSpd = cal.P_device;
   
    s = load(fundfilename);
    fns = fieldnames(s);
    P_device = SplineSpd(SToWls(cal.S_device), cal.P_device, SToWls([380 5 81]));
    gl.cal.fundamentals = eval(['s.',fns{1}])';
    gl.cal.M = gl.cal.fundamentals'*P_device;
    gl.cal.invM = inv(gl.cal.M);
    
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
    gl.screenCenterXpix = screenwidthpix/2;
    gl.screenCenterYpix = screenheightpix/2;
    
    pixpercm = gl.screenWidthpix/gl.screenWidthcm;
    theta = atan2(gl.screenWidthcm/2, gl.mondistcm)*180/pi;
    cmperdeg = gl.screenWidthcm/(2*theta);
    gl.pixperdeg = pixpercm*cmperdeg;

    gl.bkgndrgb = [gl.cal.gammaTable(gl.bkgndRGB(1)+1,1);...
                  gl.cal.gammaTable(gl.bkgndRGB(2)+1,2);...
                  gl.cal.gammaTable(gl.bkgndRGB(3)+1,3)];  
              
  
    gl.stim.on = 0;
    gl.fp.on = 0;

    HideCursor;

end
%%
% Precomputing the space/time gabor so that later on all we have to do is
% color it.  Note sure how much time this is really saving us.  In the
% future, might want to combine PrepareGabor and ShowStim into a single
% function that is called immediately before each stimulus presentation
% (if we want to change the space/time parameters of the Gabor between 
% stimulus presentations).
function PrepareGabor(nframesplateau, nframesramp, theta, sf, phi, sigma, gamma, tf, nsigmas)
    global gl;
 
    lambda = 1/sf;

    ramp = linspace(0,1,nframesramp);
    plateau = ones(1, nframesplateau);
    temporalprofile = [ramp plateau fliplr(ramp)];
    nframes = length(temporalprofile);
    stimsizeindeg = sigma*nsigmas; % This is half stim size (we we go out +/- nsigmas)
    stimsizeinpix = round(stimsizeindeg*gl.pixperdeg);  % full stim size in doublewide pixels
    [x,y] = meshgrid(stimsizeindeg*linspace(-1,1,stimsizeinpix), stimsizeindeg*linspace(-1,1,2*stimsizeinpix));
    % x and y are in dva

    xprime = x*cos(-theta) + y*sin(-theta);
    yprime = -x*sin(-theta) + y*cos(-theta);
    
    deltaphase = tf*2*pi/gl.framerate;
    gabor = zeros(2*stimsizeinpix,stimsizeinpix,nframes);
    for a = 1:length(temporalprofile)
        phase = phi+(a-1)*deltaphase;
        gabor(:,:,a) = temporalprofile(a)*exp(-(xprime.^2 + gamma.^2 .* yprime.^2)./ (2*sigma^2)).*cos(2*pi*yprime./lambda+phase); 
    end
    gl.stim.template = gabor;
end

%%
% Making a half-wave rectified Gabor
function PrepareHalfRectGabor(nframesplateau, nframesramp, theta, sf, phi, sigma, gamma, tf, nsigmas)
    global gl;
       
    PrepareGabor(nframesplateau, nframesramp, theta, sf, phi, sigma, gamma, tf, nsigmas)

    gl.stim.template(gl.stim.template < 0) = 0;
end
%%
% "Show" functions are called from REX.  They set up the appropriate fields 
% in the "gl" structure (and set the "...on" toggle field to 1). 
function ShowFP(x, y, size, fpr, fpg, fpb)
    global gl;

    gl.fp.x = x/10;
    gl.fp.y = y/10;
    gl.fp.size = size/10;
    gl.fp.rgb = [fpr, fpg, fpb];
    
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
% 
function ShowStim(stimx, stimy, l, m, s)

    global gl;
    
    NGAMMASTEPS = size(gl.cal.invgammaTable,1);

    stimx = stimx/10;
    stimy = stimy/10;

    stimconecontrast = [l m s];
    
    % Creating the drawing window
    x = stimx*gl.pixperdeg;
    y = stimy*gl.pixperdeg;
    stimsizeinpix = size(gl.stim.template,2); % in doublewide pix
    gl.drawrect = round([gl.screenCenterXpix+x-stimsizeinpix ...
                gl.screenCenterYpix-y-stimsizeinpix ...
                gl.screenCenterXpix+x+stimsizeinpix...
                gl.screenCenterYpix-y+stimsizeinpix]);

            
    if(rem(gl.drawrect(1), 2)) %if the rectangle starts on an odd pixel
        gl.drawrect(1) = gl.drawrect(1) - 1;
        gl.drawrect(3) = gl.drawrect(3) - 1;
    end
    
    if (gl.drawrect(3)-gl.drawrect(1) ~= 2*stimsizeinpix || gl.drawrect(4)-gl.drawrect(2) ~= 2*stimsizeinpix)
        sca;
        keyboard;
        % Draw window is the wrong size
    end

    % Getting RGBs
    bkgndlms = gl.cal.M*gl.bkgndrgb;
    gl.stim.rgb = gl.cal.invM*(bkgndlms.*(1+stimconecontrast'));

   
    % Checking for out of gamut errors
    % Right now just squeezing it back into the gamut - need to send a
    % message to REX that this was done.
    lims = min([gl.bkgndrgb'; 1-gl.bkgndrgb']);
    if (any(abs(gl.stim.rgb-gl.bkgndrgb) > lims'))
        scalefactor = max(abs(gl.stim.rgb-gl.bkgndrgb)./lims');  % +(1/NGAMMASTEPS);
        gl.stim.rgb = (gl.stim.rgb-gl.bkgndrgb)./scalefactor+gl.bkgndrgb;
    end
    % Precomputing all of the textures.
    % Not sure that this is necessary or even a good idea.
%     for a = 1:size(gl.stim.template,3)  % looping over frames
%         im = zeros(stimsizeinpix*2, stimsizeinpix, 3);
%         frame = gl.stim.template(:,:,a);
%         for plane = 1:3
%             tmp = frame.*(gl.stim.rgb(plane)-gl.bkgndrgb(plane)) + gl.bkgndrgb(plane);
%             tmp = round(tmp*(NGAMMASTEPS-1))+1;
%             tmp = gl.cal.invgammaTable(tmp, plane);
%             tmp = round(tmp*(NGAMMASTEPS-1));
%             im(:,:,plane) = reshape(tmp, stimsizeinpix*2, stimsizeinpix);
%         end
%         gl.stim.tex(a)=Screen('MakeTexture', gl.windowPtr, TranslateToColourModeMex(im));
%     end
    
    gl.framecounter = 0;
    gl.framecountermax = size(gl.stim.template,3);
    gl.stim.on = 1;
 end
%%
function DrawStim()

    global gl;

   % Line below is for use when we're precomputing the textures
   % Screen('DrawTexture',gl.windowPtr,gl.stim.tex(gl.framecounter+1), [], gl.drawrect,[],0);

   % The code below is for use when we're computing textures on a frame by
   % frame basis
   NGAMMASTEPS = size(gl.cal.invgammaTable,1);
   im = zeros(size(gl.stim.template,1), size(gl.stim.template,2), 3);
   frame = gl.stim.template(:,:,gl.framecounter+1);
   for plane = 1:3
       tmp = frame.*(gl.stim.rgb(plane)-gl.bkgndrgb(plane)) + gl.bkgndrgb(plane);
       tmp = round(tmp*(NGAMMASTEPS-1))+1;
       try
        tmp = gl.cal.invgammaTable(tmp, plane);
       catch
           sca
           keyboard
       end
       im(:,:,plane) = reshape(tmp, size(im(:,:,1)));
   end
   tex=Screen('MakeTexture', gl.windowPtr, TranslateToColourMode(im, [], gl.ccmode), [], [], 2);
   Screen('DrawTexture',gl.windowPtr,tex, [], gl.drawrect,[],0);
   Screen('Close',tex);
   gl.fliprequest = 1;
end

%%
function HideStim()
    global gl;

    gl.stim.on = 0;
    gl.fliprequest = 1;
     
end
%%  
function AllOff()
    global gl;
   
    gl.stim.on = 0;
    gl.fp.on = 0;
    gl.fliprequest = 1;
    
end

%%
function DealWithMessage(msgSize)
    global udpCom;
    global gl;  % needs to be here for functions evaled by this one
    message = pnet(udpCom.sock, 'read', msgSize, 'char')
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


