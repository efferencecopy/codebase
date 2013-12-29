function SMurray
% SMurray.m
%
%    Slave program for displaying rings of various widths.
%
% GDLH 4/18/09
%
%   Also displays Gabors (flankers + target in RF)
%
% GDLH 12/7/10
%
% Experimenting with bitmaps in the background
% GDLH 11/7/11

    % Communication variables
    global udpCom;
    udpCom.port = 6665;
    udpCom.sock = nan;
    udpCom.rexip = '192.168.1.120';

    % Keypress constants
    global KEY;
    KEY.ESC = 41;
    
    % Variables that are destructively modified by subfunctions
    % Because we are guaranteed to be running this code after having
    % already 
    % setup the global gl structure (via the white noise paradigm)
    % there's no sense in reassigning everything.
    global gl;
    gl.usecolourmode = 0;  % Are we using the Bits++
    gl.usebitmapbkgnd = 1; % For debugging
%    gl.mondistcm = 0;       % We can rely on inheriting this stuff
%    gl.screenWidthcm = 0;
%    gl.screenWidthpix = 0;
%    gl.screenHeightpix = 0;
%    gl.screenCenterXpix = 0;
%    gl.screenCenterYpix = 0;
%    gl.bkgndrgb = [0 0 0];
    gl.bkgnd.tex = [];
    gl.bkgnd.on = 0;
%    gl.cal.gammaTable = []; 
%    gl.cal.monSpd = [];
%    gl.cal.fundamentals = [];
%    gl.cal.M = zeros(3);
%    gl.cal.invM = zeros(3);
    
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
    gl.bkgndtype = [];
   
    % For ring (annular) stimuli
    gl.ring.on = 0;
    gl.ring.tex = [];
    gl.ring.near_tex = [];
    gl.ring.far_tex = [];
    gl.ring.drawrect_near = [];
    gl.ring.drawrect_far = [];
    gl.ring.targets.on = 0;
    gl.ring.targets.correct = 0;
    gl.ring.targets.rgb_near = zeros(1,3);
    gl.ring.targets.rgb_far = zeros(1,3);

    % For Gabor patches
    gl.gabor.template1 = []; % orientation 1
    gl.gabor.template2 = []; % orientation 2 (orthogonal to orientation 1)

    gl.gabor.rgb = [];
    gl.gabor.tf = [];
    gl.target.on = 0;
    gl.flankers.on = 0;
    gl.null.on = 0;
    gl.blockers.on = 0;
    gl.blockers.rects = [];

    [udpCom.sock, Success] = pnetStart(udpCom.port);
    if ~(Success)
        returns
    end
    pnet(udpCom.sock, 'setreadtimeout', 0);
    pnet(udpCom.sock, 'setwritetimeout', 0);
    disp('In SMurray');

    messageIsAvailable = 0;
    while ~messageIsAvailable  % endless loop.
        messageIsAvailable = pnet(udpCom.sock, 'readpacket', 1000, 'noblock');
        if (messageIsAvailable)
            DealWithMessage(messageIsAvailable);
            messageIsAvailable = 0;   
        end     
 %       if (gl.windowPtr > 0)
            if (gl.bkgnd.on)
                DrawBkgnd();
            end
            if (gl.blockers.on)
                DrawBlockers();
            end
            if (gl.fp.on)
                DrawFP();
            end
            if (gl.ring.on)
                DrawRing();
            end
            if (gl.flankers.on || gl.target.on)
                DrawGabors();
            end
            if gl.ring.targets.on
                DrawTargets();
            end
            if (gl.fliprequest)
                DoFlip();
            end
  %      end
        if (gl.timetoleave)
            if (~isempty(gl.ring.tex))
                Screen('Close',gl.ring.tex);
                gl.ring.tex = [];
                Screen('Close',gl.bkgnd.tex);
                gl.bkgnd.tex = [];
            end
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
    if (gl.ring.on || gl.target.on || gl.flankers.on || gl.null.on)
        if(gl.framecounter == 0)
           pnet(udpCom.sock, 'write', 'MACSTIMON>> >>');
           pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
        end
        if (gl.framecounter == gl.framecountermax-1)
            pnet(udpCom.sock, 'write', 'MACSTIMOFF>> >>');
            pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
            gl.ring.on = 0;
            gl.target.on = 0;
            gl.flankers.on = 0;
            gl.null.on = 0;
        end
        gl.framecounter = gl.framecounter + 1;
    end
    gl.fliprequest = 0;    
end

%%
% DOES NOT STAND ALONE  Must inherit some parameters (e.gm fundamentals)
% from another function.
function InitDisplay(mondist, screenwidth, calfilename, bkgndtype)
    global gl;
    load(calfilename);
    cal = cals{end};
    gl.mondistcm = mondist;
    gl.screenWidthcm = screenwidth;
    gl.bkgndRGB = round(255*cal.bgColor)';
    gl.cal.gammaTable = cal.gammaTable;
    gl.cal.monSpd = cal.P_device;
    gl.cal.invgammaTable = InvertGamma(gl.cal.gammaTable, 1);
    
           
    if (isempty(Screen('windows')))
        [gl.windowPtr,gl.windowRect] = Screen('OpenWindow',0, gl.bkgndRGB);
%        gl.windowPtr = Screen('OpenWindow',0, [0 0 0]);
    else
        gl.windowPtr = Screen('windows');
        gl.windowPtr = gl.windowPtr(1);
    end
    % Gamma correction is done in software so just set hardware gamma
    % lookup to the unity line.
    clut = repmat(linspace(0,1,256),3,1)';
    Screen('LoadNormalizedGammaTable', gl.windowPtr, clut);
    
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
              
  
    gl.ring.on = 0;
    gl.fp.on = 0;
    gl.blockers.on = 0;
    
    gl.bkgndtype = bkgndtype;
    
    % ***********
    % GDLH 11/7/11 trying to add a background
    % ***********
    if (gl.usebitmapbkgnd)
        if (bkgndtype == 0)
            im = imread('/Users/horwitzlab/Desktop/SMurray pics/BlankSmall.jpg','jpg');
            gl.bkgnd.tex = Screen('MakeTexture', gl.windowPtr, im);
            Screen('DrawTexture',gl.windowPtr,gl.bkgnd.tex, [], [],[],0);
            gl.bkgnd.on = 1;
        elseif (bkgndtype == 1)
            im = imread('/Users/horwitzlab/Desktop/SMurray pics/2dWall.jpg','jpg');
            gl.bkgnd.tex = Screen('MakeTexture', gl.windowPtr, im);
            Screen('DrawTexture',gl.windowPtr,gl.bkgnd.tex, [], [],[],0);
            gl.bkgnd.on = 1;
        elseif (bkgndtype == 2)
            Screen('FillRect',gl.windowPtr,gl.bkgndRGB);
        elseif (bkgndtype == 3)
            im = imread('/Users/horwitzlab/Desktop/SMurray pics/BlankSmall.jpg','jpg');
            gl.bkgnd.tex = Screen('MakeTexture', gl.windowPtr, im);
            Screen('FillRect', gl.windowPtr, gl.bkgndRGB);
        end
        gl.fliprequest = 1;
    end
    % GDLH end of debugging stuff
    
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
    gl.fp.rgb = [fpr, fpg, fpb];
    
    fpsizeinpix = round(gl.pixperdeg*gl.fp.size);

    gl.fp.drawrect = [gl.screenCenterXpix+(gl.fp.x*gl.pixperdeg)-floor(fpsizeinpix/2)+1 ...
                gl.screenCenterYpix-(gl.fp.y*gl.pixperdeg)-floor(fpsizeinpix/2)+1 ...
                gl.screenCenterXpix+(gl.fp.x*gl.pixperdeg)+ceil(fpsizeinpix/2)...
                gl.screenCenterYpix-(gl.fp.y*gl.pixperdeg)+ceil(fpsizeinpix/2)];

    
    gl.fp.on = 1;
    
end

%%
function ShowBlockers(fp, rf)
    global gl

    nrows = size(fp, 1);    
    % might not be square because blockers bounded by fixation point    
    quadrant = xor(rf(1) < 0, rf(2) < 0) + 2 * (rf(2) < 0) + 1;
    rfpix_rel = round(gl.pixperdeg * (fp + repmat(rf, nrows, 1)) / 10);
    fppix_rel = round(gl.pixperdeg * fp / 10);
    centerpix = repmat([gl.screenCenterXpix gl.screenCenterYpix], nrows, 1);
    radius = round(5 * repmat(gl.pixperdeg, nrows, 1));
    
    width = repmat(gl.screenWidthpix, nrows, 1);
    height = repmat(gl.screenHeightpix, nrows, 1);
    
    rfpix = [centerpix(:,1) + rfpix_rel(:,1) centerpix(:,2) - rfpix_rel(:,2)];
    fppix = [centerpix(:,1) + fppix_rel(:,1) centerpix(:,2) - fppix_rel(:,2)];
    
    visible_rect = zeros(nrows, 4);
    % row 1 -> visible "near" rectangle, row 2 -> "far"
    switch quadrant
        case 1
            visible_rect = [min( max( max( fppix(:,1), rfpix(:,1) - radius ), 0 ), width ) ...
                min( max( rfpix(:,2) - radius, 0 ), height ) ...
                max( min( rfpix(:,1) + radius, width ), 0 ) ...
                max( min( min( fppix(:,2), rfpix(:,2) + radius ), height ), 0 )];
        case 2
            visible_rect = [min( max( rfpix(:,1) - radius, 0 ), width ) ...
                min( max( rfpix(:,2) - radius, 0 ), height ) ...
                max( min( min( fppix(:,1), rfpix(:,1) + radius ), width ), 0 ) ...
                max( min( min( fppix(:,2), rfpix(:,2) + radius ), height ), 0 )];            
        case 3
            visible_rect = [min( max( rfpix(:,1) - radius, 0 ), width ) ...
                max( min( max( fppix(:,2), rfpix(:,2) - radius ), height ), 0 ) ...
                max( min( min( fppix(:,1), rfpix(:,1) + radius ), width ), 0 ) ...
                max( min( rfpix(:,2) + radius, height ), 0 )];            
        case 4
            visible_rect = [min( max( max( fppix(:,1), rfpix(:,1) - radius ), 0 ), width ) ...
                max( min( max( fppix(:,2), rfpix(:,2) - radius ), height ), 0 ) ...
                max( min( rfpix(:,1) + radius, width ), 0 ) ...
                max( min( rfpix(:,2) + radius, height ), 0 )];
    end
    
    gl.blockers.rects = MakeBlockerRects(visible_rect);        

    gl.blockers.on = 1;
end

%%
function DrawBlockers()
    global gl
    
    Screen('DrawTexture', gl.windowPtr, gl.bkgnd.tex, [], [], [], 0);
    Screen('FillRect', gl.windowPtr, gl.bkgndRGB, gl.blockers.rects);
    gl.fliprequest = 1;
end

%%
function rects = MakeBlockerRects(visible_rect)
    global gl
    
    if size(visible_rect, 1) == 1 % only one stimulus
        rects = [0 0 visible_rect(1) gl.screenHeightpix;
            visible_rect(1) 0 gl.screenWidthpix visible_rect(2);
            visible_rect(3) visible_rect(2) gl.screenWidthpix visible_rect(4);
            visible_rect(1) visible_rect(4) gl.screenWidthpix gl.screenHeightpix]';
    else
        near_blocker3_right = max((visible_rect(1,2) >= visible_rect(2,4)) * gl.screenWidthpix, visible_rect(2,1));
        near_blocker4_right = max((visible_rect(1,4) >= visible_rect(2,4)) * gl.screenWidthpix, visible_rect(2,1));
        rects = [...
            0 0 visible_rect(1,1) gl.screenHeightpix;
            visible_rect(1,1) 0 visible_rect(2,1) visible_rect(1,2);
            visible_rect(1,3) visible_rect(1,2) near_blocker3_right visible_rect(1,4);
            visible_rect(1,1) visible_rect(1,4) near_blocker4_right gl.screenHeightpix;
%             now the remaining rectangles around the "far" stimulus
            visible_rect(2,1) 0 gl.screenWidthpix visible_rect(2,2);
            visible_rect(2,3) visible_rect(2,2) gl.screenWidthpix visible_rect(2,4);
            visible_rect(2,1) visible_rect(2,4) gl.screenWidthpix gl.screenHeightpix ...
]';
    end
end

%%
% "Draw" functions get called on every screen refresh so long as the
% corresponding graphical object is to be displayed.  These are the
% functions that actually draw the object to the screen (and they make a
% fliprequest so that the drawn object will actually appear).
function DrawFP()
    global gl;
    
    Screen('FillRect', gl.windowPtr, gl.fp.rgb, gl.fp.drawrect);
    
    gl.fliprequest = 1;
end

%%
function HideFP()
    global gl;
   
    gl.fp.on = 0;
    gl.blockers.on = 0;
    gl.fliprequest = 1;
end

%%
% For drawing annuli.
% Many of the passed parameters aren't used, but it's useful to have them
% passed because someday we might want to have REX change them from trial
% to trial.
function ShowRing(stimx, stimy, stimid, stimod, config, type, r, g, b, nframes)
    global gl

    if length(stimx) == 1
        stimx = stimx/10+gl.fp.x;
        stimy = stimy/10+gl.fp.y;
        stimod = stimod/10;
        stimid = stimid/10;
        rgb = [r g b];
        if (~isempty(gl.ring.tex))
            Screen('Close',gl.ring.tex);
            gl.ring.tex = [];
        end

        gl.framecountermax = nframes;

        x = stimx*gl.pixperdeg;
        y = stimy*gl.pixperdeg;
        stimsizeinpix = round(stimod*gl.pixperdeg);
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


        % Making the annulus
        interval1 = linspace(-stimod, stimod, 2*stimsizeinpix);
        interval2 = linspace(-stimod, stimod, stimsizeinpix);
        [X, Y] = meshgrid(interval2,interval1);
        dist = sqrt(X.^2+Y.^2);
        annulus = dist < stimod & dist > stimid;
        im = zeros(2*stimsizeinpix, stimsizeinpix, 3);

        for plane = 1:3
            tmp = (65535/255)*(annulus.*(rgb(plane)-gl.bkgndRGB(plane))+gl.bkgndRGB(plane));
            im(:,:,plane) = reshape(tmp, 2*stimsizeinpix, stimsizeinpix);
        end
        im(:,:,4) = 255*annulus;
        Screen('BlendFunction',gl.windowPtr, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        gl.ring.tex=Screen('MakeTexture', gl.windowPtr, im);
        % GH ignoring Bits++
        %gl.ring.tex=Screen('MakeTexture', gl.windowPtr, TranslateToColourMode(im,1));
        gl.framecounter = 0;
        gl.ring.on = 1;
    else % SMurray_p, first element = near, second = far
        stimx = stimx/10;
        stimy = stimy/10;
        stimod = stimod/10;
        stimid = stimid/10;
        rgb = [r g b];

        if ~isempty(gl.ring.near_tex)
            Screen('Close', gl.ring.near_tex);
            gl.ring.near_tex = [];
        end
        if ~isempty(gl.ring.far_tex)
            Screen('Close', gl.ring.far_tex);
            gl.ring.far_tex = [];
        end

        gl.framecountermax = nframes;

        x = stimx*gl.pixperdeg;
        y = stimy*gl.pixperdeg;
        stimsizeinpix = round(stimod*gl.pixperdeg);

        gl.ring.drawrect_near = round([gl.screenCenterXpix+x(1)-stimsizeinpix(1) ...
            gl.screenCenterYpix-y(1)-stimsizeinpix(1) ...
            gl.screenCenterXpix+x(1)+stimsizeinpix(1)...
            gl.screenCenterYpix-y(1)+stimsizeinpix(1)]);

        gl.ring.drawrect_far = round([gl.screenCenterXpix+x(2)-stimsizeinpix(2) ...
            gl.screenCenterYpix-y(2)-stimsizeinpix(2) ...
            gl.screenCenterXpix+x(2)+stimsizeinpix(2)...
            gl.screenCenterYpix-y(2)+stimsizeinpix(2)]);

        if rem(gl.ring.drawrect_near(1), 2) %if the rectangle starts on an odd pixel
            gl.ring.drawrect_near(1) = gl.ring.drawrect_near(1) - 1;
            gl.ring.drawrect_near(3) = gl.ring.drawrect_near(3) - 1;
        end

        if rem(gl.ring.drawrect_far(1), 2) %if the rectangle starts on an odd pixel
            gl.ring.drawrect_far(1) = gl.ring.drawrect_far(1) - 1;
            gl.ring.drawrect_far(3) = gl.ring.drawrect_far(3) - 1;
        end

        % 		if (gl.drawrect(3)-gl.drawrect(1) ~= 2*stimsizeinpix || gl.drawrect(4)-gl.drawrect(2) ~= 2*stimsizeinpix)
        % 			sca;
        % 			keyboard;
        % 			% Draw window is the wrong size
        % 		end

        % Making the annulus
        interval_near1 = linspace(-stimod(1), stimod(1), 2*stimsizeinpix(1));
        interval_near2 = linspace(-stimod(1), stimod(1), stimsizeinpix(1));
        interval_far1 = linspace(-stimod(2), stimod(2), 2*stimsizeinpix(2));
        interval_far2 = linspace(-stimod(2), stimod(2), stimsizeinpix(2));

        [X_near, Y_near] = meshgrid(interval_near2,interval_near1);
        dist_near = sqrt(X_near.^2+Y_near.^2);
        annulus_near = dist_near < stimod(1) & dist_near > stimid(1);
        im_near = zeros(2*stimsizeinpix(1), stimsizeinpix(1), 3);

        [X_far, Y_far] = meshgrid(interval_far2,interval_far1);
        dist_far = sqrt(X_far.^2+Y_far.^2);
        annulus_far = dist_far < stimod(2) & dist_far > stimid(2);
        im_far = zeros(2*stimsizeinpix(2), stimsizeinpix(2), 3);

        for plane = 1:3
            tmp = (65535/255)*(annulus_near.*(rgb(plane)-gl.bkgndRGB(plane))+gl.bkgndRGB(plane));
            im_near(:,:,plane) = reshape(tmp, 2*stimsizeinpix(1), stimsizeinpix(1));
        end
        im_near(:,:,4) = 255*annulus_near;
        Screen('BlendFunction',gl.windowPtr, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        for plane = 1:3
            tmp = (65535/255)*(annulus_far.*(rgb(plane)-gl.bkgndRGB(plane))+gl.bkgndRGB(plane));
            im_far(:,:,plane) = reshape(tmp, 2*stimsizeinpix(2), stimsizeinpix(2));
        end
        im_far(:,:,4) = 255*annulus_far;
        Screen('BlendFunction',gl.windowPtr, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        gl.ring.near_tex = Screen('MakeTexture', gl.windowPtr, im_near);
        gl.ring.far_tex = Screen('MakeTexture', gl.windowPtr, im_far);
        % GH ignoring Bits++
        %gl.ring.tex=Screen('MakeTexture', gl.windowPtr, TranslateToColourMode(im,1));
        gl.framecounter = 0;
        gl.ring.on = 1;
    end
end

%%
% 
function PrepareGabor(theta, sf, phi, sigma, nsigmas, gamma, stimx, stimy, contrast, config, flankerdistinsigmas)
    global gl;
    
    gl.gabor.config = config;
 
    lambda = 1/sf;
    stimsizeindeg = sigma*nsigmas; % This is half stim size (we we go out +/- nsigmas)
     if (gl.usecolourmode)
         stimsizeinpix = round(stimsizeindeg*gl.pixperdeg);  % full stim size in doublewide pixels
         [x,y] = meshgrid(stimsizeindeg*linspace(-1,1,stimsizeinpix), stimsizeindeg*linspace(-1,1,2*stimsizeinpix));
         % x and y are in dva
     else
         stimsizeinpix = round(2*stimsizeindeg*gl.pixperdeg);  % full stim size in singlewide pixels
         [x,y] = meshgrid(stimsizeindeg*linspace(-1,1,stimsizeinpix), stimsizeindeg*linspace(-1,1,stimsizeinpix));
         % x and y are in dva
     end
     
    xprime = x*cos(-theta) + y*sin(-theta);
    yprime = -x*sin(-theta) + y*cos(-theta);
    gl.gabor.template1 = exp(-(xprime.^2 + gamma.^2 .* yprime.^2)./ (2*sigma^2)).*cos(2*pi*yprime./lambda+phi); 

    xprime = x*cos(pi/2-theta) + y*sin(pi/2-theta);
    yprime = -x*sin(pi/2-theta) + y*cos(pi/2-theta);
    gl.gabor.template2 = exp(-(xprime.^2 + gamma.^2 .* yprime.^2)./ (2*sigma^2)).*cos(2*pi*yprime./lambda+phi); 

    
    stimx = stimx/10;
    stimy = stimy/10;
    
    % Creating the drawing window
    x = stimx*gl.pixperdeg;
    y = stimy*gl.pixperdeg;
    
    if ~gl.usecolourmode
        stimsizeinpix = stimsizeinpix/2;
    end
    gl.drawrect = round([gl.screenCenterXpix+x-stimsizeinpix ...
                gl.screenCenterYpix-y-stimsizeinpix ...
                gl.screenCenterXpix+x+stimsizeinpix...
                gl.screenCenterYpix-y+stimsizeinpix]);

            
    if (gl.usecolourmode)
       if(rem(gl.drawrect(1), 2)) %if the rectangle starts on an odd pixel
           gl.drawrect(1) = gl.drawrect(1) - 1;
           gl.drawrect(3) = gl.drawrect(3) - 1;
       end
       
       if (gl.drawrect(3)-gl.drawrect(1) ~= 2*stimsizeinpix || gl.drawrect(4)-gl.drawrect(2) ~= 2*stimsizeinpix)
           sca;
           keyboard;
           % Draw window is the wrong size
       end
    end
    % Preparing drawing windows for flankers
    % Not sure if I have the units right
    flankerdistinpix = sigma*gl.pixperdeg*flankerdistinsigmas;

    gl.drawrect(2,:) = gl.drawrect(1,:)+[0 flankerdistinpix 0 flankerdistinpix];  % Flanker 1
    gl.drawrect(3,:) = gl.drawrect(1,:)-[0 flankerdistinpix 0 flankerdistinpix];  % Flanker 2
    gl.drawrect(4,:) = gl.drawrect(1,:)+2*[0 flankerdistinpix 0 flankerdistinpix];  % Flanker-of-flanker 1
    gl.drawrect(5,:) = gl.drawrect(1,:)-2*[0 flankerdistinpix 0 flankerdistinpix];  % Flanker-of-flanker 2
    for i = 2:5
        if(rem(gl.drawrect(i,1), 2)) %if the rectangle starts on an odd pixel
            gl.drawrect(i,1) = gl.drawrect(i,1) - 1;
            gl.drawrect(i,3) = gl.drawrect(i,3) - 1;
        end
    end

    % Getting RGBs
    stimconecontrast = [contrast contrast contrast];
    bkgndlms = gl.cal.M*gl.bkgndrgb;
    gl.gabor.rgb = gl.cal.invM*(bkgndlms.*(1+stimconecontrast'));

    % Checking for out of gamut errors
    % Right now just squeezing it back into the gamut - need to send a
    % message to REX that this was done.
    lims = min([gl.bkgndrgb'; 1-gl.bkgndrgb']);
    if (any(abs(gl.gabor.rgb-gl.bkgndrgb) > lims'))
        scalefactor = max(abs(gl.gabor.rgb-gl.bkgndrgb)./lims');  % +(1/NGAMMASTEPS);
        gl.gabor.rgb = (gl.gabor.rgb-gl.bkgndrgb)./scalefactor+gl.bkgndrgb;
    end
end

%%
function ShowGabor(in, nframes)

    global gl;

    disp('In ShowGabor');
    
    gl.framecounter = 0;
    if (in == 1)
        gl.flankers.on = 1;
    elseif (in == 2)
        gl.target.on = 1;
    end
    gl.framecountermax = nframes;
    if (in == 0)
        gl.null.on = 1;
    end
end

%%
function ShowTargets(x, y, targsize, contrast, correct_targ)
    global gl    
    
    if correct_targ % when it's the far one
        gl.ring.targets.rgb_near = gl.bkgndRGB .* (1 - contrast/100);
        gl.ring.targets.rgb_far = zeros(1,3);
    else
        gl.ring.targets.rgb_far = gl.bkgndRGB .* (1 - contrast/100);
        gl.ring.targets.rgb_near = zeros(1,3);
    end
    x = x/10;
    y = y/10;
    targsize = targsize/10;
    
    sizeinpix = round(gl.pixperdeg*targsize);

    gl.ring.targets.drawrect_near = [gl.screenCenterXpix+(x(1)*gl.pixperdeg)-floor(sizeinpix/2)+1 ...
                gl.screenCenterYpix-(y(1)*gl.pixperdeg)-floor(sizeinpix/2)+1 ...
                gl.screenCenterXpix+(x(1)*gl.pixperdeg)+ceil(sizeinpix/2)...
                gl.screenCenterYpix-(y(1)*gl.pixperdeg)+ceil(sizeinpix/2)];
            
    gl.ring.targets.drawrect_far = [gl.screenCenterXpix+(x(2)*gl.pixperdeg)-floor(sizeinpix/2)+1 ...
                gl.screenCenterYpix-(y(2)*gl.pixperdeg)-floor(sizeinpix/2)+1 ...
                gl.screenCenterXpix+(x(2)*gl.pixperdeg)+ceil(sizeinpix/2)...
                gl.screenCenterYpix-(y(2)*gl.pixperdeg)+ceil(sizeinpix/2)];
            
    gl.ring.targets.on = 1;
end

%%
function DrawTargets()
    global gl
    Screen('FillRect', gl.windowPtr, gl.ring.targets.rgb_near, gl.ring.targets.drawrect_near);
    Screen('FillRect', gl.windowPtr, gl.ring.targets.rgb_far, gl.ring.targets.drawrect_far);
    gl.fliprequest = 1;
end

%%
function HideTargets()
    global gl
    gl.ring.targets.on = 0;
    gl.fliprequest = 1;
end

%%
function DrawRing()

    global gl;
    
    if ~isempty(gl.ring.tex)
        Screen('DrawTexture',gl.windowPtr,gl.ring.tex, [], gl.drawrect,[],0);
    elseif ~isempty(gl.ring.near_tex) && ~isempty(gl.ring.far_tex)
        Screen('DrawTexture',gl.windowPtr,gl.ring.near_tex, [], gl.ring.drawrect_near,[],0);
        Screen('DrawTexture',gl.windowPtr,gl.ring.far_tex, [], gl.ring.drawrect_far,[],0);
    end
    
    gl.fliprequest = 1;
end
%%
function DrawBkgnd()

    global gl;
    
    Screen('DrawTexture',gl.windowPtr,gl.bkgnd.tex,[],[],[],0);
    
    gl.fliprequest = 1;
end


%%
function DrawGabors()

    global gl;
    
    NGAMMASTEPS = size(gl.cal.invgammaTable,1);
 
    im = zeros(size(gl.gabor.template1,1), size(gl.gabor.template1,2), 3);
    %sca
    %keyboard
    for i = 1:2
        for plane = 1:3
            if (i == 1)
                tmp = gl.gabor.template1;
            else
                tmp = gl.gabor.template2;
            end     
            tmp = tmp.*(gl.gabor.rgb(plane)-gl.bkgndrgb(plane)) + gl.bkgndrgb(plane);
            tmp = round(tmp*(NGAMMASTEPS-1))+1;
            try
                tmp = gl.cal.invgammaTable(tmp, plane);
            catch
                sca
                keyboard
            end
            if (gl.usecolourmode)
                tmp = round(tmp*(NGAMMASTEPS-1));
            else
                tmp = round(tmp*(255-1));  % JUST FOR DEBUGING
            end
            im(:,:,plane) = reshape(tmp, size(im,1), size(im,2));
        end
        if (gl.usecolourmode)
            tex(i)=Screen('MakeTexture', gl.windowPtr, TranslateToColourModeMex(im));
        else
            tex(i)=Screen('MakeTexture', gl.windowPtr, im);
        end
    end

    if (gl.target.on)
        Screen('DrawTexture',gl.windowPtr,tex(1), [], gl.drawrect(1,:));
    end

    if (gl.flankers.on)
        if (ismember(gl.gabor.config, [2 4 6 7 9 10]))  % Flankers orthogonal to target
            Screen('DrawTexture',gl.windowPtr,tex(2), [], gl.drawrect(2,:));
            Screen('DrawTexture',gl.windowPtr,tex(2), [], gl.drawrect(3,:));
        elseif (ismember(gl.gabor.config, [1 3 5 8])) % Flankers same orientation as target
            Screen('DrawTexture',gl.windowPtr,tex(1), [], gl.drawrect(2,:));
            Screen('DrawTexture',gl.windowPtr,tex(1), [], gl.drawrect(3,:));            
        end
        if (ismember(gl.gabor.config, [6 9]))  % Flankers-of-Flankers orthogonal to target
            Screen('DrawTexture',gl.windowPtr,tex(2), [], gl.drawrect(4,:));
            Screen('DrawTexture',gl.windowPtr,tex(2), [], gl.drawrect(5,:));
        elseif (ismember(gl.gabor.config, [5 7 8 10])) % Flankers-of-Flankers same orientation as target
            Screen('DrawTexture',gl.windowPtr,tex(1), [], gl.drawrect(4,:));
            Screen('DrawTexture',gl.windowPtr,tex(1), [], gl.drawrect(5,:));            
        end
        
    end
    Screen('Close',tex(1));
    Screen('Close',tex(2));

    gl.fliprequest = 1;
end

%%
function HideRing()
    global gl;

    gl.ring.on = 0;
    gl.flankers.on = 0;
    gl.target.on = 0;
    
    gl.fliprequest = 1;
     
end
%%  
function AllOff()
    global gl;
   
    gl.ring.on = 0;
    gl.flankers.on = 0;
    gl.target.on = 0;
    gl.fp.on = 0;
    gl.ring.targets.on = 0;
    gl.fliprequest = 1;
    
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


