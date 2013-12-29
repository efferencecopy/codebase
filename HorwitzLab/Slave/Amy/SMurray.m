%#ok<*DEFNU>
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
function SMurray
% Communication variables
    global udpCom;
    udpCom.port = 6665;
    udpCom.sock = nan;
    udpCom.rexip = '192.168.1.120';

    % Keypress constants
    global KEY;
    KbName('UnifyKeyNames');
    KEY.ESC = KbName('escape');

    % Variables that are destructively modified by subfunctions
    % Because we are guaranteed to be running this code after having
    % already
    % setup the global gl structure (via the white noise paradigm)
    % there's no sense in reassigning everything.
    global gl;
    gl.usebitmapbkgnd = 1; % For debugging
    gl.bkgnd.tex = [];
    gl.minibkgndnear.tex = [];
    gl.minibkgndfar.tex = [];
    gl.bkgnd.on = 0;

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
    
    gl.vpixx = ConnectedToViewPixx();
    gl.ccmode = 1;

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
    gl.blockers.angle = 0;
    gl.blockers.rect = [];

    gl.rfpix = [0 0];
    gl.fppix = [0 0];

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
        if gl.fp.on
            DrawFP();
        end
        if (gl.blockers.on)
            DrawBlockers();
        end
        if (gl.ring.on)
            DrawRing();
        end
        if gl.fp.on
            DrawFP();
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
                if ~isempty(gl.bkgnd.tex)
                    Screen('Close',gl.bkgnd.tex);
                    gl.bkgnd.tex = [];
                end
            end
            if ~isempty(gl.minibkgndnear.tex)
                Screen('Close', gl.minibkgndnear.tex);
                gl.minibkgndnear.tex = [];
            end
            if ~isempty(gl.minibkgndfar.tex)
                Screen('Close', gl.minibkgndfar.tex);
                gl.minibkgndfar.tex = [];
            end
            return
        end
        [keyisdown,secs,keycode] = KbCheck();
        if (keyisdown && keycode(KEY.ESC))
            pnet(udpCom.sock, 'close');
            sca();
            return
        end
    end
end

%%
% DoFlip() taken from WhiteNoise
function DoFlip()
    global gl udpCom

    Screen('Flip',gl.windowPtr);
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
            gl.blockers.on = 0;
        end
        gl.framecounter = gl.framecounter + 1;
    end
    gl.fliprequest = 0;
end

%%
% DOES NOT STAND ALONE  Must inherit some parameters (e.g., fundamentals)
% from another function.
function InitDisplay(mondist, screenwidth, calfilename, bkgndtype)
    global gl
    load(calfilename);
    cal = cals{end};
    gl.mondistcm = mondist;
    gl.screenWidthcm = screenwidth;
%     gl.bkgndRGB = round(255*cal.bgColor)';
    gl.bkgndrgb = cal.bgColor;
    gl.cal.gammaTable = cal.gammaTable;
    gl.cal.monSpd = cal.P_device;
    gl.cal.invgammaTable = InvertGammaTable(cal.gammaInput, cal.gammaTable, 2^16);
    
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
    
    % For 'DrawTexture'
    Screen('BlendFunction', gl.windowPtr, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    gl.framerate = Screen('NominalFrameRate', gl.windowPtr, 1);
    [screenwidthpix, screenheightpix]  = Screen('WindowSize', gl.windowPtr);
    gl.screenWidthpix = screenwidthpix;
    gl.screenHeightpix = screenheightpix;
    gl.screenCenterXpix = screenwidthpix/2;
    gl.screenCenterYpix = screenheightpix/2;
    
    pixpercm = gl.screenWidthpix/gl.screenWidthcm;
    theta = atan2(gl.screenWidthcm/2, gl.mondistcm)*180/pi;
    cmperdeg = gl.screenWidthcm/(2*theta);
    gl.pixperdeg = pixpercm*cmperdeg;
    
%     gl.bkgndrgb = [gl.cal.gammaTable(gl.bkgndRGB(1)+1,1); ...
%         gl.cal.gammaTable(gl.bkgndRGB(2)+1,2); ...
%         gl.cal.gammaTable(gl.bkgndRGB(3)+1,3)];
    
    gl.ring.on = 0;
    gl.fp.on = 0;
    gl.blockers.on = 0;
    
    gl.bkgndtype = bkgndtype;
    
    if gl.usebitmapbkgnd
        if (bkgndtype == 0 || bkgndtype == 3)
            im = imread('BlankSmall.jpg');
            gl.bkgnd.tex = Screen('MakeTexture', gl.windowPtr, im);
            Screen('DrawTexture',gl.windowPtr,gl.bkgnd.tex, [], [],[],0);
            gl.bkgnd.on = 1;
        elseif (bkgndtype == 1)
            im = imread('2dWall.jpg');
            gl.bkgnd.tex = Screen('MakeTexture', gl.windowPtr, im);
            Screen('DrawTexture',gl.windowPtr,gl.bkgnd.tex, [], [],[],0);
            gl.bkgnd.on = 1;
        elseif (bkgndtype == 2)
            Screen('FillRect', gl.windowPtr, gl.bkgndrgb);
        elseif (bkgndtype == 4)
            im = imread('LB_jpeg_3.jpg');
            gl.bkgnd.tex = Screen('MakeTexture', gl.windowPtr, im);
            Screen('DrawTexture', gl.windowPtr, gl.bkgnd.tex, [], [], [], 0);
            gl.bkgnd.on = 1;
        elseif (bkgndtype == 5)
            gl.bkgndimage = imread('BlankSmall.jpg');
        end
        gl.fliprequest = 1;
    end
    HideCursor();
end

%%
% "Show" functions are called from REX.  They set up the appropriate fields
% in the "gl" structure (and set the "...on" toggle field to 1).
function ShowFP(x, y, size, fpr, fpg, fpb)
    global gl

    gl.fp.x = x/10;
    gl.fp.y = y/10;
    gl.fp.size = size/10;
    gl.fp.rgb = [fpr, fpg, fpb]/255;

    fpsizeinpix = round(gl.pixperdeg*gl.fp.size);

    gl.fp.drawrect = [gl.screenCenterXpix+(gl.fp.x*gl.pixperdeg)-floor(fpsizeinpix/2)+1 ...
        gl.screenCenterYpix-(gl.fp.y*gl.pixperdeg)-floor(fpsizeinpix/2)+1 ...
        gl.screenCenterXpix+(gl.fp.x*gl.pixperdeg)+ceil(fpsizeinpix/2) ...
        gl.screenCenterYpix-(gl.fp.y*gl.pixperdeg)+ceil(fpsizeinpix/2)];

    gl.fp.on = 1;
end

%%
function ShowBlockers(fp, rf)
    global gl

    rfpix_rel = round(gl.pixperdeg * (fp + rf) / 10);
    fppix_rel = round(gl.pixperdeg * fp / 10);
    centerpix = [gl.screenCenterXpix gl.screenCenterYpix];

    gl.rfpix = [centerpix(:,1) + rfpix_rel(:,1) centerpix(:,2) - rfpix_rel(:,2)];
    gl.fppix = [centerpix(:,1) + fppix_rel(:,1) centerpix(:,2) - fppix_rel(:,2)];

    gl.blockers.angle = atan2(gl.rfpix(2) - gl.fppix(2), gl.rfpix(1) - gl.fppix(1)) * 180 / pi;
    gl.blockers.rect = round(1.66 * gl.pixperdeg / 2) * [-1 -1 1 1]; % approx. size of RF
    gl.blockers.rect = CenterRectOnPoint(gl.blockers.rect, gl.rfpix(1), gl.rfpix(2));

    if rem(gl.blockers.rect(1), 2)
        gl.blockers.rect([1 3]) = gl.blockers.rect([1 3]) - 1;
    end

    gl.blockers.on = 1;
end

%%
function DrawBlockers()
    global gl

    blockers_rgb = [178 88 98]/255;
    if all(abs(gl.blockers.rect) == abs(gl.blockers.rect(1))) % circle: no rotation required
        Screen('FillOval', gl.windowPtr, blockers_rgb, gl.blockers.rect);
    else
        Screen('glPushMatrix', gl.windowPtr);
        Screen('glTranslate', gl.windowPtr, gl.rfpix(1), gl.rfpix(2));
        Screen('glRotate', gl.windowPtr, gl.blockers.angle);
        Screen('glTranslate', gl.windowPtr, -gl.rfpix(1), -gl.rfpix(2));
        Screen('FillOval', gl.windowPtr, blockers_rgb, gl.blockers.rect);
        Screen('glPopMatrix', gl.windowPtr);
    end

    gl.fliprequest = 1;
end

%%
% "Draw" functions get called on every screen refresh so long as the
% corresponding graphical object is to be displayed.  These are the
% functions that actually draw the object to the screen (and they make a
% fliprequest so that the drawn object will actually appear).
function DrawFP()
    global gl
    Screen('FillRect', gl.windowPtr, gl.fp.rgb, gl.fp.drawrect);
    gl.fliprequest = 1;
end

%%
function HideFP()
    global gl
    gl.fp.on = 0;
    gl.fliprequest = 1;
end

%%
% For drawing annuli.
% Many of the passed parameters aren't used, but it's useful to have them
% passed because someday we might want to have REX change them from trial
% to trial.
function ShowRing(stimx, stimy, stimir, stimor, bkgndtype, maxor, r, g, b, nframes)
    global gl

    if length(stimx) == 1
        stimx = stimx/10+gl.fp.x;
        stimy = stimy/10+gl.fp.y;
        stimor = stimor/10;
        stimir = stimir/10;
        rgb = [r g b]'/255;

        if (~isempty(gl.ring.tex))
            Screen('Close',gl.ring.tex);
            gl.ring.tex = [];
        end

        gl.framecountermax = nframes;

        x = stimx*gl.pixperdeg;
        y = stimy*gl.pixperdeg;
        stimsizeinpix = round(stimor*gl.pixperdeg);

        gl.drawrect = round([gl.screenCenterXpix+x-stimsizeinpix ...
            gl.screenCenterYpix-y-stimsizeinpix ...
            gl.screenCenterXpix+x+stimsizeinpix ...
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
        interval1 = linspace(-stimor, stimor, 2*stimsizeinpix);
        interval2 = linspace(-stimor, stimor, stimsizeinpix);
        [X, Y] = meshgrid(interval2,interval1);
        dist = sqrt(X.^2+Y.^2);
        annulus = dist < stimor & dist > stimir;
        im = zeros(2*stimsizeinpix, stimsizeinpix, 4);

        im(:,:,1:3) = bsxfun(@plus, bsxfun(@times, annulus, reshape(rgb(:)-gl.bkgndrgb(:), [1 1 3])), ...
            reshape(gl.bkgndrgb, [1 1 3]));
        im(:,:,4) = annulus;

        if bkgndtype == 5 && (isempty(gl.minibkgndnear.tex) || isempty(gl.minibkgndfar.tex))
            maxor = maxor/10;
            or_pix = round(maxor*gl.pixperdeg);
            interval1 = linspace(-maxor, maxor, 2*or_pix);
            [X,Y] = meshgrid(interval1);
            dist = sqrt(X.^2+Y.^2);
            annulus_filled = dist < maxor;

            dim_image = [size(gl.bkgndimage, 1) size(gl.bkgndimage, 2)]; % height width
            annulus_bounds = or_pix * [-1 -1 1 1];
            annulus_bounds = round(CenterRectOnPoint(annulus_bounds, dim_image(2)/2+x, dim_image(1)/2-y));

            alpha_channel = zeros(dim_image(1), dim_image(2));
            alpha_channel(annulus_bounds(2):annulus_bounds(4)-1,annulus_bounds(1):annulus_bounds(3)-1) = ...
                annulus_filled;

            if stimx > 0 % far
                gl.minibkgndfar.tex = Screen('MakeTexture', gl.windowPtr, cat(3,gl.bkgndimage,alpha_channel), [], [], 2);
            else % near
                gl.minibkgndnear.tex = Screen('MakeTexture', gl.windowPtr, cat(3,gl.bkgndimage,alpha_channel), [], [], 2);
            end
        end

        gl.ring.tex=Screen('MakeTexture', gl.windowPtr, im, [], [], 2);
        gl.framecounter = 0;
        gl.ring.on = 1;
    else % SMurray_p, first element = near, second = far
        stimx = stimx/10;
        stimy = stimy/10;
        stimor = stimor/10;
        stimir = stimir/10;
        rgb = [r g b]/255;

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
        stimsizeinpix = round(stimor*gl.pixperdeg);

        gl.ring.drawrect_near = round([gl.screenCenterXpix+x(1)-stimsizeinpix(1) ...
            gl.screenCenterYpix-y(1)-stimsizeinpix(1) ...
            gl.screenCenterXpix+x(1)+stimsizeinpix(1)...
            gl.screenCenterYpix-y(1)+stimsizeinpix(1)]);

        gl.ring.drawrect_far = round([gl.screenCenterXpix+x(2)-stimsizeinpix(2) ...
            gl.screenCenterYpix-y(2)-stimsizeinpix(2) ...
            gl.screenCenterXpix+x(2)+stimsizeinpix(2)...
            gl.screenCenterYpix-y(2)+stimsizeinpix(2)]);

        if rem(gl.ring.drawrect_near(1), 2) % if the rectangle starts on an odd pixel
            gl.ring.drawrect_near([1 3]) = gl.ring.drawrect_near([1 3]) - 1;
        end

        if rem(gl.ring.drawrect_far(1), 2) % if the rectangle starts on an odd pixel
            gl.ring.drawrect_far([1 3]) = gl.ring.drawrect_far([1 3]) - 1;
        end

        % 		if (gl.drawrect(3)-gl.drawrect(1) ~= 2*stimsizeinpix || gl.drawrect(4)-gl.drawrect(2) ~= 2*stimsizeinpix)
        % 			sca;
        % 			keyboard;
        % 			% Draw window is the wrong size
        % 		end

        % Making the annulus
        interval_near1 = linspace(-stimor(1), stimor(1), 2*stimsizeinpix(1));
        interval_near2 = linspace(-stimor(1), stimor(1), stimsizeinpix(1));
        interval_far1 = linspace(-stimor(2), stimor(2), 2*stimsizeinpix(2));
        interval_far2 = linspace(-stimor(2), stimor(2), stimsizeinpix(2));

        [X_near, Y_near] = meshgrid(interval_near2,interval_near1);
        dist_near = sqrt(X_near.^2+Y_near.^2);
        annulus_near = dist_near < stimor(1) & dist_near > stimir(1);
        im_near = zeros(2*stimsizeinpix(1), stimsizeinpix(1), 4);

        [X_far, Y_far] = meshgrid(interval_far2,interval_far1);
        dist_far = sqrt(X_far.^2+Y_far.^2);
        annulus_far = dist_far < stimor(2) & dist_far > stimir(2);
        im_far = zeros(2*stimsizeinpix(2), stimsizeinpix(2), 4);

        im_near(:,:,1:3) = bsxfun(@plus, bsxfun(@times, annulus_near, reshape(rgb(:)-gl.bkgndrgb(:), [1 1 3])), ...
            reshape(gl.bkgndrgb, [1 1 3]));
        im_near(:,:,4) = annulus_near;

        im_far(:,:,1:3) = bsxfun(@plus, bsxfun(@times, annulus_far, reshape(rgb(:)-gl.bkgndrgb(:), [1 1 3])), ...
            reshape(gl.bkgndrgb, [1 1 3]));
        im_far(:,:,4) = annulus_far;

        gl.ring.near_tex = Screen('MakeTexture', gl.windowPtr, im_near, [], [], 2);
        gl.ring.far_tex = Screen('MakeTexture', gl.windowPtr, im_far, [], [], 2);

        gl.framecounter = 0;
        gl.ring.on = 1;
    end
end

%%
%
function PrepareGabor(theta, sf, phi, sigma, nsigmas, gamma, stimx, stimy, contrast, config, flankerdistinsigmas)
    global gl

    gl.gabor.config = config;

    lambda = 1/sf;
    stimsizeindeg = sigma*nsigmas; % This is half stim size (we we go out +/- nsigmas)
    stimsizeinpix = round(stimsizeindeg*gl.pixperdeg);  % full stim size in doublewide pixels
    [x,y] = meshgrid(stimsizeindeg*linspace(-1,1,stimsizeinpix), stimsizeindeg*linspace(-1,1,2*stimsizeinpix));

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
    global gl

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
        gl.ring.targets.rgb_near = (gl.bkgndrgb .* (1 - contrast/100));
        gl.ring.targets.rgb_far = zeros(1,3);
    else
        gl.ring.targets.rgb_far = (gl.bkgndrgb .* (1 - contrast/100));
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
    global gl

    if gl.bkgndtype == 5
        if gl.fp.x < 0 && ~isempty(gl.minibkgndnear.tex)
            Screen('DrawTexture', gl.windowPtr, gl.minibkgndnear.tex, [], [], [], 0);
        elseif gl.fp.x > 0 && ~isempty(gl.minibkgndfar.tex)
            Screen('DrawTexture', gl.windowPtr, gl.minibkgndfar.tex, [], [], [], 0);
        end
    end

    if ~isempty(gl.ring.tex)
        Screen('DrawTexture', gl.windowPtr, gl.ring.tex, [], gl.drawrect, [], 0);
    elseif ~isempty(gl.ring.near_tex) && ~isempty(gl.ring.far_tex)
        Screen('DrawTexture', gl.windowPtr, gl.ring.near_tex, [], gl.ring.drawrect_near, [], 0);
        Screen('DrawTexture', gl.windowPtr, gl.ring.far_tex, [], gl.ring.drawrect_far, [], 0);
    end

    gl.fliprequest = 1;
end
%%
function DrawBkgnd()
    global gl
    Screen('DrawTexture',gl.windowPtr,gl.bkgnd.tex,[],[],[],0);
    gl.fliprequest = 1;
end

%%
function DrawGabors()
    global gl

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
            try
                tmp = gl.cal.invgammaTable(tmp, plane);
            catch
                sca
                keyboard
            end
            im(:,:,plane) = reshape(tmp, size(im,1), size(im,2));
        end
        tex(i)=Screen('MakeTexture', gl.windowPtr, TranslateToColourModeMex(im, [], gl.ccmode), [], [], 2);
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
    global gl

    gl.ring.on = 0;
    gl.flankers.on = 0;
    gl.target.on = 0;

    gl.fliprequest = 1;
end

%%
function AllOff()
    global gl
    gl.ring.on = 0;
    gl.blockers.on = 0;
    gl.flankers.on = 0;
    gl.target.on = 0;
    gl.fp.on = 0;
    gl.ring.targets.on = 0;
    gl.fliprequest = 1;
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
        eval(message)
    catch ME
        fprintf('Ignoring uninterpretable message: "%s"\n',message);
%         error = lasterror;
%         disp(error.message);
%         disp(error.identifier);
%         disp(error.stack);
        disp(getReport(ME));
    end
end

