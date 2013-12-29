function GridLMPlane_Slave
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

    gl.fp.on = 0;
    gl.fp.x = 0;
    gl.fp.y = 0;
    gl.fp.size = 0;
    gl.fp.rgb = [0 0 0];
    gl.fp.drawrect = [0 0 0 0];

    gl.stim.on = 0;
    gl.stim.drawrect = [0 0 0 0];
    
    [udpCom.sock, Success] = pnetStart(udpCom.port);
    if ~Success, return, end
    
    pnet(udpCom.sock, 'setreadtimeout', 0);
    pnet(udpCom.sock, 'setwritetimeout', 0);

    messageIsAvailable = 0;
    while ~messageIsAvailable  % endless loop.
        messageIsAvailable = pnet(udpCom.sock, 'readpacket', 1000, 'no block');
        
        if messageIsAvailable
            DealWithMessage(messageIsAvailable);
            messageIsAvailable = 0;
        end

        if gl.stim.on
            DrawGaussianBar();
        end
        if gl.fp.on
            DrawFP();
        end
        if gl.fliprequest
            DoFlip();
        end

        if gl.timetoleave, return, end
        
        [keyisdown,junk,keycode] = KbCheck();
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
    global gl
    
    Screen('Flip', gl.windowPtr, 0, 0);
    gl.fliprequest = 0;
end

% DOES NOT STAND ALONE  Must inherit some parameters (e.g m fundamentals)
% from another function.
function InitDisplay(mondist, screenwidth, calfilename, fundfilename)
    global gl
    
    load(calfilename);
    cal = cals{end};
    gl.mondistcm = mondist;
    gl.screenWidthcm = screenwidth;
    gl.bkgndRGB = round(255*cal.bgColor)';
    gl.cal.gammaTable = cal.gammaTable;
    gl.cal.monSpd = cal.P_device;
    gl.cal.invgammaTable = InvertGamma(gl.cal.gammaTable, 1);

    s = load(fundfilename);
    fns = fieldnames(s);
    P_device = SplineSpd(SToWls(cal.S_device), cal.P_device, SToWls([380 5 81]));
    gl.cal.fundamentals = eval(['s.',fns{1}])';
    gl.cal.M = gl.cal.fundamentals'*P_device;
    gl.cal.invM = inv(gl.cal.M);
    
    gl.bkgndrgb = [gl.cal.gammaTable(gl.bkgndRGB(1)+1,1),...
                   gl.cal.gammaTable(gl.bkgndRGB(2)+1,2),...
                   gl.cal.gammaTable(gl.bkgndRGB(3)+1,3)]';

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

    gl.fp.on = 0;
    gl.targ.on = 0;  % probably want to change this

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

% "Draw" functions get called on every screen refresh so long as the
% corresponding graphical object is to be displayed.  These are the
% functions that actually draw the object to the screen (and they make a
% fliprequest so that the drawn object will actually appear).
function DrawFP()
    global gl

    Screen('Fillrect', gl.windowPtr, gl.fp.rgb, gl.fp.drawrect);
    gl.fliprequest = 1;
    
end

function HideFP()
    global gl

    gl.fp.on = 0;
    gl.fliprequest = 1;
    
end


function PrepareGaussianBar(pos_x, pos_y, theta, sig_x, sig_y, driftRate, nstd, nexpanse, L, M, S)
    global gl udpCom
    
    % User Defined Variables
    sintoGauss = 8; % aka, 4 Gaussian standard deviations under 1/2 a sin wave

    % Create Stimulus Space
    windowsize = round(gl.pixperdeg*sig_x*nstd/2);
    colidx = -windowsize:2:windowsize-1;
    rowidx = -windowsize:windowsize-1;
    [X, Y] = meshgrid(colidx, rowidx);
    
    % Create bar stimulus (Gaussian)
    xprime = X .* cos(-theta) + Y .* sin(-theta);
    yprime = -X .* sin(-theta) + Y .* cos(-theta);
    gaussian = exp(-xprime.^2 ./ (2* (gl.pixperdeg*sig_x).^2) - yprime.^2./(2* (gl.pixperdeg*sig_y).^2));
    stimspatialprofile =  gaussian;
    
    % Defining a drawing rectangle where the stimulus will be shown
    % [left, top, right, bottom]
    gl.stim.drawrect = round([ ...
        gl.screenCenterXpix + (pos_x*gl.pixperdeg) + min(colidx) - (sin(theta) * (nexpanse)/2*gl.pixperdeg) ...
        gl.screenCenterYpix - (pos_y*gl.pixperdeg) + min(rowidx) - (cos(theta) * (nexpanse)/2*gl.pixperdeg) ...
        gl.screenCenterXpix + (pos_x*gl.pixperdeg) - min(colidx) - (sin(theta) * (nexpanse)/2*gl.pixperdeg) ...
        gl.screenCenterYpix - (pos_y*gl.pixperdeg) - min(rowidx) - (cos(theta) * (nexpanse)/2*gl.pixperdeg) ]);

    % Shifting the drawing rectangle if it happens to start on an
    % even-numbered pixel. This is only important when using the Bits++
    if rem(gl.stim.drawrect(1), 2)
        gl.stim.drawrect(1) = gl.stim.drawrect(1) - 1;
        gl.stim.drawrect(3) = gl.stim.drawrect(3) - 1;
    end
    
    % Creating an NxMx3 image to be displayed.
    im = zeros(numel(rowidx), numel(colidx), 3); % preallocating space
    bkgndlms = gl.cal.M * gl.bkgndrgb;
    stimulusrgb = gl.cal.invM*(bkgndlms'.* (1 + [L M S]))'; % computing rgbs at the peak

    gl.stim.nframes = round(gl.framerate / (driftRate * sintoGauss * sig_y) * nexpanse);
    %keyboard; sca;
    %gl.stim.nframes = round(gl.framerate / driftRate * nexpanse);
    dvaadvX = ((driftRate * sintoGauss * sig_y) / gl.framerate * sin(theta)) * (0:gl.stim.nframes);
    dvaadvY = ((driftRate * sintoGauss * sig_y) / gl.framerate * cos(theta)) * (0:gl.stim.nframes);
    gl.stim.advvectX = 2 * round(dvaadvX * gl.pixperdeg /2);
    gl.stim.advvectY = round(dvaadvY * gl.pixperdeg);
    
    stim.nframesramp = round(gl.stim.nframes/4);
    stim.rampup = linspace(0,1,stim.nframesramp);
    stim.plat = ones(1,gl.stim.nframes - 2*size(stim.rampup,2));
    stim.ramponmovie = nan(1,size(stim.rampup,2));
    for i=1:size(stim.rampup,2)
        for gun = 1:3
            % rescaling stimspatialprofile to go between bkgndrgb and stimulusrgb (instead of between 0 and 1)
            tmp = (stimulusrgb(gun)-gl.bkgndrgb(gun))*(stimspatialprofile*stim.rampup(i))+gl.bkgndrgb(gun);
            %tmp = (stimulusrgb(gun)-gl.bkgndrgb(gun))*(stimspatialprofile)+gl.bkgndrgb(gun);
            % Rounding the intensity values at each pixel so that we can use
            % these values as indices into the inverse gamma table
            tmp = round(tmp*size(gl.cal.invgammaTable,1))+1;
            if max(max(tmp)) > size(gl.cal.invgammaTable,1)
                disp('Out of gammut!')
                sca; keyboard;
            end
            % converting intensities to voltages
            tmp = gl.cal.invgammaTable(tmp, gun);
            tmp = round(tmp*(size(gl.cal.invgammaTable,1)-1));
            im(:,:,gun) = reshape(tmp, numel(rowidx), numel(colidx));
        end
        stim.ramponmovie(i) = Screen('MakeTexture', gl.windowPtr, TranslateToColourMode(im,1));
    end
    for gun = 1:3
        % rescaling stimspatialprofile to go between bkgndrgb and stimulusrgb (instead of between 0 and 1)
        tmp = (stimulusrgb(gun)-gl.bkgndrgb(gun))*(stimspatialprofile)+gl.bkgndrgb(gun);
        % Rounding the intensity values at each pixel so that we can use
        % these values as indices into the inverse gamma table
        tmp = round(tmp*size(gl.cal.invgammaTable,1))+1;
        % converting intensities to voltages
        tmp = gl.cal.invgammaTable(tmp, gun);
        tmp = round(tmp*(size(gl.cal.invgammaTable,1)-1));
        im(:,:,gun) = reshape(tmp, numel(rowidx), numel(colidx));
    end
    stim.platmovie = Screen('MakeTexture', gl.windowPtr, TranslateToColourMode(im,1)) * stim.plat;
    gl.stim.textureidx = [stim.ramponmovie stim.platmovie fliplr(stim.ramponmovie)];
    gl.stim.framecounter = 1;
    
    pnet(udpCom.sock, 'write', 'macdone>> >>');
    pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
    %disp('macdone sent from prep')
    
end


function DisplayGaussianBar()
    global gl
    gl.stim.on = 1;
end


function DrawGaussianBar()
    global gl udpCom
    
    i = gl.stim.framecounter;
    Screen('DrawTexture',gl.windowPtr,gl.stim.textureidx(i),[],gl.stim.drawrect+[gl.stim.advvectX(i) gl.stim.advvectY(i) gl.stim.advvectX(i) gl.stim.advvectY(i)],[],0);
    
    if i < gl.stim.nframes
       gl.stim.framecounter = gl.stim.framecounter + 1;
       gl.fliprequest = 1;
    else
        gl.stim.on = 0;
        textureidx = unique(gl.stim.textureidx);
        for n = 1:numel(textureidx)
            Screen('Close',textureidx(n))
        end
        gl.stim.textureidx = [];
        pnet(udpCom.sock, 'write', 'macdone>> >>');
        pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
        %disp('macdone sent from draw')
    end
    
    
end

function AllOff()
global gl udpCom

    gl.stim.on = 0;
    gl.fp.on = 0;
    gl.fliprequest = 1;
    
    if ~isempty(gl.stim.textureidx)
        textureidx = unique(gl.stim.textureidx);
        for n = 1:numel(textureidx)
            Screen('Close',textureidx(n))
        end
    end

    gl.stim.textureidx = [];
    
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
        disp('in GridLMPlane_Slave')
        disp(ME.message);
        disp(ME.identifier);
        disp(ME.stack);
    end
end