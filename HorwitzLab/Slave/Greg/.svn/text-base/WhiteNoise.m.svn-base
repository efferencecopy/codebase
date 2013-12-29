function WhiteNoise
% WhiteNoise.m
%
% Generate white noise stimuli.
%
% GDLH 11/29/07
%
% Porting over to work with pnet instead of matlabUDP which
% is crashing surprisingly often for some reason.
%
% GDLH 12/22/07

    % Communication variables
    global udpCom
    udpCom.port = 6665;
    udpCom.sock = nan;
    udpCom.rexip = '192.168.1.120';

    % Keypress constants
    global KEY
    KbName('UnifyKeyNames');
    KEY.ESC = KbName('escape');

    % Variables that are destructively modified by subfunctions
    global gl
    gl.mondistcm = 0;
    gl.screenWidthcm = 0;
    gl.screenWidthpix = 0;
    gl.screenHeightpix = 0;
    gl.screenCenterXpix = 0;
    gl.screenCenterYpix = 0;
    gl.pixperdeg = 0;
    gl.bkgndRGB = [0 0 0];
    gl.windowPtr = 0;
    gl.fliprequest = 0;
    gl.gaussgamma = [];
    gl.framerate = 0;
    gl.lmsbinaryrgbmat = [];
    gl.timetoleave = 0;
    
    gl.vpixx = ConnectedToViewPixx();
    gl.ccmode = 1;

    gl.fp.on = 0;
    gl.fp.x = 0;
    gl.fp.y = 0;
    gl.fp.size = 0;
    gl.fp.rgb = [0 0 0];
    gl.fp.drawrect = [0 0 0 0];

    gl.stim.on = 0;
    gl.stim.x = 0;
    gl.stim.y = 0;
    gl.stim.seed = 1; % A bit of hack
    gl.stim.npixperstix = 0;
    gl.stim.nstixperside = 0;
    gl.stim.type = 'RGB';
    gl.stim.mu = [0 0 0];
    gl.stim.sigma = [1 1 1];
    gl.stim.drawrect = [0 0 0 0];
    gl.stim.framecounter = 0;
    gl.stim.template = [];

    gl.synthimage.im = [];
    gl.synthimage.on = 0;

    gl.gabor.rawrgb = [0 0 0];
    gl.gabor.theta = 0;
    gl.gabor.lambda = 0;
    gl.gabor.phi = 0;
    gl.gabor.sigma = 0;
    gl.gabor.gamma = 0;
    gl.gabor.xoffset = 0;
    gl.gabor.yoffset = 0;

    gl.cal.gammaTable = [];
    gl.cal.invGamma = [];
    gl.cal.monSpd = [];
    gl.cal.fundamentals = [];
    gl.cal.M = zeros(3);
    gl.cal.invM = zeros(3);

    %Start by opening a connectionless udp socket
    udpCom.sock = pnetStart(udpCom.port);
    if udpCom.sock == -1
        error('UDP not properly initialized')
    end
    pnet(udpCom.sock, 'setreadtimeout', 0);
    pnet(udpCom.sock, 'setwritetimeout', 0);

    messageIsAvailable = 0;
    while ~messageIsAvailable  % endless loop.
        messageIsAvailable = pnet(udpCom.sock, 'readpacket', 350, 'noblock');
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
            if gl.synthimage.on
                DrawImage();
            end
            if gl.fliprequest
                DoFlip();
            end
        end
        if gl.timetoleave
            return
        end
        [keyisdown,secs,keycode] = KbCheck();
        if keyisdown && keycode(KEY.ESC)
            ShowCursor();
            pnet(udpCom.sock, 'close');
            sca();
            return
        end
    end
end

%%
function DoFlip()
    global gl udpCom

    Screen('Flip',gl.windowPtr,0,0);
    if gl.stim.on || gl.synthimage.on
        if gl.stim.framecounter == 0
            pnet(udpCom.sock, 'write', 'MACSTIMON>> >>');
            pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
        end
        gl.stim.framecounter = gl.stim.framecounter + 1;
    end
    gl.fliprequest = 0;
end

%%
function InitDisplay(mondist, screenwidth, varargin) % varargin{1} and {2} are defined to be the calibration and fundamentals
    global gl

    calfilename = varargin{1};
    load(calfilename);
    cal = cals{end}; %#ok<USENS>
    gl.mondistcm = mondist;
    gl.screenWidthcm = screenwidth;
    gl.bkgndRGB = round(255*cal.bgColor)';
    gl.cal.gammaTable = cal.gammaTable;
    gl.cal.monSpd = cal.P_device;

    gl.cal.invGamma = InvertGammaTable(cal.gammaInput, gl.cal.gammaTable, 2^16);
    
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
    gl.screenWidthpix = screenwidthpix;
    gl.screenHeightpix = screenheightpix;
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
    global gl

    Screen('FillRect', gl.windowPtr, gl.fp.rgb, gl.fp.drawrect);
    gl.fliprequest = 1;
end

%%
function ShowStim(rfx, rfy, seed, npixperstix, nstixperside, mu1, mu2, mu3, sigma1, sigma2, sigma3, gausslocut, gausshicut, noisetype)
    global gl

    NGAMMASTEPS = size(gl.cal.invGamma,1);

    gl.stim.x = rfx/10;
    gl.stim.y = rfy/10;
    gl.bar.xy = [gl.stim.x gl.stim.y]; % Updating the RF location estimate in case we adjust it manually

    if isempty(gl.synthimage.im)
        disp('Initializing synthimage');
        gl.synthimage.im = zeros(nstixperside, nstixperside, 3);
    end

    gl.stim.seed = seed;
    gl.stim.npixperstix = npixperstix;
    gl.stim.nstixperside = nstixperside;

    gl.stim.sigma = [sigma1 sigma2 sigma3]/1000;

    gl.stim.mu = [gl.cal.gammaTable(gl.bkgndRGB(1)+1,1),...
        gl.cal.gammaTable(gl.bkgndRGB(2)+1,2),...
        gl.cal.gammaTable(gl.bkgndRGB(3)+1,3)]...
        +[mu1 mu2 mu3]/1000;

    switch noisetype % keep this in sync with the REX paradigm code
        case 1
            gl.stim.type = 'RGB';
        case {2,3}
            gl.stim.type = 'LMS';
        otherwise
            sca();
            error('unrecognized noise type');
    end

    stimsizeinpix = gl.stim.nstixperside*gl.stim.npixperstix;
    gl.stim.drawrect = round([gl.screenCenterXpix+(gl.stim.x*gl.pixperdeg)-floor(stimsizeinpix) ...
        gl.screenCenterYpix-(gl.stim.y*gl.pixperdeg)-floor(stimsizeinpix)-1 ...
        gl.screenCenterXpix+(gl.stim.x*gl.pixperdeg)+ceil(stimsizeinpix) ...
        gl.screenCenterYpix-(gl.stim.y*gl.pixperdeg)+ceil(stimsizeinpix)]);

    if rem(gl.stim.drawrect(1), 2) %if the rectangle starts on an odd pixel
        gl.stim.drawrect(1) = gl.stim.drawrect(1) - 1;
        gl.stim.drawrect(3) = gl.stim.drawrect(3) - 1;
    end

    if strcmp(gl.stim.type,'RGB') % Gaussian gun noise
        % Preparing the combined Gaussian/Gamma functions
        if max(norminv(gausshicut/1000) * gl.stim.sigma + gl.stim.mu) > 1 || ...
           min(norminv(gausshicut/1000) * gl.stim.sigma + gl.stim.mu) < .5 / NGAMMASTEPS
            sca();
            error('the requested colors are out of gamut');
        end
        x = linspace(gausslocut/1000,gausshicut/1000,NGAMMASTEPS)';
        invnorm = bsxfun(@plus, bsxfun(@times, norminv(x), gl.stim.sigma), gl.stim.mu);
        gl.gaussgamma = bsxfun(@(x,y) gl.cal.invGamma(x,y), round(NGAMMASTEPS*invnorm), 1:3);
    elseif strcmp(gl.stim.type,'LMS') % Binary cone noise
        switch noisetype
            case 3 % L-M binary cone noise
                colordirlms = [1 -1 0; -1 1 0];
            otherwise % Preparing the 8x3 matrix of rgb values
                colordirlms = sign(fullfact([2,2,2])-1.5);
        end

        lms = bsxfun(@plus, gl.cal.M*gl.stim.mu',...
            bsxfun(@times, colordirlms, gl.stim.sigma)');
        rgbmat = gl.cal.invM*lms;

        if any(rgbmat(:) > 1) || any(rgbmat(:) < 0)
            sca();
            error('the requested colors are out of gamut');
        end

        % Converting to voltages (using bsxfun to get linear indices into the inv gamma table)
        idxs = bsxfun(@plus, round((NGAMMASTEPS-1)*rgbmat')+1, 0:size(gl.cal.invGamma,1):numel(gl.cal.invGamma)-1);
        gl.lmsbinaryrgbmat = gl.cal.invGamma(idxs);
    end

    % Preparing a template filled with stixel indices for putting the random
    % numbers into the image matrix.
    stixidxs = reshape(1:gl.stim.nstixperside^2,gl.stim.nstixperside,gl.stim.nstixperside);
    pixextents = ones(2*gl.stim.npixperstix,gl.stim.npixperstix); % takes into account the aspect ratio in colour mode
    gl.stim.template = kron(stixidxs,pixextents);

    gl.stim.on = 1;
    gl.stim.framecounter = 0; % Gotta reinitalize framecounter somewhere at the trial start
end

%%
% Taking random numbers, putting them through the combined gaussgamma
% function and creating an image.  The order of elements is:
% red(1,1) red(2,1)... red(n,1)... red(n,n) green(1,1)...
function DrawStim()
    global gl

    imgsize = [2*gl.stim.nstixperside*gl.stim.npixperstix gl.stim.nstixperside*gl.stim.npixperstix 3];
    rgbs = zeros(gl.stim.nstixperside^2,3);
    if strcmp(gl.stim.type,'RGB')
        [randnums,gl.stim.seed] = getEJrandnums3(3*gl.stim.nstixperside^2,gl.stim.seed);
        for gun = 1:3
            idxs = (gun-1)*gl.stim.nstixperside^2+(1:gl.stim.nstixperside^2);
            % First n^2 randnums are for red, second n^2 randnums are for green, etc.
            rgbs(:,gun) = gl.gaussgamma(randnums(idxs)+1,gun);
        end
    elseif strcmp(gl.stim.type,'LMS')
        [randnums,gl.stim.seed] = getEJrandnums3(gl.stim.nstixperside^2,gl.stim.seed);
        randnums = mod(randnums, size(gl.lmsbinaryrgbmat,1))+1;
        rgbs = gl.lmsbinaryrgbmat(randnums,:);
    end
    img = reshape(rgbs(gl.stim.template,:), imgsize);

    tex = Screen('MakeTexture', gl.windowPtr, TranslateToColourMode(img, [], gl.ccmode), [], [], 2);
    Screen('DrawTexture', gl.windowPtr, tex, [], gl.stim.drawrect, [], 0);
    Screen('Close', tex);

    gl.fliprequest = 1;
end

%%
function ShowImage()
    global gl

    gl.stim.framecounter = 0;
    gl.synthimage.on = 1;
end

%%
function ImageOff()
    global gl

    gl.synthimage.on = 0;
end

%%
function DrawImage()
    global gl

    img = zeros([2*gl.stim.nstixperside*gl.stim.npixperstix gl.stim.nstixperside*gl.stim.npixperstix 3]);
    for gun = 1:3
        plane = squeeze(gl.synthimage.im(:,:,gun));
        img(:,:,gun) = plane(gl.stim.template);
    end

    tex = Screen('MakeTexture', gl.windowPtr, TranslateToColourMode(img, [], gl.ccmode), [], [], 2);
    Screen('DrawTexture', gl.windowPtr, tex, [], gl.stim.drawrect, [], 0);
    Screen('Close', tex);
    gl.fliprequest = 1;
end

%%
function AllOff()
    global gl

    gl.fp.on = 0;
    gl.stim.on = 0;
    gl.synthimage.on = 0;
    gl.fliprequest = 1;
end

%%
function DealWithMessage(msgSize)
    global udpCom gl
    message = pnet(udpCom.sock, 'read', msgSize, 'char');
    if strncmp(message,'return',6)
        a = dbstack;  % Check whether called from another function or from command line
        if ~strcmp(a(end).name, mfilename)
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

%%
% One day maybe put in a check for gl.stim.on so we don't get stuck here
% when we're supposed to be showing a stimulus
function imageReceiver(numElements)
    global udpCom gl

    REXBUFF = 8500;
    imagevect = nan(numElements, 1); % column vector will line up with output of hex2img
    startInd = 1;
    numReceived = 0;
    while numReceived < numElements
        msgSize = pnet(udpCom.sock, 'readpacket', REXBUFF, 'noblock');
        if msgSize
            inMessage = pnet(udpCom.sock, 'read', msgSize, 'char');
            tmp = myhex2num(inMessage);
            stopInd = startInd + length(tmp)-1;
            imagevect(startInd:stopInd) = tmp;
            numReceived = sum(isfinite(imagevect));
            startInd = numReceived + 1;
        end
    end
    % reshape the image. Assume it is NxNx3.
    % alternatively, could use gl.nstixperside
    n = sqrt(numElements / 3);
    gl.synthimage.im = reshape(imagevect, [n n 3]);
end

%%
% Get the fitted gabor parameters from REX and put them in the global gl
% structure so that they'll be there when we invoke GaborEdge.m
function gaborReceiver(numElements)
    global gl udpCom

    REXBUFF = 8500;
    paramvect = nan(numElements, 1);
    startInd = 1;
    numReceived = 0;
    while numReceived < numElements
        msgSize = pnet(udpCom.sock, 'readpacket', REXBUFF, 'noblock');
        if msgSize
            inMessage = pnet(udpCom.sock, 'read', msgSize, 'char');
            tmp = myhex2num(inMessage);
            stopInd = startInd + length(tmp)-1;
            paramvect(startInd:stopInd) = tmp;
            numReceived = sum(isfinite(paramvect));
            startInd = numReceived + 1;
        end
    end
    gl.gabor.rawrgb = paramvect(1:3);
    gl.gabor.theta = paramvect(4);
    gl.gabor.lambda = paramvect(5)*2*gl.stim.npixperstix/gl.pixperdeg;  % converting from stix to dva
    gl.gabor.phi = paramvect(6);
    gl.gabor.sigma = paramvect(7)*2*gl.stim.npixperstix/gl.pixperdeg;  % converting from stix to dva
    gl.gabor.gamma = paramvect(8);
    gl.gabor.xoffset = paramvect(9)*2*gl.stim.npixperstix/gl.pixperdeg;  % converting from stix to dva
    gl.gabor.yoffset = paramvect(10)*2*gl.stim.npixperstix/gl.pixperdeg;  % converting from stix to dva
end
