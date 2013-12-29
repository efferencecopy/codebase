function WhiteNoise2
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
%
% Experimenting with creating non-white noise by creating 
% independent Gaussian RVs and multiplying them by the square root
% of the desired covariance matrix.
%
% GDLH 5/31/09

    % Communication variables
    global udpCom;
    udpCom.port = 6665;
    udpCom.sock = nan;
    udpCom.rexip = '192.168.1.120';
   
    % Keypress constants
    global KEY;
    KEY.ESC = 41;

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
    gl.windowPtr = 0;
    gl.fliprequest = 0;
    gl.gaussgamma = zeros(65536,3);
    gl.framerate = 0;
    gl.lmsbinaryrgbmat = [];
    gl.timetoleave = 0;
    
    gl.fp.on = 0;
    gl.fp.x = 0;
    gl.fp.y = 0;
    gl.fp.size = 0;
    gl.fp.rgb = [0 0 0];
    gl.fp.drawrect = [0 0 0 0];
    
    gl.stim.on = 0;
    gl.stim.x = 0;
    gl.stim.y = 0;
    gl.stim.seed = 1;  % A bit of hack
    gl.stim.npixperstix = 0;
    gl.stim.nstixperside = 0;
    gl.stim.type = 'RGB';
    gl.stim.mu = [0 0 0];
    gl.stim.sigma = [1 1 1];
    gl.stim.drawrect = [0 0 0 0];
    gl.stim.framecounter = 0;
    gl.stim.template = [];
    gl.stim.gausslims = [0 0];
    gl.stim.absgausslims = [-6 6];
   
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
    if (udpCom.sock == -1)
        error('UDP not properly initialized')
    end
    pnet(udpCom.sock, 'setreadtimeout', 0);
    pnet(udpCom.sock, 'setwritetimeout', 0);

    messageIsAvailable = 0;    
    while ~messageIsAvailable  % endless loop.
        messageIsAvailable = pnet(udpCom.sock, 'readpacket', 200, 'no block');
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
            if (gl.synthimage.on)
               DrawImage();
            end
            if (gl.fliprequest)
               DoFlip();
            end
        end
        if (gl.timetoleave)
            out = gl;
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
function DoFlip()
    global gl;
    global udpCom;
    
    Screen('Flip',gl.windowPtr,0,0);
    if (gl.stim.on || gl.synthimage.on)
        if(gl.stim.framecounter == 0)
           pnet(udpCom.sock, 'write', 'MACSTIMON>> >>');
           pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
        end
        gl.stim.framecounter = gl.stim.framecounter + 1;
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
    gl.cal.invgamma = InvertGamma(gl.cal.gammaTable, 1);
    gl.cal.monSpd = cal.P_device;
    gl.stim.on = 0;
    gl.fp.on = 0;

    s = load(fundfilename);
    fns = fieldnames(s);
    P_device = SplineSpd(SToWls(cal.S_device), cal.P_device, SToWls([380 5 81]));
    gl.cal.fundamentals = eval(['s.',fns{1}])';
    gl.cal.M = gl.cal.fundamentals'*P_device;
    gl.cal.invM = inv(gl.cal.M);

    
    % Gamma correction is done in software so just set hardware gamma
    % lookup to the unity line.
    clut = repmat(linspace(0,1,256),3,1)';
    Screen('LoadNormalizedGammaTable', gl.windowPtr, clut);
        
    if (isempty(Screen('windows')))
        gl.windowPtr = Screen('OpenWindow',0, 255*cal.bgColor);
    else
        gl.windowPtr = Screen('windows');
        Screen('FillRect',gl.windowPtr,255*cal.bgColor);
        disp('Screen already open');
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
function ShowStim(rfx, rfy, seed, npixperstix, nstixperside, mu1, mu2, mu3, sigma1, sigma2, sigma3, gausslocut, gausshicut, noisetype, nbasisvects)
    global gl;

    NGAMMASTEPS = size(gl.cal.invgamma,1);
    
    % GDLH adding stuff to test out adaptive covariance stuff
    %gl.covmat = eye(3*nstixperside*nstixperside);
    ndims = nbasisvects;
    gl.covmat = normrnd(0,1,3*nstixperside^2,ndims);
    gl.covmat = gl.covmat*gl.covmat';
    sds = sqrt(diag(gl.covmat));
    gl.covmat = gl.covmat.*((1./sds)*(1./sds)');  % normalizing the variances    
    gl.sqrtcovmat = real(sqrtm(gl.covmat));
    % End of covariance stuff
    

   gl.stim.x = rfx/10;
   gl.stim.y = rfy/10;
   gl.bar.xy = [gl.stim.x gl.stim.y]; % Updating the RF locaiton estimate in case we adjust it manually
   if (isempty(gl.synthimage.im))
       disp('Initializing synthimage');
       gl.synthimage.im = zeros([nstixperside, nstixperside, 3]);
   end
    
    gl.stim.seed = seed;
    gl.stim.npixperstix = npixperstix;
    gl.stim.nstixperside = nstixperside;
    if (noisetype == 1)
        gl.stim.type = 'RGB';
    elseif (noisetype == 2)
        gl.stim.type = 'LMS';
    end
    stimsizeinpix = gl.stim.nstixperside*gl.stim.npixperstix; 
    gl.stim.drawrect = round([gl.screenCenterXpix+(gl.stim.x*gl.pixperdeg)-floor(stimsizeinpix) ...
                gl.screenCenterYpix-(gl.stim.y*gl.pixperdeg)-floor(stimsizeinpix)-1 ...
                gl.screenCenterXpix+(gl.stim.x*gl.pixperdeg)+ceil(stimsizeinpix) ...
                gl.screenCenterYpix-(gl.stim.y*gl.pixperdeg)+ceil(stimsizeinpix)]);

    if(rem(gl.stim.drawrect(1), 2)) %if the rectangle starts on an odd pixel
        gl.stim.drawrect(1) = gl.stim.drawrect(1) - 1;
        gl.stim.drawrect(3) = gl.stim.drawrect(3) - 1;
    end
    gl.stim.mu = [gl.cal.gammaTable(gl.bkgndRGB(1)+1,1),...
                  gl.cal.gammaTable(gl.bkgndRGB(2)+1,2),...
                  gl.cal.gammaTable(gl.bkgndRGB(3)+1,3)]+[mu1 mu2 mu3]/1000;
    gl.stim.sigma = [sigma1 sigma2 sigma3]/1000;
    gl.stim.gausslims = [gausslocut gausshicut]/1000;
    if (gl.stim.type == 'RGB')  % Gaussian gun noise
        % Preparing the combined Gaussian/Gamma functions
        x = linspace(gl.stim.gausslims(1),gl.stim.gausslims(2),NGAMMASTEPS);
        
        % New stuff for dealing with converting standard
        % normals to gun voltages for which the intensities are truncated
        % normals.  GDLH 6/2/09 

        % Making a truncated normal inverse CDF by using the standard normal
        % inverse CDF and then stretching it slightly in x so that
        % f(0) = minval and f(1) = maxval.
        truncnorm.inv = norminv(x);
        truncnorm.x = linspace(0,1,NGAMMASTEPS);

        % Finding a function that converts standard normals into truncated normals.
        % with the appropriate mean and variance.
        normconvert.x = linspace(gl.stim.absgausslims(1),gl.stim.absgausslims(2),NGAMMASTEPS)';  % Standard normals essentially never leave the range [-6:6]
        y = interp1(truncnorm.x, truncnorm.inv, normcdf(normconvert.x))';
        for gun = 1:3
            normconvert.y(:,gun) = y*gl.stim.sigma(gun)+gl.stim.mu(gun);
        end

        % Now integrating the inverse gamma functions with the normconvert
        % function.
        for gun = 1:3
            gl.gaussgamma(:,gun) = round(gl.cal.invgamma(round(normconvert.y(:,gun)*NGAMMASTEPS),gun)*NGAMMASTEPS);
        end
        
        if (any(gl.gaussgamma(:) < 0) | any(gl.gaussgamma(:) > NGAMMASTEPS-1))
            sca;
            error('requested colors out of gamut');
        end
        
        % GDLH 5/31/09 for the new adaptive stuff, precomputing the normal
        % inverse cdf and the inverse gamma functions separately.
        gl.invnorm = norminv(x);
        
    elseif (gl.stim.type == 'LMS')  % Binary cone noise
        % Preparing the 8x3 matrix of rgb values
        colordirlms = sign(fullfact([2,2,2])-1.5);
        rgbmat = [];
        for i = 1:size(colordirlms,1)
            lms = (gl.cal.M*gl.stim.mu')+(colordirlms(i,:).*gl.stim.sigma)';
            rgbmat = [rgbmat, gl.cal.invM*lms];
        end
        if (any(rgbmat(:) > 1 | rgbmat(:) < 0))
            sca;
            error('requested colors out of gamut');
        end
        % Converting to voltages
        for gun = 1:3
            gl.lmsbinaryrgbmat(:,gun) = gl.cal.invgamma(round(65535*rgbmat(gun,:)')+1,gun);
        end
        gl.lmsbinaryrgbmat = round(65535*gl.lmsbinaryrgbmat);
    end
    % Preparing a template filled with indices for putting the random
    % numbers into the image matrix.
    % Order is : [1 2... n; n+1 n+2... 2n; ... ]
    template = [];
    for i = 1:gl.stim.nstixperside
        tmp = [1:gl.stim.nstixperside]+(gl.stim.nstixperside*(i-1));
        tmp = repmat(tmp,2*(gl.stim.npixperstix^2),1);
        template = [template, reshape(tmp,gl.stim.npixperstix,2*gl.stim.nstixperside*gl.stim.npixperstix)'];
    end
    gl.stim.template = template;
    gl.stim.on = 1;
    gl.stim.framecounter = 0;  % Gotta reinitalize framecounter somewhere at the trial start
end

%%
% Taking random numbers, putting them through the combined gaussgamma
% function and creating an image.  The order of elements is:
% red(1,1) red(2,1)... red(n,1)... red(n,n) green(1,1)... 
function DrawStim()
    global gl;

    
    NGAMMASTEPS = size(gl.cal.invgamma,1);
    image = zeros([2*gl.stim.nstixperside*gl.stim.npixperstix gl.stim.nstixperside*gl.stim.npixperstix 3]);
        
    if (gl.stim.type == 'RGB');
        [randnums,gl.stim.seed] = getEJrandnums3(3*gl.stim.nstixperside^2,gl.stim.seed);
        % First, we need to convert these randnums into standard normal
        % variates
        randnums = norminv((randnums+1)./(NGAMMASTEPS+1)); % quantized, truncated normals [-4.1696:4.1696]
        % The above line takes a lot of time.  Could make this faster by
        % precaculating a norminv lookup table.  GDLH 6/3/09
        randnums = (gl.sqrtcovmat*randnums)';
        slope = NGAMMASTEPS/(gl.stim.absgausslims(2)-gl.stim.absgausslims(1));
        intercept = 2^15;   % 2^16/2
        plane = round(slope*randnums+intercept);        
        
        for gun = 1:3
            idxs = (gun-1)*gl.stim.nstixperside^2+[1:gl.stim.nstixperside^2];
            try
                tmp = gl.gaussgamma(plane(idxs),gun);
            catch
                sca
                keyboard
            end
            image(:,:,gun) = tmp(gl.stim.template);
        end        
    elseif (gl.stim.type == 'LMS')
        [randnums,gl.stim.seed] = getEJrandnums3(gl.stim.nstixperside^2,gl.stim.seed);
        randnums = mod(randnums, size(gl.lmsbinaryrgbmat,1))+1;
        rgbs = gl.lmsbinaryrgbmat(randnums,:);
        for gun = 1:3
            tmp = rgbs(:,gun);
            image(:,:,gun) = tmp(gl.stim.template);
        end
    end
        
    tex=Screen('MakeTexture', gl.windowPtr, TranslateToColourModeMex(image));
    Screen('DrawTexture',gl.windowPtr,tex, [], gl.stim.drawrect,[],0);
    Screen('Close',tex);

    gl.fliprequest = 1;
end

%%
function ShowImage()
    global gl;

    gl.stim.framecounter = 0;  
    gl.synthimage.on = 1;
end

%%
function ImageOff()
    global gl;

    gl.synthimage.on = 0;
end

%% 
function DrawImage()
    global gl;

    image = zeros([2*gl.stim.nstixperside*gl.stim.npixperstix gl.stim.nstixperside*gl.stim.npixperstix 3]);
    for gun = 1:3
        plane = squeeze(gl.synthimage.im(:,:,gun));
        image(:,:,gun) = plane(gl.stim.template);
    end
    tex=Screen('MakeTexture', gl.windowPtr, TranslateToColourModeMex(image));
    Screen('DrawTexture',gl.windowPtr,tex, [], gl.stim.drawrect,[],0);
    Screen('Close',tex);
    gl.fliprequest = 1;
end

%%
function AllOff()
    global gl;
   
    gl.fp.on = 0;
    gl.stim.on = 0;
    gl.synthimage.on = 0;
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

%%
% One day maybe put in a check for gl.stim.on so we don't get stuck here
% when we're supposed to be showing a stimulus
function imageReceiver(numElements)
    global udpCom;
    global gl;

    REXBUFF = 8000;
    imagevect = nan(numElements, 1); %column vector will line up with output of hex2img
    startInd = 1;
    numReceived = 0;
    while (numReceived < numElements)
        msgSize = pnet(udpCom.sock, 'readpacket', REXBUFF, 'noblock');
        if(msgSize)
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
    n = sqrt(numElements ./ 3);
    gl.synthimage.im = reshape(imagevect, [n, n, 3]);
end

%%
% Get the fitted gabor parameters from REX and put them in the global gl
% structure so that they'll be there when we invoke GaborEdge.m
function gaborReceiver(numElements)
    global udpCom;
    global gl;
    
    REXBUFF = 8000;
    paramvect = nan(numElements, 1);
    startInd = 1;
    numReceived = 0;
    while (numReceived < numElements)
        msgSize = pnet(udpCom.sock, 'readpacket', REXBUFF, 'noblock');
        if(msgSize)
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


