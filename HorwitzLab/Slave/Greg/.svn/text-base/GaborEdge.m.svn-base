function GaborEdge
% GaborEdge.m
%
%    Slave program for displaying Gabors with superimposed edges of various
%    variety.
%
% GDLH 8/8/08

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
    % The stuff below we can plan on inheriting
    % but commenting it in for debugging purposes
    %gl.mondistcm = 0;       % We can rely on inheriting this stuff
   % gl.screenWidthcm = 0;
   % gl.screenWidthpix = 0;
   % gl.screenHeightpix = 0;
   % gl.screenCenterXpix = 0;
   % gl.screenCenterYpix = 0;
   % gl.pixperdeg = 0;
   % gl.bkgndrgb = [0 0 0];
   % gl.windowPtr = 0;
   % gl.cal.gammaTable = []; 
   % gl.cal.monSpd = [];
   % gl.cal.fundamentals = [];
   % gl.cal.M = zeros(3);
   % gl.cal.invM = zeros(3);
   % gl.stim.x = -5;  % Need this around so that REX can query it's value.
   % gl.stim.y = 0;  % Don't zero it out.
   % gl.stim.nstixperside = 10;
   % gl.stim.npixperstix = 3;
   % gl.gabor.theta = 1;
   % gl.gabor.lambda = .4;
   % gl.gabor.phi = 0;
   % gl.gabor.sigma = .4;
   % gl.gabor.gamma = 1;
   % gl.gabor.xoffset = 0;
   % gl.gabor.yoffset = 0;
   % gl.gabor.rawrgb = [-.3 .3 -.3]';
    
    % End of stuff that we can rely on inheriting
    
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
    
    gl.ep.x = 0;
    gl.ep.y = 0;

    gl.stim.on = 0;
    gl.stim.mat = zeros(1024/2, 768/2, 3);  % One full quadrant of the screen

    % These are computed and sent back to REX
    gl.gabor.cc = [0 0 0];  % cone contrast
    gl.gabor.rgb = [0 0 0];  % delta from background
    gl.edge.cc = [0 0 0];   % cone contrast
    gl.edge.rgb = [0 0 0];  % delta from background

    [udpCom.sock, Success] = pnetStart(udpCom.port);
    if ~(Success)
        return
    end
    pnet(udpCom.sock, 'setreadtimeout', 0);
    pnet(udpCom.sock, 'setwritetimeout', 0);
    disp('In GaborEdge');

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
% This function barely does anything because the assumption is that
% everything has already been set up by WhiteNoise.m
% Gabor parameters are manipulated so that they are in degrees of visual angle 
% when REX requests them.
function InitDisplay()
    global gl;
    % InitDisplay()
    % Debugging - everything from here to the comment below
    % starting with "debugging" can be removed.
   % load(calfilename);
   % cal = cals{end};
   % gl.mondistcm = mondist;
   % gl.screenWidthcm = screenwidth;
   % gl.bkgndRGB = round(255*cal.bgColor)';
   % gl.cal.gammaTable = cal.gammaTable;
   % gl.cal.monSpd = cal.P_device;

%    s = load(fundfilename);
%    fns = fieldnames(s);
%    P_device = SplineSpd(SToWls(cal.S_device), cal.P_device, SToWls([380 5 81]));
%    gl.cal.fundamentals = eval(['s.',fns{1}])';
%    gl.cal.M = gl.cal.fundamentals'*P_device;
%    gl.cal.invM = inv(gl.cal.M);
    
    % Gamma correction is done in software so just set hardware gamma
    % lookup to the unity line.
%    clut = repmat(linspace(0,1,256),3,1)';
%    Screen('LoadNormalizedGammaTable', gl.windowPtr, clut);
        
%    if (isempty(Screen('windows')))
%        gl.windowPtr = Screen('OpenWindow',0, 255*gl.cal.bgColor);
%    else
%        gl.windowPtr = Screen('windows');
%        Screen('FillRect',gl.windowPtr,255*gl.cal.bgColor);
%    end
%    gl.framerate = Screen('FrameRate', gl.windowPtr, 1);
%    [screenwidthpix, screenheightpix]  = Screen('WindowSize', gl.windowPtr);
%    gl.screenWidthpix = screenwidthpix;
%    gl.screenHeightpix = screenheightpix;
%    gl.screenCenterXpix = screenwidthpix/2;
%    gl.screenCenterYpix = screenheightpix/2;
    
%    pixpercm = gl.screenWidthpix/gl.screenWidthcm;
%    theta = atan2(gl.screenWidthcm/2, gl.mondistcm)*180/pi;
%    cmperdeg = gl.screenWidthcm/(2*theta);
%    gl.pixperdeg = pixpercm*cmperdeg;
 
    % End of stuff for debugging 
    
    gl.cal.invgamma = InvertGamma(gl.cal.gammaTable, 1);
    gl.bkgndrgb = [gl.cal.gammaTable(gl.bkgndRGB(1)+1,1);...
                  gl.cal.gammaTable(gl.bkgndRGB(2)+1,2);...
                  gl.cal.gammaTable(gl.bkgndRGB(3)+1,3)];
    gl.stim.on = 0;
    gl.fp.on = 0;
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
function HideFP()
    global gl;
   
    gl.fp.on = 0;
    gl.fliprequest = 1;

end

%%
% Many of the passed parameters aren't used, but it's useful to have them
% passed because someday we might want to have REX change them from trial
% to trial.
function ShowStim(rfx, rfy, nframes, gausslim, gtheta, lambda, phi, sigma, gamma, xoff, yoff, gcolortype, gcont, grawr, grawg, grawb, etheta, edisp, econt, ecolortype)

    global gl;
    
    gl.framecountermax = nframes;

    bkgndlms = gl.cal.M*gl.bkgndrgb;
    cc = getEdgeColor(econt, ecolortype);  
    gl.edge.cc = cc;
    gl.edge.rgb = gl.cal.invM*(gl.edge.cc.*bkgndlms);

    if (gcolortype == 6)
        cc = gl.cal.invM'*[grawr; grawg; grawb];
    else
        cc = (gl.cal.M*[grawr; grawg; grawb])./bkgndlms;       
    end
    cc = cc.*(sqrt(gcont.^2)./sqrt(mean(cc.^2))); % normalizing RMS contrast
    if (gcont ~= 0)
        cc = cc*sign(gcont);
    end
    if (gcolortype > 0)
        cc = ManipulateGaborColor(cc, gcont, gcolortype);
    end
    gl.gabor.cc = cc;
    gl.gabor.rgb = gl.cal.invM*(gl.gabor.cc.*bkgndlms);
    
    rgbextents = [gl.bkgndrgb+gl.gabor.rgb+gl.edge.rgb,...
                    gl.bkgndrgb-gl.gabor.rgb-gl.edge.rgb];

    if (any(rgbextents(:) > 1 | rgbextents(:) < 0))
        gl.gabor.cc = [0 0 0]';
        gl.gabor.rgb = [0 0 0]';
        gl.edge.cc = [0 0 0]';
        gl.edge.rgb = [0 0 0]';
        
    end

    sigma = sigma*gl.pixperdeg/2;  % In doublewide pixels
    lambda = lambda*gl.pixperdeg/2;  % In doublewide pixels

    transformedlims = round(norminv([1-gausslim/1000 gausslim/1000],0,1)*abs(sigma)/min([1 gamma]));  % In doublewide pixels
    transformedlims = sort(transformedlims);  % In case sigma or gamma is negative
    stimsizeinpix = transformedlims(2)-transformedlims(1)+1;  % In doublewide pixels
    % drawrect(3)-drawrect(1) should be the x size of the image in pixels
    % and drawrect(4)-drawrect(2) should be the y size
    x = (gl.stim.x+xoff+gl.ep.x)*gl.pixperdeg;
    y = (gl.stim.y+yoff+gl.ep.y)*gl.pixperdeg;
    gl.drawrect = round([gl.screenCenterXpix+x-stimsizeinpix ...
                gl.screenCenterYpix-y-stimsizeinpix ...
                gl.screenCenterXpix+x+stimsizeinpix...
                gl.screenCenterYpix-y+stimsizeinpix]);
            % GH: use rfx and y instead of gl.stim.x and y?
            
    gl.ep.x = 0;
    gl.ep.y = 0;
    
    if(rem(gl.drawrect(1), 2)) %if the rectangle starts on an odd pixel
        gl.drawrect(1) = gl.drawrect(1) - 1;
        gl.drawrect(3) = gl.drawrect(3) - 1;
    end
    
    if (gl.drawrect(3)-gl.drawrect(1) ~= 2*stimsizeinpix || gl.drawrect(4)-gl.drawrect(2) ~= 2*stimsizeinpix)
        sca;
        keyboard;
        % Draw window is the wrong size
    end
    
    % Making the Gabor
    interval1 = linspace(transformedlims(1),transformedlims(2),2*stimsizeinpix);
    interval2 = linspace(transformedlims(1),transformedlims(2),stimsizeinpix);
    [X, Y] = meshgrid(interval2,interval1);
    xprime = X.*cos(-gtheta)+Y.*sin(-gtheta);
    yprime = -X.*sin(-gtheta)+Y.*cos(-gtheta);
    fittedgabor = exp(-(xprime.^2+gamma.^2.*yprime.^2)./(2.*sigma.^2)).*cos(2.*pi.*yprime./lambda-phi);

    % Making the luminance edge
    orientation = gtheta + etheta;
    % Cleaning up rounding errors - if we;re close to a multiple of pi/2,
    % make it pi/2
    if (mod(orientation, pi/2) < 1e-5 || mod(orientation, pi/2) > pi/2-1e-5) 
        orientation = pi/2*round(orientation/(pi/2));
    end

    edgeim = cosd(orientation*180/pi)*(Y+edisp*sigma*cosd(orientation*180/pi))+...
            sind(orientation*180/pi)*(X+edisp*sigma*sind(orientation*180/pi));
    edgeim(edgeim == 0) = eps;
    edgeim = sign(edgeim).*exp(-(xprime.^2+gamma.^2.*yprime.^2)./(2.*sigma.^2));

    im = zeros(2*stimsizeinpix, stimsizeinpix, 3);
    for plane = 1:3
        tmp = (fittedgabor.*gl.gabor.rgb(plane))+(edgeim.*gl.edge.rgb(plane))+gl.bkgndrgb(plane);
        tmp = round(tmp*size(gl.cal.invgamma,1)-1)+1;
        tmp = gl.cal.invgamma(tmp, plane);
        tmp = round(tmp*(size(gl.cal.invgamma,1)-1));
        im(:,:,plane) = reshape(tmp, 2*stimsizeinpix, stimsizeinpix);
    end
    gl.stim.mat = im;
    
    gl.framecounter = 0;
    gl.stim.on = 1;
    
    %%
    % A subfunction for manipulating the color of the gabor
    % for now just projecting it orthogonal to various mechanisms
    % defined by directions in cone contrast space 
    % 0 = No manipulation
    % 1 = Project color orthogonal to [1 1 0] (LMS)
    % 2 = Project color orthogonal to [1 1 .1] (LMS)
    % 3 = Project color orthogonal to [2 1 0] (LMS)
    % 4 = Achromatic
    
    function outcc = ManipulateGaborColor(incc, rmscontrast, arg)     
        if (~ismember(arg, [0 1 2 3 4 5 6]))
           error(['Unknown colortype: ',num2str(arg)]); 
        end
        
        if (arg == 1)
            lum4isolum = [1 1 0]; % in cone contrast units
        elseif (arg == 2)
            lum4isolum = [1 1 .1]; % in cone contrast units
        elseif (arg == 3)
            lum4isolum = [2 1 0]; % in cone contrast units
        end

        if (arg < 4)
            lum4isolum = lum4isolum./norm(lum4isolum);
            cc = (incc' - (lum4isolum*incc)*lum4isolum)';
        elseif (arg == 4)  % Achromatic
            if (incc(1)+incc(2) > 0) % If L+M > 0 make it bright
                cc = [1 1 1]';
            else
                cc = [-1 -1 -1]';
            end
        elseif (arg == 5) % L-M
            if (incc(1)-incc(2) > 0)
                cc = [1 -1 0]';
            else
                cc = [-1 1 0]';
            end    
        end
        if (norm(cc) > 0)
            outcc = cc.*(sqrt(rmscontrast.^2)./sqrt(mean(cc.^2)));
        else
            outcc = [0;0;0];
        end
    end
    
    %%
    % A subfunction for getting the cone contrast of the edge
    % We are using different definitions for the contrast of the
    % edge and the gabor:
    % gabor_contrast = sqrt(mean(cc.^2))
    % edge_contrast = sqrt(mean(cc([1 2]).^2)
    % This is a little awkward, but the rationale is that rms is a natural 
    % definition of contrast for the Gabor, and we're using almost the same
    % definition for the edge but ignoring the S-cone component (which
    % isn't supposed to contribute to luminance).  This definintion agrees
    % with the standard definition definition of luminance contrast 
    % as long as L and M cone contrasts are the same.
    function out = getEdgeColor(contrast, arg)   
        if (~ismember(arg, [0 1 2 3]))
           error(['Unknown colortype: ',num2str(arg)]); 
        end

        if (arg == 0)
            cc = [1 1 1]'; 
        elseif(arg == 1)
            cc = [1 1 0]'; 
        elseif(arg == 2)
            cc = [1 1 .1]';
        elseif(arg == 3)
            cc = [1 .5 0]'; % "contrast" is not strictly accurate here
        end
        out = cc.*(sqrt(contrast.^2)./sqrt(mean(cc([1 2]).^2))); % normalizing RMS contrast
        out = out.*sign(contrast);
    end
 end
%%
function DrawStim()

    global gl;
  
    tex=Screen('MakeTexture', gl.windowPtr, TranslateToColourMode(gl.stim.mat,1));
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
function eyepos(x, y)
    global gl;

    gl.ep.x = x/40;
    gl.ep.y = y/40;
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


