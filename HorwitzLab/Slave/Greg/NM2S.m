function NM2S
% NM2S.m
%
%    Slave program for running the non-match to sample discrimination
%    paradigm.
%
% GDLH 11/2/08

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
    gl.bkgndrgb = [0 0 0];
    gl.bkgnd.tex = [];
    
    gl.bkgndgrad.tex = [];
    gl.bkgndgrad.on = 0;
    
    gl.windowPtr = 0;
    gl.cal.gammaTable = []; 
    gl.cal.monSpd = [];
    gl.cal.fundamentals = [];
    gl.cal.M = zeros(3);
    gl.cal.invM = zeros(3);
    gl.cal.invgamma = [];
    
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
    
    gl.sactarg1.on = 0;
    gl.sactarg1.x = 0;
    gl.sactarg1.y = 0;
    gl.sactarg1.size = 0;
    gl.sactarg1.rgb = [0 0 0];
    gl.sactarg1.drawrect = [0 0 0 0];
    
    gl.sactarg2.on = 0;
    gl.sactarg2.x = 0;
    gl.sactarg2.y = 0;
    gl.sactarg2.size = 0;
    gl.sactarg2.rgb = [0 0 0];
    gl.sactarg2.drawrect = [0 0 0 0];

    gl.fp.on = 0;
    gl.fp.x = 0;
    gl.fp.y = 0;
    gl.fp.size = 0;
    gl.fp.rgb = [0 0 0];
    gl.fp.drawrect = [0 0 0 0];
    
    gl.sample.on = 0;
    gl.sample.x = 0;
    gl.sample.y = 0;
    gl.sample.drawrect = [0 0 0 0];
    gl.sample.cc = [0 0 0];
    gl.sample.rgb = [0 0 0];
    gl.sample.tex = [];
    
    gl.match.on = 0;
    gl.match.theta = 0;
    gl.match.drawrect = [0 0 0 0];
    
    gl.nonmatch.on = 0;
    gl.nonmatch.theta = 0;
    gl.nonmatch.drawrect = [0 0 0 0];
    gl.nonmatch.tex = [];

    gl.sample.cc = [0 0 0];
    gl.sample.rgb = [0 0 0];
    
    gl.ep.x = 0;
    gl.ep.y = 0;

    gl.targs.on = 0;
 
    [udpCom.sock, Success] = pnetStart(udpCom.port);
    if ~(Success)
        return
    end
    pnet(udpCom.sock, 'setreadtimeout', 0);
    pnet(udpCom.sock, 'setwritetimeout', 0);
    disp('In NM2S');

    messageIsAvailable = 0;
    while ~messageIsAvailable  % endless loop.
        messageIsAvailable = pnet(udpCom.sock, 'readpacket', 1000, 'no block');
        if (messageIsAvailable)
            DealWithMessage(messageIsAvailable);
            messageIsAvailable = 0;   
        end     
        if (gl.windowPtr > 0)
            if (gl.bkgndgrad.on)
                DrawBkgndGrad();
            end
            if (gl.sample.on)
                DrawSample();
            end
            if (gl.match.on)
                DrawMatch();
            end
            if (gl.nonmatch.on)
                DrawNonmatch();
            end
            if (gl.sactarg1.on)
                DrawSacTarg1();
            end
            if (gl.sactarg2.on)
                DrawSacTarg2();
            end
            if (gl.fp.on)
                DrawFP();
            end
            if (gl.fliprequest)
                DoFlip();
            end
        end
        if (gl.timetoleave)
            WrapItUp;
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
    if (gl.sample.on || gl.nonmatch.on)
        if(gl.framecounter == 0)
           pnet(udpCom.sock, 'write', 'MACSTIMON>> >>');
           pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
        end
        if (gl.framecounter == gl.framecountermax-1)
            pnet(udpCom.sock, 'write', 'MACSTIMOFF>> >>');
            pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
            gl.match.on = 0;
            gl.nonmatch.on = 0;
            gl.bkgndgrad.on = 0;
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
    gl.cal.monSpd = cal.P_device;

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
        gl.windowPtr = Screen('OpenWindow',0, gl.bkgndRGB);
    else
        gl.windowPtr = Screen('windows');
        gl.windowPtr = gl.windowPtr(1);
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

    gl.cal.invgamma = InvertGamma(gl.cal.gammaTable, 1);
    gl.bkgndrgb = [gl.cal.gammaTable(gl.bkgndRGB(1)+1,1);...
                  gl.cal.gammaTable(gl.bkgndRGB(2)+1,2);...
                  gl.cal.gammaTable(gl.bkgndRGB(3)+1,3)];  
              
    % Making a texture for the uniform background
    im = zeros(screenheightpix, screenwidthpix/2, 3);
    for plane = 1:3
        tmp = gl.bkgndrgb(plane);
        tmp = round(tmp*size(gl.cal.invgamma,1)-1)+1;
        tmp = gl.cal.invgamma(tmp, plane);
        tmp = round(tmp*(size(gl.cal.invgamma,1)-1));
        im(:,:,plane) = repmat(tmp, screenheightpix, screenwidthpix/2);
    end
    gl.bkgnd.tex=Screen('MakeTexture', gl.windowPtr, TranslateToColourMode(im,1));
    Screen('DrawTexture',gl.windowPtr,gl.bkgnd.tex,[],[],[],0);
    Screen('Flip',gl.windowPtr,0,0);
  
    gl.sample.on = 0;
    gl.match.on = 0;
    gl.nonmatch.on = 0;
    gl.fp.on = 0;
    gl.bkgndgrad.on = 0;
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

    gl.fp.drawrect = [gl.screenCenterXpix+(gl.fp.x*gl.pixperdeg)-floor(fpsizeinpix/2)...
                gl.screenCenterYpix-(gl.fp.y*gl.pixperdeg)-floor(fpsizeinpix/2)...
                gl.screenCenterXpix+(gl.fp.x*gl.pixperdeg)+ceil(fpsizeinpix/2)-1,...
                gl.screenCenterYpix-(gl.fp.y*gl.pixperdeg)+ceil(fpsizeinpix/2)-1];
    
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
function ShowSample(x, y, stimsize, nframes, l, m, s, theta)

    global gl;
    
    gl.sample.x = x;
    gl.sample.y = y;
    gl.sample.size = stimsize/10;
    gl.framecountermax = nframes;

    if (~isempty(gl.sample.tex))
        Screen('Close',gl.sample.tex);
        gl.sample.tex = [];
    end
    
    bkgndlms = gl.cal.M*gl.bkgndrgb;
    gl.sample.cc = [l; m; s];
    gl.sample.rgb = gl.cal.invM*(gl.sample.cc.*bkgndlms)+gl.bkgndrgb;
  
    if (any(gl.sample.rgb(:) > 1 | gl.sample.rgb(:) < 0))
        gl.sample.cc = [0 0 0]';
        gl.sample.rgb = [0 0 0]';
    end

  
    % drawrect(3)-drawrect(1) should be the x size of the image in pixels
    % and drawrect(4)-drawrect(2) should be the y size
    stimsizeinpix = round(gl.pixperdeg*gl.sample.size/2);
  %  drawrect = round([gl.screenCenterXpix+(gl.sample.x*gl.pixperdeg)-stimsizeinpix ...
  %             gl.screenCenterYpix-(gl.sample.y*gl.pixperdeg)-stimsizeinpix ...
  %              gl.screenCenterXpix+(gl.sample.x*gl.pixperdeg)+stimsizeinpix...
  %              gl.screenCenterYpix-(gl.sample.y*gl.pixperdeg)+stimsizeinpix]);
    drawrect = round([gl.screenCenterXpix+(gl.sample.x*gl.pixperdeg)-stimsizeinpix*sqrt(2)...
                gl.screenCenterYpix-(gl.sample.y*gl.pixperdeg)-stimsizeinpix*sqrt(2)...
                gl.screenCenterXpix+(gl.sample.x*gl.pixperdeg)+stimsizeinpix*sqrt(2)...
                gl.screenCenterYpix-(gl.sample.y*gl.pixperdeg)+stimsizeinpix*sqrt(2)]);
            
    if(rem(drawrect(1), 2)) %if the rectangle starts on an odd pixel
        drawrect(1) = drawrect(1) - 1;
        drawrect(3) = drawrect(3) - 1;
    end
    
   % if (drawrect(3)-drawrect(1) ~= 2*stimsizeinpix || drawrect(4)-drawrect(2) ~= 2*stimsizeinpix)
   %     sca;
   %     keyboard;
        % Draw window is the wrong size
   % end
    gl.sample.drawrect = drawrect;

    % Texture for the sample
    im = zeros(round(2*stimsizeinpix*sqrt(2)), round(stimsizeinpix*sqrt(2)), 3);
    tmprows = linspace(-sqrt(2)/2, sqrt(2)/2, round(2*stimsizeinpix*sqrt(2)));
    tmpcols = linspace(-sqrt(2)/2, sqrt(2)/2, round(stimsizeinpix*sqrt(2)));
    [x,y] = meshgrid(tmpcols, tmprows);
    rotx = cos(-theta)*x+sin(-theta)*y;
    roty = -sin(-theta)*x+cos(-theta)*y;
    L = rotx < 0.5 & rotx > -0.5 & roty < 0.5 & roty > -0.5;
    tmpplane = zeros(length(tmprows), length(tmpcols));

    for plane = 1:3
        tmp1 = gl.sample.rgb(plane);
        tmp1 = round(tmp1*size(gl.cal.invgamma,1)-1)+1;
        tmp1 = gl.cal.invgamma(tmp1, plane);
        tmp1 = round(tmp1*(size(gl.cal.invgamma,1)-1));
        tmp0 = gl.bkgndrgb(plane);
        tmp0 = round(tmp0*size(gl.cal.invgamma,1)-1)+1;
        tmp0 = gl.cal.invgamma(tmp0, plane);
        tmp0 = round(tmp0*(size(gl.cal.invgamma,1)-1));
        tmpplane(L) = tmp1;
        tmpplane(~L) = tmp0;
        im(:,:,plane) = tmpplane;
        %im(:,:,plane) = repmat(tmp, 2*stimsizeinpix, stimsizeinpix);
    end
    gl.sample.tex=Screen('MakeTexture', gl.windowPtr, TranslateToColourMode(im,1));
 
    gl.framecounter = 0;
    gl.sample.on = 1;
 end
%%
function DrawSample()
    global gl;
  
    Screen('DrawTexture',gl.windowPtr,gl.sample.tex, [], gl.sample.drawrect,[],0);
    gl.fliprequest = 1;
end

%%
function HideSample()
    global gl;

    gl.sample.on = 0;
    gl.fliprequest = 1;
     
end
%%
% All cone contrasts (l, m, s, gradl, gradm, grads) are with respect to the
% standard background (on which the calibration was done) given by
% gl.bkgndrgb.  This will hopefully make the comparisons easier to make.
% The luminance pedestal is included in l, m, s, which specify the
% correct (nonmatch) stimulus.
%
% To ensure a local match (i.e the local cone contrast across the edge of
% the two stimuli is matched) the following condition must be satisfied:
%
% (nonmatchcc+1)/(gradcc+1) = samplecc+1
%
% In other words, if the nonmatch S-cone contrast were 0.26 and the same
% S-cone contrast were 0.20, the gradient S-cone contrast would have to be
% 0.05 to produce a local match.
function ShowDiscrim(matchx, matchy, nonmatchx, nonmatchy, theta, stimsize, nframes, l, m, s, gradl, gradm, grads, gradpos)

    global gl;
    
    gl.match.x = matchx;
    gl.match.y = matchy;
    gl.match.theta = theta;
    gl.nonmatch.x = nonmatchx;
    gl.nonmatch.y = nonmatchy;
    gl.nonmatch.theta = theta;
    gl.framecountermax = nframes;
    if (matchx == 0 && matchy == 0)
        return;
    end

    bkgndlms = gl.cal.M*gl.bkgndrgb;
    gl.nonmatch.cc = [l; m; s];
    gl.nonmatch.rgb = gl.cal.invM*(gl.nonmatch.cc.*bkgndlms)+gl.bkgndrgb;
  
    if (any(gl.nonmatch.rgb(:) > 1 | gl.nonmatch.rgb(:) < 0))
        gl.nonmatch.cc = [0 0 0]';
        gl.nonmatch.rgb = [0 0 0]';
    end

    if (~isempty(gl.nonmatch.tex))
        Screen('Close',gl.nonmatch.tex);
        gl.nonmatch.tex = [];
    end
    if (~isempty(gl.bkgndgrad.tex))
        Screen('Close',gl.bkgndgrad.tex);
        gl.bkgndgrad.tex = [];
    end

    % Making background screen with the gradient
    gradendpointsinpix = round([gl.screenCenterXpix-(abs(gradpos)*gl.pixperdeg),...
                            gl.screenCenterXpix+(abs(gradpos)*gl.pixperdeg)]./2);
    
    gradccs = [0 0 0; gradl, gradm, grads]';  % assuming nonmatch is on the right
    if (gradpos < 0); % nonmatch on left
        gradccs = fliplr(gradccs);
    end
    for cone = 1:3
        tmp = [repmat(gradccs(cone,1),1,gradendpointsinpix(1)),...
                      linspace(gradccs(cone,1),gradccs(cone,2),gradendpointsinpix(2)-gradendpointsinpix(1)),...
                      repmat(gradccs(cone,2),1,gl.screenWidthpix/2-gradendpointsinpix(2))];
        gradientlms(cone,:) = bkgndlms(cone)*(1+tmp);
    end

    gradientrgb = gl.cal.invM*gradientlms;
    im = zeros(gl.screenHeightpix, gl.screenWidthpix/2, 3);
    for plane = 1:3
        tmp = gradientrgb(plane,:);
        tmp = round(tmp*size(gl.cal.invgamma,1)-1)+1;
        tmp = gl.cal.invgamma(tmp, plane);
        tmp = round(tmp*(size(gl.cal.invgamma,1)-1));
        im(:,:,plane) = repmat(tmp', gl.screenHeightpix, 1);
    end
    gl.bkgndgrad.tex=Screen('MakeTexture', gl.windowPtr, TranslateToColourMode(im,1));
  
    % drawrect(3)-drawrect(1) should be the x size of the image in pixels
    % and drawrect(4)-drawrect(2) should be the y size
    stimsize = stimsize/10;  % Now in degrees
    stimsizeinpix = round(gl.pixperdeg*stimsize/2);
   
    targpos = [gl.match.x gl.match.y; gl.nonmatch.x gl.nonmatch.y];
    targpos = targpos + [gl.ep.x gl.ep.y; gl.ep.x gl.ep.y];
    gl.ep.x = 0;
    gl.ep.y = 0;
    
    for j = 1:2
       % drawrect = round([gl.screenCenterXpix+(targpos(j,1)*gl.pixperdeg)-stimsizeinpix*sqrt(2)...
       %     gl.screenCenterYpix-(targpos(j,2)*gl.pixperdeg)-stimsizeinpix*sqrt(2)...
       %     gl.screenCenterXpix+(targpos(j,1)*gl.pixperdeg)+stimsizeinpix*sqrt(2)...
       %     gl.screenCenterYpix-(targpos(j,2)*gl.pixperdeg)+stimsizeinpix*sqrt(2)]);
        drawrect = gl.sample.drawrect +...
            round([targpos(j,1) -targpos(j,2) targpos(j,1) -targpos(j,2)]*gl.pixperdeg);
         
       if(rem(drawrect(1), 2)) %if the rectangle starts on an odd pixel
            drawrect(1) = drawrect(1) - 1;
            drawrect(3) = drawrect(3) - 1;
        end
        if (j == 1)
            gl.match.drawrect = drawrect;
        else
            gl.nonmatch.drawrect = drawrect;
        end
    end

    % Texture for the distractor
    im = zeros(round(2*stimsizeinpix*sqrt(2)), round(stimsizeinpix*sqrt(2)), 3);
    tmprows = linspace(-sqrt(2)/2, sqrt(2)/2, round(2*stimsizeinpix*sqrt(2)));
    tmpcols = linspace(-sqrt(2)/2, sqrt(2)/2, round(stimsizeinpix*sqrt(2)));
    [x,y] = meshgrid(tmpcols, tmprows);
    rotx = cos(-theta)*x+sin(-theta)*y;
    roty = -sin(-theta)*x+cos(-theta)*y;
    L = rotx < 0.5 & rotx > -0.5 & roty < 0.5 & roty > -0.5;
    tmpplane = zeros(length(tmprows), length(tmpcols));

    for plane = 1:3
        tmp1 = gl.nonmatch.rgb(plane);
        tmp1 = round(tmp1*size(gl.cal.invgamma,1)-1)+1;
        tmp1 = gl.cal.invgamma(tmp1, plane);
        tmp1 = round(tmp1*(size(gl.cal.invgamma,1)-1));
        if (gradpos < 0); % nonmatch on left
            tmp0 = gradientrgb(plane,1);
        else
            tmp0 = gradientrgb(plane,end);           
        end
        tmp0 = round(tmp0*size(gl.cal.invgamma,1)-1)+1;
        tmp0 = gl.cal.invgamma(tmp0, plane);
        tmp0 = round(tmp0*(size(gl.cal.invgamma,1)-1));
        tmpplane(L) = tmp1;
        tmpplane(~L) = tmp0;
        im(:,:,plane) = tmpplane;
    end
    
    %%%%%%%
    gl.nonmatch.tex=Screen('MakeTexture', gl.windowPtr, TranslateToColourMode(im,1));

    gl.framecounter = 0;
    gl.match.on = 1;
    gl.nonmatch.on = 1;
    gl.bkgndgrad.on = 1;
end

%%
function DrawMatch()
    global gl;

%    Screen('DrawTexture',gl.windowPtr,gl.sample.tex, [], gl.match.drawrect,gl.match.theta*180/pi,0);
    Screen('DrawTexture',gl.windowPtr,gl.sample.tex, [], gl.match.drawrect,[],0);
    
    gl.fliprequest = 1;
end

%%
function DrawNonmatch()
    global gl;
    
    %Screen('DrawTexture',gl.windowPtr,gl.nonmatch.tex, [], gl.nonmatch.drawrect,gl.nonmatch.theta*180/pi,0);
    Screen('DrawTexture',gl.windowPtr,gl.nonmatch.tex, [], gl.nonmatch.drawrect, [],0);
   
    gl.fliprequest = 1;
end

%%
function DrawBkgndGrad()
    global gl;
    
    Screen('DrawTexture',gl.windowPtr,gl.bkgndgrad.tex, [], [], [],0);
    gl.fliprequest = 1;
end


%%
function HideDiscrim()
    global gl;

    gl.match.on = 0;
    gl.nonmatch.on = 0;
    gl.bkgndgrad.on = 0;
    gl.fliprequest = 1;
     
end

%%
function ShowTargs(x1, y1, x2, y2, size, r, g, b)
    global gl;

    gl.sactarg1.x = x1/10;
    gl.sactarg1.y = y1/10;
    gl.sactarg1.size = size/10;
    gl.sactarg1.rgb = [r, g, b];
    
    gl.sactarg2.x = x2/10;
    gl.sactarg2.y = y2/10;
    gl.sactarg2.size = size/10;
    gl.sactarg2.rgb = [r, g, b];
    
    sizeinpix = round(gl.pixperdeg*size/10);

    gl.sactarg1.drawrect = [gl.screenCenterXpix+(gl.sactarg1.x*gl.pixperdeg)-floor(sizeinpix/2)...
                gl.screenCenterYpix-(gl.sactarg1.y*gl.pixperdeg)-floor(sizeinpix/2)...
                gl.screenCenterXpix+(gl.sactarg1.x*gl.pixperdeg)+ceil(sizeinpix/2)-1,...
                gl.screenCenterYpix-(gl.sactarg1.y*gl.pixperdeg)+ceil(sizeinpix/2)-1];

    gl.sactarg2.drawrect = [gl.screenCenterXpix+(gl.sactarg2.x*gl.pixperdeg)-floor(sizeinpix/2)...
                gl.screenCenterYpix-(gl.sactarg2.y*gl.pixperdeg)-floor(sizeinpix/2)...
                gl.screenCenterXpix+(gl.sactarg2.x*gl.pixperdeg)+ceil(sizeinpix/2)-1,...
                gl.screenCenterYpix-(gl.sactarg2.y*gl.pixperdeg)+ceil(sizeinpix/2)-1];
 
    gl.sactarg1.on = 1;
    gl.sactarg2.on = 1;
    
end
%%
function DrawSacTarg1()
    global gl;

    Screen('Fillrect', gl.windowPtr, gl.sactarg1.rgb, gl.sactarg1.drawrect);
    gl.fliprequest = 1;
end

%%
function DrawSacTarg2()
    global gl;

    Screen('Fillrect', gl.windowPtr, gl.sactarg2.rgb, gl.sactarg2.drawrect);
    gl.fliprequest = 1;
end
%%  
function AllOff()
    global gl;
   
    gl.sample.on = 0;
    gl.match.on = 0;
    gl.nonmatch.on = 0;
    gl.fp.on = 0;
    gl.sactarg1.on = 0;
    gl.sactarg2.on = 0;
    gl.bkgndgrad.on = 0;
    
    gl.fliprequest = 1;
    
end

%%
function WrapItUp()
    global gl;
    
    gl.windowPtr
    
    if (~isempty(gl.bkgnd.tex))
        disp('Closing bkgnd tex');
        Screen('Close',gl.bkgnd.tex);
    end
    if (~isempty(gl.bkgndgrad.tex))
        disp('Closing bkgnd gradient tex');
        Screen('Close',gl.bkgndgrad.tex);
    end
    if (~isempty(gl.nonmatch.tex))
        disp('Closing nonmatch tex');
        Screen('Close',gl.nonmatch.tex);
    end
    if (~isempty(gl.sample.tex))
        disp('Closing sample tex');
        Screen('Close',gl.sample.tex);
    end
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
