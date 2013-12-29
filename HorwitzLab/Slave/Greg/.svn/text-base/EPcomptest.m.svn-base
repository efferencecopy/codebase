function EPcomptest

% Experimental code for testing whther we can query
% REX via UDP for the current eye postiion on every screen refresh

   % Communication variables
    global udpCom;
    udpCom.port = 6665;
    udpCom.sock = nan;
    udpCom.rexip = '192.168.1.120';
    
    global gl;
    gl.mondist = 0;
    gl.screenwidth = 0;
    gl.pixperdeg = 0;
    gl.bkgnd.rgb = [0 0 0];
    gl.windowPtr = 0;
    gl.fliprequest = 0;
    
    gl.stim.on = 0;
    gl.stim.x = 0;
    gl.stim.y = 0;
    gl.stim.rgb = [0 0 0];
    gl.stim.size = 1;
    gl.stimvect = zeros(1000,1);
    gl.stimvectidx = 1;


    [udpCom.sock, Success] = pnetStart(udpCom.port);
    if ~(Success)
        return
    end
    pnet(udpCom.sock, 'setreadtimeout', 0);
    pnet(udpCom.sock, 'setwritetimeout', 0);

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
            if (gl.fliprequest)           
                Screen('Flip',gl.windowPtr,0,0);
                gl.fliprequest = 0;
            end
        end
        [keyisdown,secs,keycode] = KbCheck();
        if (keyisdown && keycode(41))
            ShowCursor;
            pnet(udpCom.sock, 'close');
            Screen('CloseAll');
            return;
        end
    end
end



function InitDisplay(mondist, screenwidth)
    global gl;
    
    gl.mondistcm = mondist;
    gl.screenWidthcm = screenwidth;
    gl.stim.on = 0;

         
    gl.windowPtr = Screen('OpenWindow',0, [190 190 190]);
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
function DrawStim()
    global gl;
    global udpCom

    tic
    outmsg = 'Eye pos request';
    pnet(udpCom.sock, 'write', outmsg);
    pnet(udpCom.sock, 'writepacket', '10.208.108.60', 6665);
 
    % Wait right here for the message to go through
    
    messageIsAvailable = 0;
    while ~messageIsAvailable  % endless loop.
        messageIsAvailable = pnet(udpCom.sock, 'readpacket', 1000, 'no block');
        if (messageIsAvailable)
            DealWithMessage(messageIsAvailable);
         end     
    end
    gl.stimvect(gl.stimvectidx) = toc;
    if (gl.stimvectidx < length(gl.stimvect))
        gl.stimvectidx = gl.stimvectidx + 1; 
    end
    w = gl.windowPtr;
    
    stimsizeinpix = round(gl.pixperdeg*gl.stim.size);
    [ScreenWidth, ScreenHeight] = Screen('WindowSize',w);
    screenCenterX = ScreenWidth/2;
    screenCenterY = ScreenHeight/2;
    drawrect = [screenCenterX+(gl.stim.x*gl.pixperdeg)-floor(stimsizeinpix/2)...
                screenCenterY-(gl.stim.y*gl.pixperdeg)-floor(stimsizeinpix/2)...
                screenCenterX+(gl.stim.x*gl.pixperdeg)+ceil(stimsizeinpix/2)...
                screenCenterY-(gl.stim.y*gl.pixperdeg)+ceil(stimsizeinpix/2)];
    Screen('FillOval', gl.windowPtr, gl.stim.rgb, drawrect);

    gl.fliprequest = 1;
    

end

function ShowStim()
    global gl;

    gl.stim.on = 1;
end

%%
function DealWithMessage(msgSize)
    global udpCom;
    global gl;  % needs to be here for functions evaled by this one
    message = pnet(udpCom.sock, 'read', msgSize, 'char');
  
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
function eyepos(x, y)
    global gl;

    gl.stim.x = x/40;
    gl.stim.y = y/40;
end