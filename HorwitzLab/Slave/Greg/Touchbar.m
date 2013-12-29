function Touchbar
%
%
%This is a script that allows the TouchBar Rex program to present a
%fixation point. This program was adapted from RFmapper (GH).
%
%CAH 11/07


    % Communication variables
    global udpCom;
    udpCom.PORT = 6665;
    udpCom.SOCK = nan;
    udpCom.REXIP = '192.168.1.120';

    % Variables that are destructively modified by subfunctions
    global gl;
    gl.mondist = 0;
    gl.screenwidth = 0;
    gl.pixperdeg = 0;
    gl.bkgnd.rgb = [0 0 0];
    gl.windowPtr = 0;  
    gl.fp.on = 0;
    gl.fp.RGB = [0 0 0];
    gl.fp.x = 0;
    gl.fp.y = 0;
    gl.fp.size = 0;
    
    gl.square.on = 0;
    gl.square.RGB = [0 0 0];
    gl.square.x = 0;
    gl.square.y = 0;
    gl.square.size = 0;
    
    %stuff to send to rex
    global toRex;
    toRex.calStructPlease = 'calStructPlease';
    toRex.gammaTable = nan(256,3);
    toRex.monSpect = nan(101, 3);
    toRex.M_mtx = nan(3,3);
    
    %load the calibration structure and set up the calibration data to send
    %to rex.
    load('/Monitor Calibration/Monitor data/Dell 4/Dell4BitsCal');
    load('T_cones_smj.mat');
    calStruct = cals{end};
    
    fundamentals = T_cones_smj;
    MonSpd = SplineSpd(SToWls(calStruct.S_device), calStruct.P_device, SToWls(S_cones_smj));
    toRex.M_mtx = fundamentals * MonSpd;
    toRex.gammaTable = calStruct.gammaTable;
    toRex.monSpect = calStruct.P_device;
    
    %establish a connectionless UDP socket. Do this without shaking hands
    %with the remote host(s)
    [udpCom.SOCK, Success] = pnetStart(udpCom.PORT);
    if ~(Success)
        return
    end
    pnet(udpCom.SOCK, 'setreadtimeout', 0);
    pnet(udpCom.SOCK, 'setwritetimeout', 0);
    
    messageIsAvailable = 0;
    
    while (1) % endless loop.

        %start by checking for a message
         while ~(messageIsAvailable)
            messageIsAvailable = pnet(udpCom.SOCK, 'readpacket', 1000, 'noblock');
            if (messageIsAvailable)
                DealWithMessage(messageIsAvailable);
                messageIsAvailable = 0;
                break
            end

            %press ESC to shutdown
            [keyisdown,secs,keycode] = KbCheck();
            if (keycode(41)) %if you press escape
                if (gl.windowPtr>0);
                    Screen('LoadNormalizedGammaTable', gl.windowPtr, gl.OldGammaTable);
                    Screen('CloseAll');
                    ShowCursor;
                end
                %pnet('closeall');
                return
            end
        end

        % update the stimulus accordingly
        if (gl.windowPtr > 0) 
            if (gl.fp.on)
                DrawFP();
            end
            if (gl.square.on)
                DrawStim();
            end
            
            Screen('Flip',gl.windowPtr);
        end
        
    end %while
    
end%function

%%
function InitDisplay(mondist, screenwidth, bkgndRGB)
    global gl;
    
    gl.mondist = mondist;
    gl.mondist = screenwidth;
    gl.bkgnd.rgb = repmat(bkgndRGB, 1,3);
    
    % don't use a fancy gamma for this paradigm, just use a linear CLUT
    clut = repmat(linspace(0,1,256),3,1)';
    gl.OldGammaTable = Screen('LoadNormalizedGammaTable', gl.windowPtr, clut);

    
    % Open a window
    gl.windowPtr = Screen('OpenWindow',0,gl.bkgnd.rgb);
   
    HideCursor;
    [winWidth, winHeight] = Screen('WindowSize', gl.windowPtr);
    pixpercm = winWidth/screenwidth;
    theta = atan2(screenwidth/2, mondist)*180/pi;
    cmperdeg = screenwidth/(2*theta);
    gl.pixperdeg = pixpercm*cmperdeg;
end

%%
function FPcontrol(OnOff, x, y, size, fpRGB)

    global gl;
    
    %Toggel the FP on and off
    gl.fp.on = OnOff;
    
    %set up the FP size and RGB
    gl.fp.x = x/10;
    gl.fp.y = y/10;
    gl.fp.size = size/10;
    gl.fp.RGB = repmat(fpRGB, 1, 3);
end

%%
function Stimcontrol(OnOff, x, y, size, R, G, B)

    global gl;
    
    %Toggle the square on and off
    gl.square.on = OnOff;
    
    %set up the FP size and RGB
    gl.square.x = x/10;
    gl.square.y = y/10;
    gl.square.size = size/10;
    gl.square.RGB = [R,G,B];
    if (~gl.square.on)
        gl.square.RGB = gl.bkgnd.rgb;
    end        
end


%%
function DrawFP()
    global gl;
    
    fpsizeinpix = round(gl.pixperdeg*gl.fp.size);
    [ScreenWidth, ScreenHeight] = Screen('WindowSize',gl.windowPtr);
    screenCenterX = ScreenWidth/2;
    screenCenterY = ScreenHeight/2;
    drawrect = [screenCenterX+(gl.fp.x*gl.pixperdeg)-floor(fpsizeinpix/2)...
                screenCenterY-(gl.fp.y*gl.pixperdeg)-floor(fpsizeinpix/2)...
                screenCenterX+(gl.fp.x*gl.pixperdeg)+ceil(fpsizeinpix/2)...
                screenCenterY-(gl.fp.y*gl.pixperdeg)+ceil(fpsizeinpix/2)];

    Screen('FillOval', gl.windowPtr, gl.fp.RGB, drawrect);

end

%%
function DrawStim()
    global gl;
    
    fpsizeinpix = round(gl.pixperdeg*gl.square.size);
    [ScreenWidth, ScreenHeight] = Screen('WindowSize',gl.windowPtr);
    screenCenterX = ScreenWidth/2;
    screenCenterY = ScreenHeight/2;
    drawrect = [screenCenterX+(gl.square.x*gl.pixperdeg)-floor(fpsizeinpix/2)...
                screenCenterY-(gl.square.y*gl.pixperdeg)-floor(fpsizeinpix/2)...
                screenCenterX+(gl.square.x*gl.pixperdeg)+ceil(fpsizeinpix/2)...
                screenCenterY-(gl.square.y*gl.pixperdeg)+ceil(fpsizeinpix/2)];

    Screen('FillRect', gl.windowPtr, gl.square.RGB, drawrect);
end

%%
function DealWithMessage(msgSize)
    global udpCom;
    global toRex;
    
    message = pnet(udpCom.SOCK, 'read', msgSize, 'char');
    try    
         eval(message);
    catch
        fprintf('Ignoring uninterpretable message: "%s"\n',message);
        error = lasterror;
        disp(error.message);
        disp(error.identifier);
        disp(error.stack);

    end
end

%%
function sendCals(udpCom, toRex, msg);

    %      SEND REX A PTB CALIBRATION STRUCTURE
    % 
    %  EXAMPLE: sendCals(udpCom, toRex.gammaTable, message);  ==> Sends the gammaTable
    %           sendCals(udpCom, toRex.calStructPlease, message); ==> initiates transfer
    %
    % This script is designed to feed Rex the pertinent components of a
    % calibration structure. These components include the GammaTable, the
    % monitor emission spectra, and the M matrix. The message is sent as a
    % hexadecimal character array with elements of the cal structure delimited
    % by spaces. This function requires udpCom to be a structure array that
    % specifies the udp socket, port and RexIP. toRex should be an element of
    % the calibration structure that the calling function (to sendCals) has
    % access to.
    %
    % CAH 11/20/07



    %******** CONSTANTS *******%
    REXBUFF = 8000;  %The size of the UDP read buffer on the Rex machine in bytes

    %**************************%

    %The first udp message will be from Rex's main state set requesting the
    %calibration structure. Send back an acknowledgment of the request. This
    %acknowledgment alerts the appropriate sub-function in Rex's udp "listening
    %state set" that it's time to acquire data. In this case, the
    %"achknowledgment" is simply the same string that Rex sent the Mac.
    pnet(udpCom.SOCK, 'write', ['MAC', msg]);
    pnet(udpCom.SOCK, 'writepacket', udpCom.REXIP, udpCom.PORT);
    pause(0.005); %pause to let Rex catch up.


    %Turn the calibration structure elements that rex requests into a long
    %string of space delimited hexadecimal numbers.
    calElement = sprintf('%bx ', reshape(toRex, 1, size(toRex,1).*size(toRex, 2)));
    fprintf('%s\n', calElement);

    %Chop up the cal struct into the appropriate number of segments, then
    %send it off. The last segment will cause and out of bounds error unless
    %length(calElement) is an integer multiple of maxMsgLength. Avoid this
    %error by catching the last chunk and including just the last few
    %elements in the snippet.
    maxMsgLength = REXBUFF-4; %buffer size minus 4 (header + \0)
    for a = 1:maxMsgLength:length(calElement);
        try
            chunk = calElement(a:(a+maxMsgLength-1));
        catch
            chunk = calElement(a:length(calElement));
        end

        pnet(udpCom.SOCK, 'write', ['MAC',chunk])
        pnet(udpCom.SOCK, 'writepacket', udpCom.REXIP, udpCom.PORT);

        pause(0.005); %avoid Rex buffer overflow issues (necessary?)
    end

    %Tell Rex that calElement has been sent in its entirety
    pnet(udpCom.SOCK, 'write', 'MACmsgComplete');
    pnet(udpCom.SOCK, 'writepacket', udpCom.REXIP, udpCom.PORT);

end

