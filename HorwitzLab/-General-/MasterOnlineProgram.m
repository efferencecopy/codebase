function MasterOnlineProgram
    % Based on "MasterSlaveProgram" on the Mac.  GDLH 7/15/08
    
    global udpCom
    udpCom.PORT = 6665;
    udpCom.SOCK = nan;
    udpCom.REXIP = '192.168.1.120';
    
    global gl %#ok<NUSED>
    
    KbName('UnifyKeyNames');
    
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
                disp('Back in MasterOnlineProgram');
                messageIsAvailable = 0;
                break
            end
            [keyisdown,secs,keycode] = KbCheck();
            if (keyisdown && keycode(KbName('escape'))) %if you press escape
                return
            end
        end
    end
    
    %%
    function DealWithMessage(msgSize)
        
        message = pnet(udpCom.SOCK, 'read', msgSize, 'char');
        disp(['Invoking new online analysis program: ',message]);
        try
            eval(message);
        catch ME
            fprintf('Ignoring uninterpretable message: %s\n', message);
            disp(getReport(ME));
        end
    end
end
