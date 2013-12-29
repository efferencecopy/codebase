function MasterSlaveProgram
 
    global udpCom;
    udpCom.PORT = 6665;
    udpCom.SOCK = nan;
    udpCom.REXIP = '192.168.1.120';
    
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
                disp('Back in MasterSlaveProgram');   
                messageIsAvailable = 0;
                break
            end
            [keyisdown,secs,keycode] = KbCheck();
            if keyisdown && keycode(KbName('escape')) %if you press escape
                pnet(udpCom.SOCK, 'close');
                return
            end
        end
    end
    
    %%
    function DealWithMessage(msgSize)
  
        msg = pnet(udpCom.SOCK, 'read', msgSize, 'char');
        disp(['Invoking new slave function.  Message is: ',msg]);
        try    
            eval(msg);
        catch ME
            fprintf('Ignoring uninterpretable message: "%s"\n',msg);
            disp('in MasterSlaveProgram')
            disp(getReport(ME));
        end
    end
end
