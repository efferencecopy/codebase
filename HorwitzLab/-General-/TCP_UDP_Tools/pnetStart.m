function [sock Success] = pnetStart(port, hostIP, ShakeHands, TimeOut)

%    INITIALIZE UDP COMMUNICATION BETWEEN COMPUTERS
%
%   example: [sock, success] = pnetStart(port, [hostIP], [ShakeHands], [TimeOut])
%
% Initializes communication via UDP. Port: should be the local host port
% number. hostIP: should be represented as a STRING. If host IP is given
% then it is assumed that the udp socket should be "connected", otherwise
% it is assumed that the socket is "connectionless". If a hostIP is
% specified than the host machine must be using the same port as the local
% host. ShakeHands: if equal to 1 then the routine ensures that the host
% computer is able to send and recieve messages to/from the host computer.
% TimeOut is an optional argument and specifies how long (in Seconds) the
% funtion will attempt to initiate a UDP connection. If not specified, the
% funtion will try for 10sec. Outputs: Sock is the scalar socket
% identifier. Success returns 1 if the initialization was successful and 0
% otherwise.
%
%10/07 CAH

%%% CONSTANTS and setup %%%
DEFAULT_TIMEOUT = 10;   %the default timeout time in sec
Success = 0;


if nargin<1;
    error('Too few inputs, try again');
elseif nargin<=2;
    TimeOut = DEFAULT_TIMEOUT;
    ShakeHands = 0;
elseif nargin<=3;
    TimeOut = DEFAULT_TIMEOUT;
end


%Unix will be unhappy if you try to open an already functional socket, so
%begin by closing all the sockets and the opening a new one.
pnet('closeall');
sock = pnet('udpsocket', port);

%if the hostIP is specified then assume that the socket will be "connected".
if (nargin>=2)
    pnet(sock, 'udpconnect', hostIP, port);
    [boundIP, boundPort] = pnet(sock, 'getHost');
end

% if you don't want handshaking, stop here
if ~(ShakeHands)
    try
        status = pnet(sock, 'status');
        if(status)
            Success = 1; %(boundPort == port);
        end
    catch
        fprintf(' ** udp initialization failed  **\n');
    end
    return
end



%now enter the while loop that will run until TimeOut, Success, or
%indefinately if TimeOut is not specified. Remeber that pnet will bonk if
%it tries to send packets to an "unbound" port, so make sure that a
%hostport is bound. Then try sending and recieving packets. Based on what
%is received, update Success.
StartTime = GetSecs;
ElapsedTime = 0;
Expired = (ElapsedTime>TimeOut);


readySig = 1; %arbitrary signal
localReceipt = 0; %if this machine has successfully received a msg
hostReceipt = 0; %if the host has successfully received a msg
while ~(Expired || Success);
    if (boundPort == port);
        %try sending a message
        outMSG = [readySig, localReceipt, hostReceipt];
        pnet(sock, 'write', outMSG);
        pnet(sock, 'writepacket', hostIP, port);

        %try receiving a message, update Success Accordingly.
        
        msgSize = pnet(sock, 'readpacket', 30, 'noblock');
        
        if (msgSize>0)
            msgIN = pnet(sock, 'read', 30, 'double');
            localReceipt = msgIN(1); %did the local machine recieve a msg from a ready host?
            hostReceipt = msgIN(2);  %has the host machine recieved a msg?
            Success = isempty(find(msgIN==0));
        end

    else
        [boundIP, boundPort] = pnet(sock,'gethost');
    end

    ElapsedTime = GetSecs-StartTime;
    Expired = (ElapsedTime>TimeOut);
end


if ~(Success)
    pnet(sock, 'close');
    fprintf(' ** udp initialization failed  **\n');
end






