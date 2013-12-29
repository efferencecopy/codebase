function p = rewonline()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Seeing if we can give a reward based on spikes.
%
%  GDLH 4/20/08
%
%  It takes 13-16 msec from the time of the spike to the
%  time of the reward (put reward TTL into analog Plexon
%  channel and compared spike time to rising phase of TTL
%  in NeuroExplorer (demo mode).
%
% GDLH 11/26/08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global Variables
global gl;  

global udpCom;  % Needed by "DealWithMessage()" and associated subfunctions
    udpCom.sock = [];
    udpCom.port = 6665;
    udpCom.rexip = '192.168.1.120';

% Some initializations
p = InitPStruct(0); % Structure of the events from each trial
s = InitPlex();     % Link to Plexon datastream
[sock, Success] = pnetStart(6665);  % UDP communication with REX
stopnow = 0;    % Flag to break out of the endless while loop (hit ESC to set to 1)
gl_tau = 100;  % time constant of 100 ms
gl_threshold = 1;
gl_spiketimebuffer = [0];
figure(1);
set(gcf,'DefaultAxesUnits','pixels')
set(gcf,'position',[150 300 900 250]);
set(gcf,'ButtonDownFcn','drawnow');  % in case user drags figure
threshslider = uicontrol('style','slider','Min',0,'Max',5,'Value',1,'Position',[25 220 200 20]);
axes('position',[ 25 25 200 100]); hold on;

% The main loop
while (~stopnow)
    stopnow = CheckForESCKey();
    
    % Check for a message from REX
    msgSize = pnet(udpCom.sock, 'readpacket', 200, 'noblock');
    if(msgSize)
        DealWithMessage(msgSize);
    end
    [n, eventList] = PL_GetTS(s);
    if (n > 0)
       p = ProcessEventList(p, eventList);
    end
    if (~any(p.spikes{1}))
        continue;
    end
    if (p.spikes{1}(end) ~= gl_spiketimebuffer(end))
        gl_spiketimebuffer = [gl_spiketimebuffer; p.spikes{1}];
        lastspiketime = gl_spiketimebuffer(end);
        timeweightedspikes = exp((gl_spiketimebuffer-lastspiketime)./gl_tau);
        spikerate = sum(timeweightedspikes(1:end-1));
        if (spikerate > get(threshslider,'Value')) % Is it time to give a reward?
            sendToRex(udpCom, [], 'integer', 'REWON')
        elseif (timeweightedspikes(1) < 10^-3)  % getting rid of really old spikes
            gl_spiketimebuffer(timeweightedspikes < 10^-3) = []; 
        end
        plot(lastspiketime,spikerate,'ko');
        p.lastprocessed_t = lastspiketime;
    end
    p = CleanUpEvents(p);
end % big while loop
PL_Close(s);
end
