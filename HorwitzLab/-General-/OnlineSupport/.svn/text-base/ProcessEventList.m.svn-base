%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Go through the eventList returned by PL_GetTS and put the information 
% in the 'p' structure.
%
% Fields of the p structure are:
%   p.times             Vector containng the times of the events
%   p.events         	Vector containing the event codes
%   p.spikes            Cell array containing the spike times
%   p.processnowflag    1 = the p structure contains a full trial of events
%                       0 = does not yet contain a full trial of events
function p = ProcessEventList(p, eventList)

    CODETYPECOL = 1;
    CHANNELCOL = 2;
    EVENTCOL = 3;
    TIMECOL = 4;
    
    SPIKE = 1;
    EVENT = 4;
    
    eventList(:,TIMECOL) = eventList(:,TIMECOL)*1000;   % converting to msec
    L = eventList(:,TIMECOL) > p.lastprocessed_t;
    Levent = eventList(:,CODETYPECOL) == EVENT;
    Lspike = eventList(:,CODETYPECOL) == SPIKE;
    
    p.times = [p.times; eventList(L&Levent,TIMECOL)];
    p.events = [p.events; eventList(L&Levent,EVENTCOL)];
    for i = 1:length(p.spikes)
        Liso = eventList(:,EVENTCOL) == i;
        spiketimevect = eventList(L&Lspike&Liso,TIMECOL);
        p.spikes{i} = [p.spikes{i}; spiketimevect];
    end
    if (any(p.events == p.processnowcode))
       p.processnowflag = 1;
    end
end