%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p = CleanupEvents(p)
%
% Once we've crunched the data from a single trial
% we remove it from the p structure.  We don't just reinitialize
% the p structure though because it may contain some of the ecodes 
% from the next trial.
%
% If p.processnowcode is set to something non-zero, this function will
% remove all the events that preceded the earliest instantiation of that 
% code and set p.lastprocessed_t to the time of that code.  This is useful
% for online analyses that wait until a particular code shows up and then
% does some processing on (a trial's-worth) of data.
%
% Alternatively, if p.processnowcode is set to 0, we remove all the events
% that precede p.lastprocessed_t (which, in this circumstance, is expected 
% to be set by the user somewhere in the online analysis code itself).
% This is useful for online analyses that do not wait for a particular code
% to be dropped before processing data.
function p = CleanUpEvents(p)

    if (p.processnowcode ~= 0)
        pncodeidx = find(p.events == p.processnowcode,1);
        p.lastprocessed_t = p.times(pncodeidx);
    end 
    L = p.times <= p.lastprocessed_t;
    p.times(L) = [];
    p.events(L) = [];
    for i = 1:length(p.spikes)
        spiketimevect = p.spikes{i};
        p.spikes{i} = [spiketimevect(spiketimevect > p.lastprocessed_t)];
    end
    p.processnowflag = 0;
end